"""Proposed heuristic approach with the multi-MP method in Algorithm 1.

The method coordinates a central dispatch master with independent substation
topology masters and contingency feasibility/optimality subproblems.
"""


import pandas as pd
import numpy as np
import itertools
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
import time
from copy import deepcopy
from tqdm import tqdm

from .original_MIP_model import *  # Shared data preparation and AC-OPF utilities.


# Originally tested with Python 3.7 and Gurobi 10.0.3.


def _extract_linearization_points(model, Lines, cont_list):
    """Extract branch-voltage and angle linearization points from a solution."""
    Vol0_li = {}
    delta0_li = {}
    for l,i,j in Lines:
        for c in cont_list:
            suffix = f'[{l},{i},{j},{c}]'
            V2_value = model.getVarByName('V2_li' + suffix).x
            delta_value = model.getVarByName('delta_li' + suffix).x
            Vol0_li[l,i,j,c] = np.sqrt(V2_value)
            delta0_li[l,i,j,c] = delta_value
    return Vol0_li, delta0_li


def create_hmmp_line_fsp(
    data,
    cont_list=None,
    Vol0_li=None,
    delta0_li=None,
):
    """Build the non-radial line feasibility subproblem (FSP).

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        cont_list: Line-contingency identifiers represented in the model.
        Vol0_li: Optional branch-voltage linearization points; unit values are
            used when omitted.
        delta0_li: Optional branch-angle linearization points; zero values are
            used when omitted.

    Returns:
        A dictionary containing the unsolved Gurobi model and model-size
        statistics.
    """


    # Input data

    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    branch=data['branch']
    Lines=data['Lines']
    line=data['line']
    L2B=data['L2B']
    NumberL2B=data['NumberL2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']
    D2B=data['D2B']
    Gen_data=data['Gen_data']
    G=data['G']
    G2B=data['G2B']
    vmin = data['vmin']
    vmax = data['vmax']
    beta = data['beta']

    BigM_l=data['BigM_l']
    BigM_b=data['BigM_b']
    Maxdelta=data['Maxdelta']
    MaxV2=data['MaxV2']
    BigM_busbar=data['BigM_busbar']


    # Model
    model=gp.Model('BCC Benders FSP-1')

    cont_list=cont_list.copy()
    if Vol0_li is None:
        Vol0_li={(l,i,j,c):1.0 for l,i,j in Lines for c in cont_list}
    if delta0_li is None:
        delta0_li={(l,i,j,c):0.0 for l,i,j in Lines for c in cont_list}


    model.Params.OutputFlag=0


    # Variables


    Pgi=model.addVars(G,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pgi')
    Qgi=model.addVars(G,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qgi')

    Pdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pdi')
    Qdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Qdi')


    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')
    Pflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_li')
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')
    Qflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_li')
    Ploss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss')
    Qloss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss')
    Ploss_vol=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss_vol')
    Qloss_vol=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss_vol')
    Ploss_delta=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss_delta')
    Qloss_delta=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss_delta')
    epsilon=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='epsilon')

    Pflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_bus')  #From b1 to b2!
    Qflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_bus')  #From b1 to b2!


    delta_bi=model.addVars(Bus,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi')
    delta_li=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_li')
    V2_bi=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi')
    V2_li=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2_li')


    z_bus=model.addVars(Bus,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_bus')
    z_li=model.addVars(Lines,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_li')
    z_g=model.addVars(G,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_g')
    z_d=model.addVars(DemandSet,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_d')

    OF_FSP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_FSP')             # obj func var

    sp_up=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sp_up')
    sp_dn=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sp_dn')
    sq_up=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sq_up')
    sq_dn=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sq_dn')


    # Constraints


    # Active generation

    eq_Pg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar1']
                             for g in G   ), name='eq_Pg1min')
    eq_Pg1max=model.addConstrs( (    Pgi[g,'busbar1'] <=(1-z_g[g])*Gen_data.loc[g]['Pmax']
                             for g in G  ), name='eq_Pg1max')

    eq_Pg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar2']
                             for g in G   ), name='eq_Pg2min')
    eq_Pg2max=model.addConstrs( (   Pgi[g,'busbar2']  <= z_g[g]*Gen_data.loc[g]['Pmax']
                             for g in G   ), name='eq_Pg2max')

    # Reactive generation
    eq_Qg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar1']
                             for g in G   ), name='eq_Qg1min')
    eq_Qg1max=model.addConstrs( (    Qgi[g,'busbar1']<=(1-z_g[g])*Gen_data.loc[g]['Qmax']
                             for g in G  ), name='eq_Qg1max')

    eq_Qg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar2']
                             for g in G   ), name='eq_Qg2min')
    eq_Qg2max=model.addConstrs( (   Qgi[g,'busbar2'] <= z_g[g]*Gen_data.loc[g]['Qmax']
                             for g in G   ), name='eq_Pg2max')


    # Demand
    eq_Pd1=model.addConstrs( ( Pdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Pd']
                            for d in DemandSet ), name='eq_Pd1')
    eq_Pd2=model.addConstrs( ( Pdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Pd']
                            for d in DemandSet ), name='eq_Pd2')
    eq_Qd1=model.addConstrs( ( Qdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Qd']
                            for d in DemandSet ), name='eq_Qd1')
    eq_Qd2=model.addConstrs( ( Qdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Qd']
                            for d in DemandSet ), name='eq_Qd2')


    # Line-contingency constraints
    eqPflow1_cont=model.addConstrs( ( Pflow[l,i,j,c]==0
                            for l,i,j in Lines
                           for c in cont_list
                           if c==l ) , name='Eq_Pflow_cont')
    eqPflow_ei_cont=model.addConstrs( ( Pflow_li[l,i,j,m,c]==0
                            for l,i,j in Lines
                            for m in busbar
                           for c in cont_list
                           if c==l ) , name='Eq_Pflow_ei_cont')
    eqQflow1_cont=model.addConstrs( ( Qflow[l,i,j,c]==0
                            for l,i,j in Lines
                           for c in cont_list
                           if c==l ) , name='Eq_Qflow_cont')
    eqQflow_ei_cont=model.addConstrs( ( Qflow_li[l,i,j,m,c]==0
                            for l,i,j in Lines
                            for m in busbar
                           for c in cont_list
                           if c==l ) , name='Eq_Qflow_ei_cont')

    # Branch-flow limits
    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar1',c]
                             for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Pflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar2',c]
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Pflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2max')

    eq_flow=model.addConstrs((   Pflow[l,i,j,c] == Pflow_li[l,i,j,'busbar1',c]+Pflow_li[l,i,j,'busbar2',c]
                                for l,i,j in Lines for c in cont_list), name='eq_Pflow')


    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar1',c]
                             for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Qflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar2',c]
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Qflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2max')

    eq_flow=model.addConstrs((   Qflow[l,i,j,c] == Qflow_li[l,i,j,'busbar1',c]+Qflow_li[l,i,j,'busbar2',c]
                                for l,i,j in Lines for c in cont_list), name='eq_Qflow')

    # Tight linear AC limits
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij1')

    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij2')

    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij3')

    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij4')


    # Linear AC branch-flow equations

    model.addConstrs( (  Pflow[l,i,j,c] ==
    0.5*branch.loc[(l,i,j)]['g_ij']*( V2_li[l,i,j,c] - V2_li[l,j,i,c] )
    -branch.loc[(l,i,j)]['b_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) + Ploss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')

    model.addConstrs( (  Qflow[l,i,j,c] ==
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j,c] - V2_li[l,j,i,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) + Qloss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')

    # Linearized branch-loss equations
    model.addConstrs((Ploss[l,i,j,c] == Ploss_vol[l,i,j,c] + Ploss_delta[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_total')
    model.addConstrs((Qloss[l,i,j,c] == Qloss_vol[l,i,j,c] + Qloss_delta[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_total')
    model.addConstrs((0 <= Ploss[l,i,j,c] + epsilon[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_pos')
    model.addConstrs((Ploss_vol[l,i,j,c] ==
                      branch.loc[(l,i,j)]['g_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                      -0.5*branch.loc[(l,i,j)]['g_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_vol')
    model.addConstrs((Qloss_vol[l,i,j,c] ==
                      -branch.loc[(l,i,j)]['b_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                      +0.5*branch.loc[(l,i,j)]['b_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_vol')
    model.addConstrs((Ploss_delta[l,i,j,c] ==
                      branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                      -0.5*branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_delta')
    model.addConstrs((Qloss_delta[l,i,j,c] ==
                      -branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                      +0.5*branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_delta')


    # Voltage magnitude
    model.addConstrs( ( V2_bi[b,i,c] <= vmax**2 for b in Bus for i in busbar for c in cont_list) , name='eqvmaxb')
    model.addConstrs( ( V2_bi[b,i,c] >= vmin**2  for b in Bus for i in busbar for c in cont_list) , name='eqvminb')
    model.addConstrs( ( V2_li[l,i,j,c] <= vmax**2 for l,i,j in Lines for c in cont_list) , name='eqvmaxl')
    model.addConstrs( ( V2_li[l,i,j,c] >= vmin**2  for l,i,j in Lines for c in cont_list) , name='eqvminl')


    eq_delta_bus1=model.addConstrs((   -BigM_b*(1-z_bus[b]) <= delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll]
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_delta_bus1')
    eq_delta_bus2=model.addConstrs((   delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll] <= BigM_b*(1-z_bus[b])
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_delta_bus2')
    model.addConstrs(( -z_li[l,i,j]*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar1',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar1',c] <= z_li[l,i,j]*Maxdelta
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar2',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*Maxdelta
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmax')

    # Squared voltage magnitude
    eq_V2_bus1=model.addConstrs((   -MaxV2*(1-z_bus[b]) <= V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll]
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_V2_bus1')
    eq_V2_bus2=model.addConstrs((   V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll] <= MaxV2*(1-z_bus[b])
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_V2_bus2')
    model.addConstrs(( -z_li[l,i,j]*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar1',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar1',c] <= z_li[l,i,j]*MaxV2
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar2',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*MaxV2
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmax')


    eq_delta_ref=model.addConstrs((delta_bi[Bus[0],'busbar1',ll]==0     for ll in cont_list ), name='ref_bus_angle' )


    # Busbar flow
    model.addConstrs((  Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmax')
    model.addConstrs((  -Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmin')

    # Power balance
    eq_balance_busbar1=model.addConstrs((
                quicksum( Pgi[g,m]   for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m]  for d,b in D2B.select('*',b)   ) + sp_up[b,m,ll] - sp_dn[b,m,ll] ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Pflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_balance_busbar2=model.addConstrs((
                quicksum( Pgi[g,m] for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m]  for d,b in D2B.select('*',b)   ) + sp_up[b,m,ll] - sp_dn[b,m,ll] ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Pflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')


    eq_Qbalance_busbar1=model.addConstrs((
                quicksum( Qgi[g,m]   for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m]   for d,b in D2B.select('*',b)   ) + sq_up[b,m,ll] - sq_dn[b,m,ll] ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Qflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_Qbalance_busbar2=model.addConstrs((
                quicksum( Qgi[g,m] for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m]  for d,b in D2B.select('*',b)   ) + sq_up[b,m,ll] - sq_dn[b,m,ll] ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Qflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')


    model.addConstr(OF_FSP ==  quicksum( sp_up[b,m,c]+sp_dn[b,m,c]+sq_up[b,m,c]+sq_dn[b,m,c]
                                        for b in Bus for m in busbar for c in cont_list)  ,name='Eq_OF_FSP')


    model.setObjective(OF_FSP,GRB.MINIMIZE)

    model.update()

    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs


    # Model and size statistics
    return {'model':model,
            'NumBinVars':NumBinVars,
            'NumConVars':NumConVars,
            'NumConstrs':NumConstrs
            }


def solve_hmmp_line_fsp(
    data,
    TopologyMP,
    model0,
    print_result=False,
):
    """Solve a prepared line FSP for a master-problem solution.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        TopologyMP: Master-problem topology and dispatch values to impose.
        model0: Unsolved FSP created by ``create_hmmp_line_fsp``.
        print_result: Print the feasibility objective when ``True``.

    Returns:
        A dictionary containing linking-constraint duals, solve time, and the
        FSP objective value.
    """


    # Input data

    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    branch=data['branch']
    Lines=data['Lines']
    line=data['line']
    L2B=data['L2B']
    NumberL2B=data['NumberL2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']
    D2B=data['D2B']
    Gen_data=data['Gen_data']
    G=data['G']
    G2B=data['G2B']
    vmin = data['vmin']
    vmax = data['vmax']
    beta = data['beta']

    BigM_l=data['BigM_l']
    BigM_b=data['BigM_b']
    Maxdelta=data['Maxdelta']
    MaxV2=data['MaxV2']
    BigM_busbar=data['BigM_busbar']


    # Copy the reusable model and link it to the master solution.
    model=model0.copy()


    # Constraints


    model.addConstrs(  (model.getVarByName('Pgi['+str(g)+','+str(i)+']') == TopologyMP['Pgi'][g,i]     for g in G for i in busbar )  ,name='eqBendersPgi')

    model.addConstrs(  (model.getVarByName('Qgi['+str(g)+','+str(i)+']') == TopologyMP['Qgi'][g,i]     for g in G for i in busbar )  ,name='eqBendersQgi')

    model.addConstrs(  (model.getVarByName('z_bus['+str(b)+']') == (TopologyMP['bus'][b])      for b in Bus )  ,name='eqBendersZ_bus')

    model.addConstrs(  (model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')  == (TopologyMP['l_i'][l,i,j])      for l,i,j in Lines )  ,name='eqBendersZ_li')

    model.addConstrs(  (model.getVarByName('z_g['+str(g)+']') == (TopologyMP['g'][g])      for g in G )  ,name='eqBendersZ_g')

    model.addConstrs(  (model.getVarByName('z_d['+str(d)+']') == (TopologyMP['d'][d])     for d in DemandSet )  ,name='eqBendersZ_d')


    model.update()
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time  #execution time


    status = model.Status


    if status == GRB.INFEASIBLE:
        print('\n\nFSP Optimization was stopped with infeasibility!')


        # Relax bounds to diagnose an infeasible model.
        print('\n\nThe model is infeasible; relaxing the bounds\n\n')
        orignumvars = model.NumVars
        # Relax only variable bounds.
        model.feasRelaxS(0, False, True, False)

        model.optimize()

        status = model.Status
        if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
                print('The relaxed model cannot be solved \
                       because it is infeasible or unbounded')
        if status != GRB.OPTIMAL:
            print('Optimization was stopped with status %d' % status)

        # Report artificial variables introduced by the relaxation.
        print('\nSlack values:')
        slacks = model.getVars()[orignumvars:]
        for sv in slacks:
            if sv.X > 1e-9:
                print('%s = %g' % (sv.VarName, sv.X))


    # Extract linking-constraint duals.

    MuPgi={}; MuQgi={}; Muz_bus={}; Muz_li={}; Muz_g={}; Muz_d={}

    for g in G:
        for m in busbar:
            MuPgi[(g,m)]=model.getConstrByName(str('eqBendersPgi['+g+','+m+']')).pi
            MuQgi[(g,m)]=model.getConstrByName(str('eqBendersQgi['+g+','+m+']')).pi
    for b in Bus:
        Muz_bus[b]=model.getConstrByName(str('eqBendersZ_bus['+b+']')).pi

    for g in G:
        Muz_g[g]=model.getConstrByName(str('eqBendersZ_g['+g+']')).pi
    for d in DemandSet:
        Muz_d[d]=model.getConstrByName(str('eqBendersZ_d['+d+']')).pi
    for l,i,j in Lines:
        Muz_li[(l,i,j)]=model.getConstrByName(str('eqBendersZ_li['+l+','+i+','+j+']')).pi

    Mu={
        'Pgi':MuPgi,
        'Qgi':MuQgi,
        'bus':Muz_bus,
        'g':Muz_g,
        'd':Muz_d,
        'l_i':Muz_li
    }


    Mu_Pg={}; Mu_Qg={}
    for g in G:
        if model.getVarByName('z_g['+str(g)+']').x == 0:
            Mu_Pg[g]=Mu['Pgi'][g,'busbar1']
            Mu_Qg[g]=Mu['Qgi'][g,'busbar1']
        elif model.getVarByName('z_g['+str(g)+']').x == 1:
            Mu_Pg[g]=Mu['Pgi'][g,'busbar2']
            Mu_Qg[g]=Mu['Qgi'][g,'busbar2']

    Mu['Pg']=Mu_Pg
    Mu['Qg']=Mu_Qg


    if print_result==True:

        print(model.getVarByName('OF_FSP').x)


    # Results
    return {
        'Mu':Mu,
        'time':ex_time,
        'OF_FSP': model.getVarByName('OF_FSP').x
    }


def create_hmmp_substation_osp(
    data,
    cont_list=None,
    Up_redispatch=0,
    Dn_redispatch=1.0,
    Probabilistic=False,
    FixedCost=False,
    Alpha=0,
    Pg_market=None,
    Max_Sw_bus=0,
    Vol0_li=None,
    delta0_li=None,
):
    """Build the substation-contingency optimality subproblem (OSP).

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        cont_list: Contingency identifiers represented in the model.
        Up_redispatch: Available upward active-power redispatch fraction.
        Dn_redispatch: Available downward active-power redispatch fraction.
        Probabilistic: Weight shedding costs by contingency probabilities.
        FixedCost: Enforce the market-generation cost limit when ``True``.
        Alpha: Relative margin applied to the market-generation cost limit.
        Pg_market: Reference market dispatch for the fixed-cost constraint.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.
        Vol0_li: Optional branch-voltage linearization points; unit values are
            used when omitted.
        delta0_li: Optional branch-angle linearization points; zero values are
            used when omitted.

    Returns:
        A dictionary containing the unsolved Gurobi model and model-size
        statistics.
    """

    # Input data

    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    branch=data['branch']
    Lines=data['Lines']
    line=data['line']
    L2B=data['L2B']
    NumberL2B=data['NumberL2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']
    D2B=data['D2B']
    Gen_data=data['Gen_data']
    G=data['G']
    G2B=data['G2B']
    vmin = data['vmin']
    vmax = data['vmax']
    beta = data['beta']

    BigM_l=data['BigM_l']
    BigM_b=data['BigM_b']
    Maxdelta=data['Maxdelta']
    MaxV2=data['MaxV2']
    BigM_busbar=data['BigM_busbar']


    # Model
    model=gp.Model('OSP-v61')

    cont_list=cont_list.copy()
    if Vol0_li is None:
        Vol0_li={(l,i,j,c):1.0 for l,i,j in Lines for c in cont_list}
    if delta0_li is None:
        delta0_li={(l,i,j,c):0.0 for l,i,j in Lines for c in cont_list}
    model.Params.OutputFlag=0


    Pgi=model.addVars(G,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pgi')
    dPgi_up=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_up')
    dPgi_dn=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_dn')
    Qgi=model.addVars(G,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qgi')
    dQgi_up=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_up')
    dQgi_dn=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_dn')

    Pdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pdi')
    Pdi_Shed=model.addVars(DemandSet,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Pdi_Shed')
    Qdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Qdi')
    Qdi_Shed=model.addVars(DemandSet,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Qdi_Shed')

    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')
    Pflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_li')
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')
    Qflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_li')

    Ploss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss')
    Qloss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss')
    Ploss_vol=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss_vol')
    Qloss_vol=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss_vol')
    Ploss_delta=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss_delta')
    Qloss_delta=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss_delta')
    epsilon=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='epsilon')

    Pflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_bus')  #From b1 to b2!
    Qflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_bus')  #From b1 to b2!

    delta_bi=model.addVars(Bus,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi')
    delta_li=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_li')
    V2_bi=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi')
    V2_li=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2_li')


    z_bus=model.addVars(Bus,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_bus')
    z_li=model.addVars(Lines,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_li')
    z_g=model.addVars(G,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_g')
    z_d=model.addVars(DemandSet,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_d')

    OF_OSP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_OSP')             # obj func var
    RDCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='RDCost')
    ShedCost=model.addVars(cont_list,lb=0,vtype=GRB.CONTINUOUS,name='ShedCost')
    Shed_c=model.addVars(cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Shed_c')
    TotalShedCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalShedCost')


    # Constraints


    # Topology and switching constraints

    # Reliability switching
    eq_2line_busbar2=model.addConstrs(( 2*(1-z_bus[b]) <=
                                   quicksum(z_li[l,i,j] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus ), name='eq_2line_busbar2' )

    eq_2line_busbar1=model.addConstrs(( 2*(1-z_bus[b]) <=
                                   quicksum(1-z_li[l,i,j] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus ), name='eq_2line_busbar1' )

    eq_2line_busbar3=model.addConstrs(( z_bus[b] == 1
                            for b in Bus
                             if NumberL2B[b]<=3  ), name='eq_2line_busbar3' )

    # Maximum number of busbar splits
    eq_MaxSw_bus=model.addConstr((  quicksum( (1-z_bus[b]) for b in Bus) <= Max_Sw_bus ) , name='eq_MaxSw_bus')        # Eq max line sw


    # Symmetry breaking
    for b in Bus:
        lmin,b=L2B.select('*',b)[0]
        lmin,i,j=Lines.select(lmin,b,'*')[0]
        eq_symmetry1=model.addConstr((  z_li[lmin,i,j] ==0  ), name='eq_symmetry')

    eq_zbus_ref=model.addConstr( (  z_bus[Bus[0]]==1    ), name='eq_zbus_ref')


    # Active generation
    eq_Pg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar1']+dPgi_up[g,'busbar1',c]-dPgi_dn[g,'busbar1',c]
                             for g in G for c in cont_list  ), name='eq_Pg1min')
    eq_Pg1max=model.addConstrs( (    Pgi[g,'busbar1']+dPgi_up[g,'busbar1',c]-dPgi_dn[g,'busbar1',c] <=(1-z_g[g])*Gen_data.loc[g]['Pmax']
                             for g in G for c in cont_list ), name='eq_Pg1max')

    eq_Pg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar2'] +dPgi_up[g,'busbar2',c]-dPgi_dn[g,'busbar2',c]
                             for g in G for c in cont_list  ), name='eq_Pg2min')
    eq_Pg2max=model.addConstrs( (   Pgi[g,'busbar2']+dPgi_up[g,'busbar2',c]-dPgi_dn[g,'busbar2',c]  <= z_g[g]*Gen_data.loc[g]['Pmax']
                             for g in G for c in cont_list  ), name='eq_Pg2max')

    # Reactive generation
    eq_Qg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c]
                             for g in G for c in cont_list  ), name='eq_Qg1min')
    eq_Qg1max=model.addConstrs( (    Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c] <=(1-z_g[g])*Gen_data.loc[g]['Qmax']
                             for g in G for c in cont_list ), name='eq_Qg1max')

    eq_Qg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar2'] +dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]
                             for g in G for c in cont_list  ), name='eq_Qg2min')
    eq_Qg2max=model.addConstrs( (   Qgi[g,'busbar2']+dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]  <= z_g[g]*Gen_data.loc[g]['Qmax']
                             for g in G for c in cont_list  ), name='eq_Pg2max')


    # Corrective redispatch

    Eq_dPgUp_res=model.addConstrs(( dPgi_up[g,i,c] <= Gen_data.loc[g]['Pmax']*Up_redispatch
                          for g in G for i in busbar for c in cont_list ),name='Eq_dPgUp_res')
    Eq_dPgDn_res=model.addConstrs(( dPgi_dn[g,i,c] <= Gen_data.loc[g]['Pmax']*Dn_redispatch
                            for g in G for i in busbar for c in cont_list ),name='Eq_dPgDn_res')
    model.addConstrs(( dPgi_up[g,i,c] == 0
                          for g in G for i in busbar for c in cont_list if c==0 ),name='Eq_dPg0')

    Eq_dQgUp_res=model.addConstrs(( dQgi_up[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #*Up_redispatch
                          for g in G for i in busbar for c in cont_list ),name='Eq_dQgUp_res')
    Eq_dQgDn_res=model.addConstrs(( dQgi_dn[g,i,c] <= Gen_data.loc[g]['Qmax']*1# *Dn_redispatch
                            for g in G for i in busbar for c in cont_list ),name='Eq_dQgDn_res')


    # Demand
    eq_Pd1=model.addConstrs( ( Pdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Pd']
                            for d in DemandSet ), name='eq_Pd1')
    eq_Pd2=model.addConstrs( ( Pdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Pd']
                            for d in DemandSet ), name='eq_Pd2')
    eq_Qd1=model.addConstrs( ( Qdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Qd']
                            for d in DemandSet ), name='eq_Qd1')
    eq_Qd2=model.addConstrs( ( Qdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Qd']
                            for d in DemandSet ), name='eq_Qd2')

    # Load shedding
    model.addConstrs( (  Pdi_Shed[d,'busbar1',c]==Pdi[d,'busbar1']
                              for b in Bus for d,b in D2B.select('*',b)
                               for c in cont_list if c==str(b+'-1') ) , name='eq_Pd1Shed')
    model.addConstrs( (  Pdi_Shed[d,'busbar2',c]==Pdi[d,'busbar2']
                              for b in Bus for d,b in D2B.select('*',b)
                               for c in cont_list if c==str(b+'-2') ) , name='eq_Pd2Shed')
    model.addConstrs( (  Pdi_Shed[d,i,c]<=Pdi[d,i]
                              for d in DemandSet for i in busbar
                               for c in cont_list) , name='eq_PdShedlimit')

    model.addConstrs( (  Qdi_Shed[d,i,c]== (Pdemand.loc[d]['Qd'])/(Pdemand.loc[d]['Pd'])*Pdi_Shed[d,i,c]
                              for d in DemandSet for i in busbar
                               for c in cont_list) , name='eq_QdShed')


    # Line-contingency constraints
    eqPflow1_cont=model.addConstrs( ( Pflow[l,i,j,c]==0
                            for l,i,j in Lines
                           for c in cont_list
                           if c==l ) , name='Eq_Pflow_cont')
    eqPflow_ei_cont=model.addConstrs( ( Pflow_li[l,i,j,m,c]==0
                            for l,i,j in Lines
                            for m in busbar
                           for c in cont_list
                           if c==l ) , name='Eq_Pflow_ei_cont')
    eqQflow1_cont=model.addConstrs( ( Qflow[l,i,j,c]==0
                            for l,i,j in Lines
                           for c in cont_list
                           if c==l ) , name='Eq_Qflow_cont')
    eqQflow_ei_cont=model.addConstrs( ( Qflow_li[l,i,j,m,c]==0
                            for l,i,j in Lines
                            for m in busbar
                           for c in cont_list
                           if c==l ) , name='Eq_Qflow_ei_cont')

    # Branch-flow limits
    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar1',c]
                             for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Pflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar2',c]
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Pflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2max')

    eq_flow=model.addConstrs((   Pflow[l,i,j,c] == Pflow_li[l,i,j,'busbar1',c]+Pflow_li[l,i,j,'busbar2',c]
                                for l,i,j in Lines for c in cont_list), name='eq_Pflow')


    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar1',c]
                             for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Qflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar2',c]
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Qflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2max')

    eq_flow=model.addConstrs((   Qflow[l,i,j,c] == Qflow_li[l,i,j,'busbar1',c]+Qflow_li[l,i,j,'busbar2',c]
                                for l,i,j in Lines for c in cont_list), name='eq_Qflow')

    # Tight linear AC limits
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij1')

    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij2')

    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij3')

    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij4')


    # Linear AC branch-flow equations

    model.addConstrs( (  Pflow[l,i,j,c] ==
    0.5*branch.loc[(l,i,j)]['g_ij']*( V2_li[l,i,j,c] - V2_li[l,j,i,c] )
    -branch.loc[(l,i,j)]['b_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) + Ploss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')

    model.addConstrs( (  Qflow[l,i,j,c] ==
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j,c] - V2_li[l,j,i,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) + Qloss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')

    # Linearized branch-loss equations
    model.addConstrs((Ploss[l,i,j,c] == Ploss_vol[l,i,j,c] + Ploss_delta[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_total')
    model.addConstrs((Qloss[l,i,j,c] == Qloss_vol[l,i,j,c] + Qloss_delta[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_total')
    model.addConstrs((0 <= Ploss[l,i,j,c] + epsilon[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_pos')
    model.addConstrs((Ploss_vol[l,i,j,c] ==
                      branch.loc[(l,i,j)]['g_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                      -0.5*branch.loc[(l,i,j)]['g_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_vol')
    model.addConstrs((Qloss_vol[l,i,j,c] ==
                      -branch.loc[(l,i,j)]['b_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                      +0.5*branch.loc[(l,i,j)]['b_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_vol')
    model.addConstrs((Ploss_delta[l,i,j,c] ==
                      branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                      -0.5*branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_delta')
    model.addConstrs((Qloss_delta[l,i,j,c] ==
                      -branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                      +0.5*branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_delta')


    model.addConstrs( ( V2_bi[b,i,c] <= vmax**2  for b in Bus for i in busbar for c in cont_list) , name='eqvmaxb')
    model.addConstrs( ( V2_bi[b,i,c] >= vmin**2  for b in Bus for i in busbar for c in cont_list) , name='eqvminb')
    model.addConstrs( ( V2_li[l,i,j,c] <= vmax**2 for l,i,j in Lines for c in cont_list) , name='eqvmaxl')
    model.addConstrs( ( V2_li[l,i,j,c] >= vmin**2  for l,i,j in Lines for c in cont_list) , name='eqvminl')


    eq_delta_bus1=model.addConstrs((   -BigM_b*(1-z_bus[b]) <= delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll]
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_delta_bus1')
    eq_delta_bus2=model.addConstrs((   delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll] <= BigM_b*(1-z_bus[b])
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_delta_bus2')
    model.addConstrs(( -z_li[l,i,j]*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar1',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar1',c] <= z_li[l,i,j]*Maxdelta
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar2',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*Maxdelta
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmax')

    # Squared voltage magnitude
    eq_V2_bus1=model.addConstrs((   -MaxV2*(1-z_bus[b]) <= V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll]
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_V2_bus1')
    eq_V2_bus2=model.addConstrs((   V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll] <= MaxV2*(1-z_bus[b])
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_V2_bus2')
    model.addConstrs(( -z_li[l,i,j]*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar1',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar1',c] <= z_li[l,i,j]*MaxV2
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar2',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*MaxV2
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmax')


    eq_delta_ref=model.addConstrs((delta_bi[Bus[0],'busbar1',ll]==0     for ll in cont_list ), name='ref_bus_angle' )


    # Busbar and coupler contingencies
    model.addConstrs((  Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmax')
    model.addConstrs((  -Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmin')

    model.addConstrs((  Pflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==b ), name='eq_CouplerCont')
    model.addConstrs((  Pflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-1') ), name='eq_bus1Cont')
    model.addConstrs((  Pflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-2') ), name='eq_bus2Cont')

    model.addConstrs( ( Pflow_li[l,i,j,'busbar1',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-1') ) , name='Eq_Pflow_ei_bus1contFr')
    model.addConstrs( ( Pflow_li[l,i,j,'busbar2',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-2') ) , name='Eq_Pflow_ei_bus2contFr')

    model.addConstrs((  Qflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarQflowmax')
    model.addConstrs((  -Qflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarQflowmin')

    model.addConstrs((  Qflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==b ), name='eq_CouplerCont')
    model.addConstrs((  Qflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-1') ), name='eq_bus1Cont')
    model.addConstrs((  Qflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-2') ), name='eq_bus2Cont')

    model.addConstrs( ( Qflow_li[l,i,j,'busbar1',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-1') ) , name='Eq_Qflow_ei_bus1contFr')
    model.addConstrs( ( Qflow_li[l,i,j,'busbar2',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-2') ) , name='Eq_Qflow_ei_bus2contFr')


    # Power balance

    eq_balance_busbar1=model.addConstrs((
                quicksum( Pgi[g,m] +dPgi_up[g,m,ll]-dPgi_dn[g,m,ll]  for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m] - Pdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Pflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_balance_busbar2=model.addConstrs((
                quicksum( Pgi[g,m]+dPgi_up[g,m,ll]-dPgi_dn[g,m,ll] for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m]- Pdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Pflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')


    eq_Qbalance_busbar1=model.addConstrs((
                quicksum( Qgi[g,m] +dQgi_up[g,m,ll]-dQgi_dn[g,m,ll]  for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m] - Qdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Qflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_Qbalance_busbar2=model.addConstrs((
                quicksum( Qgi[g,m]+dQgi_up[g,m,ll]-dQgi_dn[g,m,ll] for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m]- Qdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Qflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')


    # Objective components
    model.addConstr(RDCost == quicksum( Gen_data.loc[g]['c_res']*(dPgi_up[g,i,c]+0.0001*dPgi_dn[g,i,c])
                                       for g in G for i in busbar for c in cont_list)
                              ,name='Eq_RD')

    model.addConstrs( (Shed_c[c] == quicksum(Pdi_Shed[d,i,c]
                                              for d in DemandSet for i in busbar)     for c in cont_list) ,name='Eq_Shed_c')

    if Probabilistic==True:
        model.addConstrs( (ShedCost[c] == quicksum(Pdemand.loc[d]['ShedCost']*data['Prob_cont'][c]*Pdi_Shed[d,i,c]
                                              for d in DemandSet for i in busbar)     for c in cont_list) ,name='Eq_ShedC')
    else:
        model.addConstrs( (ShedCost[c] == quicksum(Pdemand.loc[d]['ShedCost']*Pdi_Shed[d,i,c]
                                              for d in DemandSet for i in busbar)     for c in cont_list) ,name='Eq_ShedC')


    model.addConstr(TotalShedCost ==  quicksum(ShedCost[c]  for c in cont_list)   ,name='Eq_TotalShedCost')


    model.addConstr(OF_OSP == RDCost + TotalShedCost ,name='Eq_OF_OSP')


    model.setObjective(OF_OSP,GRB.MINIMIZE)

    model.update()


    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs


    # Model and size statistics
    return {'model':model,
            'NumBinVars':NumBinVars,
            'NumConVars':NumConVars,
            'NumConstrs':NumConstrs
            }


def solve_hmmp_substation_osp(
    data,
    TopologyMP,
    model0,
    cont_list=None,
):
    """Solve a prepared OSP for a master-problem topology.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        TopologyMP: Master-problem topology and dispatch values to impose.
        model0: Unsolved OSP created by ``create_hmmp_substation_osp``.
        cont_list: Contingency identifiers represented by ``model0``.

    Returns:
        A dictionary containing linking duals, solve time, OSP cost, and
        contingency- and demand-level load shedding, together with updated
        voltage and angle linearization points.
    """

    # Input data

    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    branch=data['branch']
    Lines=data['Lines']
    line=data['line']
    L2B=data['L2B']
    NumberL2B=data['NumberL2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']
    D2B=data['D2B']
    Gen_data=data['Gen_data']
    G=data['G']
    G2B=data['G2B']
    vmin = data['vmin']
    vmax = data['vmax']
    beta = data['beta']

    BigM_l=data['BigM_l']
    BigM_b=data['BigM_b']
    Maxdelta=data['Maxdelta']
    MaxV2=data['MaxV2']
    BigM_busbar=data['BigM_busbar']


    # Copy the reusable model and link it to the master solution.
    model=model0.copy()


    # Constraints


    model.addConstrs(  (model.getVarByName('Pgi['+str(g)+','+str(i)+']') == TopologyMP['Pgi'][g,i]     for g in G for i in busbar )  ,name='eqBendersPgi')

    model.addConstrs(  (model.getVarByName('Qgi['+str(g)+','+str(i)+']') == TopologyMP['Qgi'][g,i]     for g in G for i in busbar )  ,name='eqBendersQgi')

    model.addConstrs(  (model.getVarByName('z_bus['+str(b)+']') == (TopologyMP['bus'][b])      for b in Bus )  ,name='eqBendersZ_bus')

    model.addConstrs(  (model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')  == (TopologyMP['l_i'][l,i,j])      for l,i,j in Lines )  ,name='eqBendersZ_li')

    model.addConstrs(  (model.getVarByName('z_g['+str(g)+']') == (TopologyMP['g'][g])      for g in G )  ,name='eqBendersZ_g')

    model.addConstrs(  (model.getVarByName('z_d['+str(d)+']') == (TopologyMP['d'][d])     for d in DemandSet )  ,name='eqBendersZ_d')


    model.update()
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time


    # Extract linking-constraint duals.


    MuPgi={}; MuQgi={}; Muz_bus={}; Muz_li={}; Muz_g={}; Muz_d={}

    for g in G:
        for m in busbar:
            MuPgi[(g,m)]=model.getConstrByName(str('eqBendersPgi['+g+','+m+']')).pi
            MuQgi[(g,m)]=model.getConstrByName(str('eqBendersQgi['+g+','+m+']')).pi
    for b in Bus:
        Muz_bus[b]=model.getConstrByName(str('eqBendersZ_bus['+b+']')).pi
    for g in G:
        Muz_g[g]=model.getConstrByName(str('eqBendersZ_g['+g+']')).pi
    for d in DemandSet:
        Muz_d[d]=model.getConstrByName(str('eqBendersZ_d['+d+']')).pi
    for l,i,j in Lines:
        Muz_li[(l,i,j)]=model.getConstrByName(str('eqBendersZ_li['+l+','+i+','+j+']')).pi

    Mu={
        'Pgi':MuPgi,
        'Qgi':MuQgi,
        'bus':Muz_bus,
        'g':Muz_g,
        'd':Muz_d,
        'l_i':Muz_li
    }

    Mu_Pg={}; Mu_Qg={}
    for g in G:
        if model.getVarByName('z_g['+str(g)+']').x == 0:
            Mu_Pg[g]=Mu['Pgi'][g,'busbar1']
            Mu_Qg[g]=Mu['Qgi'][g,'busbar1']
        elif model.getVarByName('z_g['+str(g)+']').x == 1:
            Mu_Pg[g]=Mu['Pgi'][g,'busbar2']
            Mu_Qg[g]=Mu['Qgi'][g,'busbar2']

    Mu['Pg']=Mu_Pg
    Mu['Qg']=Mu_Qg


    ShedCost_df = pd.DataFrame(columns=['ShedCost(c)'],index=cont_list)
    PdShed_dic = pd.DataFrame(columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,busbar,cont_list]))
    for c in cont_list:
        ShedCost_df.loc[c]=model.getVarByName('ShedCost['+str(c)+']').x
        for d in DemandSet:
            for i in busbar:
                PdShed_dic.loc[d,i,c]=model.getVarByName('Pdi_Shed['+str(d)+','+str(i)+','+str(c)+']').x

    Vol0_li, delta0_li = _extract_linearization_points(
        model,
        Lines,
        cont_list,
    )


    # Results
    return {
        'Mu':Mu,
        'time':ex_time,
        'OF_OSP': model.getVarByName('OF_OSP').x,
        'ShedCost_df':ShedCost_df,
        'PdShed_dic':PdShed_dic,
        'Vol0_li':Vol0_li,
        'delta0_li':delta0_li
    }


def create_hmmp_substation_mp(
    data,
    substation,
    Probabilistic=False,
    FixedCost=False,
    Alpha=0,
    Pg_market=None,
    Max_Sw_bus=0,
    Up_redispatch=0,
    Dn_redispatch=1.0,
    Vol0_li=None,
    delta0_li=None,
):
    """Build the second master problem (MP2) for one substation.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        substation: Substation whose topology is optimized by this MP2.
        Probabilistic: Weight shedding costs by contingency probabilities.
        FixedCost: Enforce the market-generation cost limit when ``True``.
        Alpha: Relative margin applied to the market-generation cost limit.
        Pg_market: Reference market dispatch for the fixed-cost constraint.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.
        Up_redispatch: Available upward active-power redispatch fraction.
        Dn_redispatch: Available downward active-power redispatch fraction.
        Vol0_li: Optional branch-voltage linearization points; unit values are
            used when omitted.
        delta0_li: Optional branch-angle linearization points; zero values are
            used when omitted.

    Returns:
        A dictionary containing the unsolved substation MP2 and model-size
        statistics.
    """

    # Input data

    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    branch=data['branch']
    Lines=data['Lines']
    line=data['line']
    L2B=data['L2B']
    NumberL2B=data['NumberL2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']
    D2B=data['D2B']
    Gen_data=data['Gen_data']
    G=data['G']
    G2B=data['G2B']
    vmin = data['vmin']
    vmax = data['vmax']
    beta = data['beta']

    BigM_l=data['BigM_l']
    BigM_b=data['BigM_b']
    Maxdelta=data['Maxdelta']
    MaxV2=data['MaxV2']
    BigM_busbar=data['BigM_busbar']


    # Model
    model=gp.Model('Benders_MP2_cont_v61')


    cont_list = [0,substation,str(substation+'-1'),str(substation+'-2')]  #c0 is included!
    if Vol0_li is None:
        Vol0_li={(l,i,j,c):1.0 for l,i,j in Lines for c in cont_list}
    if delta0_li is None:
        delta0_li={(l,i,j,c):0.0 for l,i,j in Lines for c in cont_list}


    model.Params.OutputFlag=0


    # Variables


    Pgi=model.addVars(G,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pgi')
    dPgi_up=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_up')
    dPgi_dn=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_dn')
    Qgi=model.addVars(G,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qgi')
    dQgi_up=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_up')
    dQgi_dn=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_dn')

    Pdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pdi')
    Pdi_Shed=model.addVars(DemandSet,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Pdi_Shed')
    Qdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Qdi')
    Qdi_Shed=model.addVars(DemandSet,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Qdi_Shed')

    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')
    Pflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_li')
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')
    Qflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_li')

    Ploss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss')
    Qloss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss')
    Ploss_vol=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss_vol')
    Qloss_vol=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss_vol')
    Ploss_delta=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss_delta')
    Qloss_delta=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss_delta')
    epsilon=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='epsilon')

    Pflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_bus')  #From b1 to b2!
    Qflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_bus')  #From b1 to b2!

    delta_bi=model.addVars(Bus,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi')
    delta_li=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_li')
    V2_bi=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi')
    V2_li=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2_li')


    z_bus=model.addVars(Bus,vtype=GRB.BINARY,name='z_bus')
    z_li=model.addVars(Lines,vtype=GRB.BINARY,name='z_li')
    z_g=model.addVars(G,vtype=GRB.BINARY,name='z_g')
    z_d=model.addVars(DemandSet,vtype=GRB.BINARY,name='z_d')

    GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost')
    RDCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='RDCost')
    ShedCost=model.addVars(cont_list,lb=0,vtype=GRB.CONTINUOUS,name='ShedCost')
    TotalShedCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalShedCost')

    OF_MP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_MP')


    # Topology and switching constraints
    eq_2line_busbar2=model.addConstrs(( 2*(1-z_bus[b]) <=
                                   quicksum(z_li[l,i,j] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus ), name='eq_2line_busbar2' )
    eq_2line_busbar1=model.addConstrs(( 2*(1-z_bus[b]) <=
                                   quicksum(1-z_li[l,i,j] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus ), name='eq_2line_busbar1' )


    eq_2line_busbar3=model.addConstrs(( z_bus[b] == 1
                            for b in Bus
                             if NumberL2B[b]<=3  ), name='eq_2line_busbar3' )

    # Maximum number of busbar splits
    eq_MaxSw_bus=model.addConstr((  quicksum( (1-z_bus[b]) for b in Bus) <= Max_Sw_bus) , name='eq_MaxSw_bus')

    # Symmetry breaking
    for b in Bus:
        lmin,b=L2B.select('*',b)[0]
        lmin,i,j=Lines.select(lmin,b,'*')[0]
        eq_symmetry1=model.addConstr((  z_li[lmin,i,j] ==0  ), name='eq_symmetry1')

    eq_zbus_ref=model.addConstr( (  z_bus[Bus[0]]==1    ), name='eq_zbus_ref')

    # Active generation

    eq_Pg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar1']+dPgi_up[g,'busbar1',c]-dPgi_dn[g,'busbar1',c]
                             for g in G for c in cont_list  ), name='eq_Pg1min')
    eq_Pg1max=model.addConstrs( (    Pgi[g,'busbar1']+dPgi_up[g,'busbar1',c]-dPgi_dn[g,'busbar1',c] <=(1-z_g[g])*Gen_data.loc[g]['Pmax']
                             for g in G for c in cont_list ), name='eq_Pg1max')

    eq_Pg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar2'] +dPgi_up[g,'busbar2',c]-dPgi_dn[g,'busbar2',c]
                             for g in G for c in cont_list  ), name='eq_Pg2min')
    eq_Pg2max=model.addConstrs( (   Pgi[g,'busbar2']+dPgi_up[g,'busbar2',c]-dPgi_dn[g,'busbar2',c]  <= z_g[g]*Gen_data.loc[g]['Pmax']
                             for g in G for c in cont_list  ), name='eq_Pg2max')

    # Reactive generation
    eq_Qg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c]
                             for g in G for c in cont_list  ), name='eq_Qg1min')
    eq_Qg1max=model.addConstrs( (    Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c] <=(1-z_g[g])*Gen_data.loc[g]['Qmax']
                             for g in G for c in cont_list ), name='eq_Qg1max')

    eq_Qg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar2'] +dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]
                             for g in G for c in cont_list  ), name='eq_Qg2min')
    eq_Qg2max=model.addConstrs( (   Qgi[g,'busbar2']+dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]  <= z_g[g]*Gen_data.loc[g]['Qmax']
                             for g in G for c in cont_list  ), name='eq_Pg2max')


    # Corrective redispatch

    Eq_dPgUp_res=model.addConstrs(( dPgi_up[g,i,c] <= Gen_data.loc[g]['Pmax']*Up_redispatch
                          for g in G for i in busbar for c in cont_list ),name='Eq_dPgUp_res')
    Eq_dPgDn_res=model.addConstrs(( dPgi_dn[g,i,c] <= Gen_data.loc[g]['Pmax']*Dn_redispatch
                            for g in G for i in busbar for c in cont_list ),name='Eq_dPgDn_res')
    model.addConstrs(( dPgi_up[g,i,c] == 0
                          for g in G for i in busbar for c in cont_list if c==0 ),name='Eq_dPg0')

    Eq_dQgUp_res=model.addConstrs(( dQgi_up[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #*Up_redispatch
                          for g in G for i in busbar for c in cont_list ),name='Eq_dQgUp_res')
    Eq_dQgDn_res=model.addConstrs(( dQgi_dn[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #Dn_redispatch
                            for g in G for i in busbar for c in cont_list ),name='Eq_dQgDn_res')


    # Demand

    eq_Pd1=model.addConstrs( ( Pdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Pd']
                            for d in DemandSet ), name='eq_Pd1')
    eq_Pd2=model.addConstrs( ( Pdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Pd']
                            for d in DemandSet ), name='eq_Pd2')
    eq_Qd1=model.addConstrs( ( Qdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Qd']
                            for d in DemandSet ), name='eq_Qd1')
    eq_Qd2=model.addConstrs( ( Qdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Qd']
                            for d in DemandSet ), name='eq_Qd2')

    # Load shedding
    model.addConstrs( (  Pdi_Shed[d,'busbar1',c]==Pdi[d,'busbar1']
                              for b in Bus for d,b in D2B.select('*',b)
                               for c in cont_list if c==str(b+'-1') ) , name='eq_Pd1Shed')
    model.addConstrs( (  Pdi_Shed[d,'busbar2',c]==Pdi[d,'busbar2']
                              for b in Bus for d,b in D2B.select('*',b)
                               for c in cont_list if c==str(b+'-2') ) , name='eq_Pd2Shed')
    model.addConstrs( (  Pdi_Shed[d,i,c]<=Pdi[d,i]
                              for d in DemandSet for i in busbar
                               for c in cont_list) , name='eq_PdShedlimit')

    model.addConstrs( (  Qdi_Shed[d,i,c]== (Pdemand.loc[d]['Qd'])/(Pdemand.loc[d]['Pd'])*Pdi_Shed[d,i,c]
                              for d in DemandSet for i in busbar
                               for c in cont_list) , name='eq_QdShed')


    # Line-contingency constraints
    eqPflow1_cont=model.addConstrs( ( Pflow[l,i,j,c]==0
                            for l,i,j in Lines
                           for c in cont_list
                           if c==l ) , name='Eq_Pflow_cont')
    eqPflow_ei_cont=model.addConstrs( ( Pflow_li[l,i,j,m,c]==0
                            for l,i,j in Lines
                            for m in busbar
                           for c in cont_list
                           if c==l ) , name='Eq_Pflow_ei_cont')
    eqQflow1_cont=model.addConstrs( ( Qflow[l,i,j,c]==0
                            for l,i,j in Lines
                           for c in cont_list
                           if c==l ) , name='Eq_Qflow_cont')
    eqQflow_ei_cont=model.addConstrs( ( Qflow_li[l,i,j,m,c]==0
                            for l,i,j in Lines
                            for m in busbar
                           for c in cont_list
                           if c==l ) , name='Eq_Qflow_ei_cont')

    # Branch-flow limits
    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar1',c]
                             for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Pflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar2',c]
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2min')

    eq_flow2max=model.addConstrs( ( Pflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2max')

    eq_flow=model.addConstrs((   Pflow[l,i,j,c] == Pflow_li[l,i,j,'busbar1',c]+Pflow_li[l,i,j,'busbar2',c]
                                for l,i,j in Lines for c in cont_list), name='eq_Pflow')


    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar1',c]
                             for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Qflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar2',c]
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Qflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines  for c in cont_list ) ,name='eq_flow2max')

    eq_flow=model.addConstrs((   Qflow[l,i,j,c] == Qflow_li[l,i,j,'busbar1',c]+Qflow_li[l,i,j,'busbar2',c]
                                for l,i,j in Lines for c in cont_list), name='eq_Qflow')

    # Tight linear AC limits
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij1')

    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij2')

    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij3')

    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij4')


    # Linear AC branch-flow equations

    model.addConstrs( (  Pflow[l,i,j,c] ==
    0.5*branch.loc[(l,i,j)]['g_ij']*( V2_li[l,i,j,c] - V2_li[l,j,i,c] )
    -branch.loc[(l,i,j)]['b_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) + Ploss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')

    model.addConstrs( (  Qflow[l,i,j,c] ==
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j,c] - V2_li[l,j,i,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) + Qloss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')

    # Linearized branch-loss equations
    model.addConstrs((Ploss[l,i,j,c] == Ploss_vol[l,i,j,c] + Ploss_delta[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_total')
    model.addConstrs((Qloss[l,i,j,c] == Qloss_vol[l,i,j,c] + Qloss_delta[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_total')
    model.addConstrs((0 <= Ploss[l,i,j,c] + epsilon[l,i,j,c]
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_pos')
    model.addConstrs((Ploss_vol[l,i,j,c] ==
                      branch.loc[(l,i,j)]['g_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                      -0.5*branch.loc[(l,i,j)]['g_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_vol')
    model.addConstrs((Qloss_vol[l,i,j,c] ==
                      -branch.loc[(l,i,j)]['b_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                      +0.5*branch.loc[(l,i,j)]['b_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_vol')
    model.addConstrs((Ploss_delta[l,i,j,c] ==
                      branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                      -0.5*branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqPloss_delta')
    model.addConstrs((Qloss_delta[l,i,j,c] ==
                      -branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                      +0.5*branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                      for l,i,j in Lines for c in cont_list if l!=c), name='eqQloss_delta')


    # Voltage magnitude
    model.addConstrs( ( V2_bi[b,i,c] <= vmax**2  for b in Bus for i in busbar for c in cont_list) , name='eqvmaxb')
    model.addConstrs( ( V2_bi[b,i,c] >= vmin**2  for b in Bus for i in busbar for c in cont_list) , name='eqvminb')
    model.addConstrs( ( V2_li[l,i,j,c] <= vmax**2 for l,i,j in Lines for c in cont_list) , name='eqvmaxl')
    model.addConstrs( ( V2_li[l,i,j,c] >= vmin**2  for l,i,j in Lines for c in cont_list) , name='eqvminl')

    eq_delta_bus1=model.addConstrs((   -BigM_b*(1-z_bus[b]) <= delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll]
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_delta_bus1')
    eq_delta_bus2=model.addConstrs((   delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll] <= BigM_b*(1-z_bus[b])
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_delta_bus2')
    model.addConstrs(( -z_li[l,i,j]*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar1',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar1',c] <= z_li[l,i,j]*Maxdelta
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar2',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*Maxdelta
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmax')

    # Squared voltage magnitude
    eq_V2_bus1=model.addConstrs((   -MaxV2*(1-z_bus[b]) <= V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll]
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_V2_bus1')
    eq_V2_bus2=model.addConstrs((   V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll] <= MaxV2*(1-z_bus[b])
                               for b in Bus for ll in cont_list
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_V2_bus2')
    model.addConstrs(( -z_li[l,i,j]*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar1',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar1',c] <= z_li[l,i,j]*MaxV2
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar2',c]
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*MaxV2
                                for l,i,j in Lines
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmax')


    eq_delta_ref=model.addConstrs((delta_bi[Bus[0],'busbar1',ll]==0     for ll in cont_list ), name='ref_bus_angle' )


    # Busbar and coupler contingencies
    model.addConstrs((  Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmax')
    model.addConstrs((  -Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmin')

    model.addConstrs((  Pflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==b ), name='eq_CouplerCont')
    model.addConstrs((  Pflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-1') ), name='eq_bus1Cont')
    model.addConstrs((  Pflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-2') ), name='eq_bus2Cont')

    model.addConstrs( ( Pflow_li[l,i,j,'busbar1',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-1') ) , name='Eq_Pflow_ei_bus1contFr')
    model.addConstrs( ( Pflow_li[l,i,j,'busbar2',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-2') ) , name='Eq_Pflow_ei_bus2contFr')

    model.addConstrs((  Qflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarQflowmax')
    model.addConstrs((  -Qflow_bus[b,ll]<=BigM_busbar*(z_bus[b])
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarQflowmin')

    model.addConstrs((  Qflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==b ), name='eq_CouplerCont')
    model.addConstrs((  Qflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-1') ), name='eq_bus1Cont')
    model.addConstrs((  Qflow_bus[b,ll]==0
                          for b in Bus for ll in cont_list if ll==str(b+'-2') ), name='eq_bus2Cont')

    model.addConstrs( ( Qflow_li[l,i,j,'busbar1',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-1') ) , name='Eq_Qflow_ei_bus1contFr')
    model.addConstrs( ( Qflow_li[l,i,j,'busbar2',c]==0
                            for b in Bus for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list
                           if c==str(b+'-2') ) , name='Eq_Qflow_ei_bus2contFr')


    eq_balance_busbar1=model.addConstrs((
                quicksum( Pgi[g,m] +dPgi_up[g,m,ll]-dPgi_dn[g,m,ll]  for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m] - Pdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Pflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_balance_busbar2=model.addConstrs((
                quicksum( Pgi[g,m]+dPgi_up[g,m,ll]-dPgi_dn[g,m,ll] for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m]- Pdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Pflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')


    eq_Qbalance_busbar1=model.addConstrs((
                quicksum( Qgi[g,m] +dQgi_up[g,m,ll]-dQgi_dn[g,m,ll]  for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m] - Qdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Qflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_Qbalance_busbar2=model.addConstrs((
                quicksum( Qgi[g,m]+dQgi_up[g,m,ll]-dQgi_dn[g,m,ll] for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m]- Qdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Qflow_bus[b,ll]
                            for b in Bus for m in busbar for ll in cont_list
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')


    # Objective components
    model.addConstr(GenCost == quicksum(Gen_data.loc[g]['b']*Pgi[g,i] for g in G for i in busbar) ,name='Eq_GC')


    model.addConstr(RDCost == quicksum( Gen_data.loc[g]['c_res']*(dPgi_up[g,i,c]+0.0001*dPgi_dn[g,i,c])
                                       for g in G for i in busbar for c in cont_list)
                              ,name='Eq_RD')


    if Probabilistic==True:
        model.addConstrs( (ShedCost[c] == quicksum(Pdemand.loc[d]['ShedCost']*data['Prob_cont'][c]*Pdi_Shed[d,i,c]
                                              for d in DemandSet for i in busbar)     for c in cont_list) ,name='Eq_ShedC')
    else:
        model.addConstrs( (ShedCost[c] == quicksum(Pdemand.loc[d]['ShedCost']*Pdi_Shed[d,i,c]
                                              for d in DemandSet for i in busbar)     for c in cont_list) ,name='Eq_ShedC')


    model.addConstr(TotalShedCost ==  quicksum(ShedCost[c]  for c in cont_list)   ,name='Eq_TotalShedCost')


    if FixedCost==True:
        model.addConstr(OF_MP ==  RDCost + TotalShedCost  ,name='Eq_OF')
        model.addConstr(GenCost <= (1+Alpha)*quicksum(Gen_data.loc[g]['b']*Pg_market[g] for g in G) ,name='Eq_GC2')
    else:
        model.addConstr(OF_MP == GenCost + RDCost + TotalShedCost  ,name='Eq_OF')


    model.setObjective(OF_MP,GRB.MINIMIZE)

    model.update()

    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs


    # Model and size statistics
    return {'model':model,
            'NumBinVars':NumBinVars,
            'NumConVars':NumConVars,
            'NumConstrs':NumConstrs
            }


def solve_hmmp_substation_mp(
    data,
    substation,
    MP2model,
    OSPmodel,
    FSPmodels,
    Zfixed=None,
    split_list=[],
    PgFix=None,
    QgFix=None,
    line_cont_list=[],
    Max_FSP_iter=20,
    FSP_criteria=0,
    Max_Sw_bus=0,
):
    """Solve one substation MP2 with local feasibility and optimality checks.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        substation: Substation whose topology is optimized.
        MP2model: Reusable MP2 created by ``create_hmmp_substation_mp``.
        OSPmodel: Reusable OSP for the selected substation.
        FSPmodels: Mapping from line contingencies to reusable FSP models.
        Zfixed: Current system topology used to fix other substations.
        split_list: Substations already selected for busbar splitting.
        PgFix: Active dispatch fixed by the central master problem.
        QgFix: Reactive dispatch fixed by the central master problem.
        line_cont_list: Line contingencies checked by the local FSP loop.
        Max_FSP_iter: Maximum number of local feasibility-cut iterations.
        FSP_criteria: Feasibility objective tolerance.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.

    Returns:
        A dictionary containing the topology, OSP duals and cost, shedding,
        solve time, objective improvement from splitting, and updated voltage
        and angle linearization points.
    """

    # Input data

    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    branch=data['branch']
    Lines=data['Lines']
    line=data['line']
    L2B=data['L2B']
    NumberL2B=data['NumberL2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']
    D2B=data['D2B']
    Gen_data=data['Gen_data']
    G=data['G']
    G2B=data['G2B']
    vmin = data['vmin']
    vmax = data['vmax']
    beta = data['beta']

    BigM_l=data['BigM_l']
    BigM_b=data['BigM_b']
    Maxdelta=data['Maxdelta']
    MaxV2=data['MaxV2']
    BigM_busbar=data['BigM_busbar']


    # Copy the reusable substation master model.
    model=MP2model.copy()


    cont_list = [0,substation,str(substation+'-1'),str(substation+'-2')]  #c0 is included!


    if Zfixed!=None:
        for b in Bus:
            if b!=substation:
                if b not in split_list: #everything assumed at default topology
                    model.addConstr(( model.getVarByName('z_bus['+str(b)+']') == 1 ), name='eqFix')
                    model.getVarByName('z_bus['+str(b)+']').Start = 1
                    for l,i,j in Lines.select('*',b,'*'):
                        model.addConstr((  model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']') == 0 ), name='eqFix')
                        model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']').Start = 0
                    for d,i in D2B.select('*',b):
                        model.addConstr((  model.getVarByName('z_d['+str(d)+']') == 0 ), name='eqFix')
                        model.getVarByName('z_d['+str(d)+']').Start = 0
                    for g,i in G2B.select('*',b):
                        model.addConstr((  model.getVarByName('z_g['+str(g)+']') == 0), name='eqFix')
                        model.getVarByName('z_g['+str(g)+']').Start = 0

                elif b in split_list:
                    model.addConstr(( model.getVarByName('z_bus['+str(b)+']') == (Zfixed['bus'][b]) ), name='eqFix')
                    model.getVarByName('z_bus['+str(b)+']').Start = (Zfixed['bus'][b])
                    for l,i,j in Lines.select('*',b,'*'):
                        model.addConstr((  model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']') == (Zfixed['l_i'][l,i,j])), name='eqFix')
                        model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']').Start = (Zfixed['l_i'][l,i,j])
                    for d,i in D2B.select('*',b):
                        model.addConstr((  model.getVarByName('z_d['+str(d)+']') == (Zfixed['d'][d]) ), name='eqFix')
                        model.getVarByName('z_d['+str(d)+']').Start = (Zfixed['d'][d])
                    for g,i in G2B.select('*',b):
                        model.addConstr((  model.getVarByName('z_g['+str(g)+']') == (Zfixed['g'][g]) ), name='eqFix')
                        model.getVarByName('z_g['+str(g)+']').Start = (Zfixed['g'][g])

            elif b==substation: #substation free to change
                model.getVarByName('z_bus['+str(b)+']').Start = 1
                for l,i,j in Lines.select('*',b,'*'):
                    model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']').Start = 0
                for d,i in D2B.select('*',b):
                    model.getVarByName('z_d['+str(d)+']').Start = 0
                for g,i in G2B.select('*',b):
                     model.getVarByName('z_g['+str(g)+']').Start = 0

    else:
        print('MP_2 requires a fixed topology from previuos iteration.')


    # Fixed central dispatch

    eq_Pg1=model.addConstrs( (   model.getVarByName('Pgi['+str(g)+',busbar1]')   == (1- model.getVarByName('z_g['+str(g)+']') )*PgFix[g]
                             for g in G ), name='eq_Pg1')

    eq_Pg2=model.addConstrs( (  model.getVarByName('Pgi['+str(g)+',busbar2]') == model.getVarByName('z_g['+str(g)+']')*PgFix[g]
                             for g in G  ), name='eq_Pg2')

    eq_Pg1=model.addConstrs( (   model.getVarByName('Qgi['+str(g)+',busbar1]')  == (1- model.getVarByName('z_g['+str(g)+']') )*QgFix[g]
                             for g in G ), name='eq_Qg1')

    eq_Pg2=model.addConstrs( (  model.getVarByName('Qgi['+str(g)+',busbar2]') == model.getVarByName('z_g['+str(g)+']')*QgFix[g]
                             for g in G  ), name='eq_Qg2')


    model._importantvars = [ model.getVarByName('z_bus['+str(substation)+']')  ] #z_bus[substation] is important for us
    model._data = [] # elements=[zbus,obj]

    model.update()

    start_time = time.time()
    model.optimize(callback=MP2_obj_cb)
    end_time = time.time()
    ex_time=end_time-start_time  #execution time

    tot_time = 0
    tot_time += ex_time


    status = model.Status

    if status == GRB.INFEASIBLE:
        print('\n\nMP2 was stopped with infeasibility!')
        print('cont list: ',cont_list)
        print(Zfixed)

        # Relax bounds to diagnose an infeasible model.
        print('\n\nThe model is infeasible; relaxing the bounds\n\n')
        orignumvars = model.NumVars
        # Relax only variable bounds.
        model.feasRelaxS(0, False, False, True) #relaxing constraints
        model.optimize()
        status = model.Status
        if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
                print('The relaxed model cannot be solved \
                       because it is infeasible or unbounded')
        if status != GRB.OPTIMAL:
            print('Optimization was stopped with status %d' % status)

        # Report artificial variables introduced by the relaxation.
        print('\nSlack values:')
        slacks = model.getVars()[orignumvars:]
        for sv in slacks:
            if sv.X > 1e-9:
                print('%s = %g' % (sv.VarName, sv.X))


    # MP2 solution and objective improvement
    start_obj = model._startobjval
    end_obj = model.objVal


    if model.getVarByName('z_bus['+str(substation)+']').x==0:
        for instance in model._data:
            if instance[0] == 1: #if z_bus==1
                if instance[1] +1e-3 < start_obj: #if the found objective for zbus=1 is lower than the MIPstart
                    start_obj = instance[1]


    delta_obj = start_obj - end_obj


    # Build the topology dictionary.
    def UpdateTopology():
        PgiMP={}; QgiMP={} ;z_busMP={}; z_gMP={}; z_dMP={}; z_liMP={}

        for g in G:
            for i in busbar:
                PgiMP[(g,i)]=model.getVarByName('Pgi['+str(g)+','+str(i)+']').x  #Pgi[g,i].x
                QgiMP[(g,i)]=model.getVarByName('Qgi['+str(g)+','+str(i)+']').x  #Qgi[g,i].x
        for b in Bus:
            z_busMP[b]=model.getVarByName('z_bus['+str(b)+']').x #z_bus[b].x
        for g in G:
            z_gMP[g]=model.getVarByName('z_g['+str(g)+']').x #z_g[g].x
        for d in DemandSet:
            z_dMP[d]=model.getVarByName('z_d['+str(d)+']').x #z_d[d].x
        for l,i,j in Lines:
            z_liMP[(l,i,j)]=model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']').x #z_li[l,i,j].x

        return {
            'Pgi':PgiMP,
            'Qgi':QgiMP,
            'bus':z_busMP,
              'g':z_gMP,
              'd':z_dMP,
              'l_i':z_liMP}

    TopologyMP=UpdateTopology()

    # Local line-feasibility loop
    k1=0        #FSP inner loop
    NumberofCuts=0

    while k1 < Max_FSP_iter:
        print('FSP check iter %i...'%k1,end='\r')
        FSP_Objc=0
        FSP_Obj=0

        for c in line_cont_list:
            if c in Find_lines_hops(data,sub=substation,hops=1):
                if c in data['line_cont_notradial']:

                    result = solve_hmmp_line_fsp(data=data,TopologyMP=TopologyMP,model0=FSPmodels[c])
                    FSP_Objc=result['OF_FSP']
                    FSP_Obj+=FSP_Objc
                    Mu=result['Mu']
                    tot_time += result['time']

                    if FSP_Objc > FSP_criteria:

                        model.addConstr(  ( FSP_Objc
                                    +quicksum( Mu['bus'][b]*( model.getVarByName('z_bus['+str(b)+']') -TopologyMP['bus'][b])  for b in Bus)
                                    +quicksum( Mu['g'][g]*( model.getVarByName('z_g['+str(g)+']') -TopologyMP['g'][g])  for g in G)
                                    +quicksum( Mu['d'][d]*( model.getVarByName('z_d['+str(d)+']') -TopologyMP['d'][d])  for d in DemandSet)
                                    +quicksum( Mu['l_i'][(l,i,j)]*( model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']') -TopologyMP['l_i'][l,i,j])  for l,i,j in Lines)
                                                <= 0   )
                                ,name=str('MP2-FSP-'+str(NumberofCuts)))
                        NumberofCuts+=1
        # FSP evaluation complete
        # Update MP2 with feasibility cuts.
        if FSP_Obj > FSP_criteria:
            k1+=1
            model.update()
            start_time = time.time()
            model.optimize()
            model.optimize(callback=MP2_obj_cb)
            end_time = time.time(); ex_time=end_time-start_time; tot_time += ex_time


            TopologyMP=UpdateTopology()

            start_obj = model._startobjval
            end_obj = model.objVal

            if model.getVarByName('z_bus['+str(substation)+']').x==0:
                for instance in model._data:
                    if instance[0] == 1: #if z_bus==1
                        if instance[1] < start_obj: #if the found objective for zbus=1 is lower than the MIPstart
                            start_obj = instance[1]
            delta_obj = start_obj - end_obj

            TopologyMP=UpdateTopology()


        else:
            break

            # End local FSP loop.
    # Local line-feasibility checks complete.

    Vol0_li, delta0_li = _extract_linearization_points(
        model,
        Lines,
        cont_list,
    )

    sub_cont_list = []
    for b in Bus:
        if b==substation:
            sub_cont_list+=[b,str(b+'-1'),str(b+'-2')]


    result = solve_hmmp_substation_osp(data=data,cont_list=sub_cont_list,TopologyMP=TopologyMP, model0=OSPmodel)
    OF_OSP=result['OF_OSP']
    Mu=result['Mu']
    tot_time += result['time']
    ShedCost_df = result['ShedCost_df']
    PdShed_dic = result['PdShed_dic']
    Vol0_li.update(result['Vol0_li'])
    delta0_li.update(result['delta0_li'])


    if model.getVarByName('z_bus['+str(substation)+']').x==0:
        print('\t\t\tbus %s openned.'%substation)


    # Results
    return {
        'TopologyMP':TopologyMP,
        'Mu':Mu,
        'OF_OSP':OF_OSP,
        'ShedCost_df':ShedCost_df,
        'PdShed_dic':PdShed_dic,
        'time':tot_time,
        'delta_obj':delta_obj,
        'Vol0_li':Vol0_li,
        'delta0_li':delta0_li
    }


def MP2_obj_cb(
    model,
    where,
):
    """Record MP2 incumbent objectives for open and closed busbar states.

    Args:
        model: Active MP2 Gurobi model with callback state initialized.
        where: Gurobi callback location code.

    Returns:
        None. Callback results are stored on the model object.
    """
    if where == GRB.Callback.MIPSOL:
        if model.cbGet(GRB.Callback.MIPSOL_SOLCNT) == 0:
            # creates new model attribute '_startobjval'
            model._startobjval = model.cbGet(GRB.Callback.MIPSOL_OBJ)

        zbus = model.cbGetSolution(model._importantvars )[0] #zbus value
        obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)

        model._data.append([zbus,obj])


def SC_SR_HMMP(
    data,
    line_cont_list=[],
    Max_iter=10,
    Max_OSP_iter=10,
    OSP_criteria=0.1,
    Max_FSP_iter=10,
    FSP_criteria=0,
    Pg_market_fix=None,
    FixedCost=False,
    Alpha=0,
    Pg_market=None,
    Probabilistic=False,
    Max_Sw_bus=0,
    Up_redispatch=0,
    Dn_redispatch=1.0,
    Vol0_li=None,
    delta0_li=None,
):
    """Run the proposed heuristic multi-MP method from Algorithm 1.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        line_cont_list: Line contingencies included in security evaluation.
        Max_iter: Maximum number of central multi-MP iterations.
        Max_OSP_iter: Reserved maximum OSP iteration count.
        OSP_criteria: Upper/lower-bound gap used as the stopping tolerance.
        Max_FSP_iter: Maximum number of feasibility-cut iterations.
        FSP_criteria: Feasibility objective tolerance.
        Pg_market_fix: Optional active dispatch fixed in the central master.
        FixedCost: Enforce the market-generation cost limit when ``True``.
        Alpha: Relative margin applied to the market-generation cost limit.
        Pg_market: Reference market dispatch for the fixed-cost constraint.
        Probabilistic: Weight shedding costs by contingency probabilities.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.
        Up_redispatch: Available upward active-power redispatch fraction.
        Dn_redispatch: Available downward active-power redispatch fraction.
        Vol0_li: Optional branch-voltage linearization points indexed by line
            endpoints and contingency; unit values are used when omitted.
        delta0_li: Optional branch-angle linearization points indexed by line
            endpoints and contingency; zero values are used when omitted.

    Returns:
        A dictionary containing the best topology and dispatch, convergence
        bounds, costs, shedding, timing, iteration details, model sizes, and
        the updated voltage and angle linearization points for that solution.
    """

    # Input data

    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    branch=data['branch']
    Lines=data['Lines']
    line=data['line']
    L2B=data['L2B']
    NumberL2B=data['NumberL2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']
    D2B=data['D2B']
    Gen_data=data['Gen_data']
    G=data['G']
    G2B=data['G2B']


    # Model
    model=gp.Model('BCC_main_SCOPF')
    model.Params.OutputFlag=0

    if Max_Sw_bus!=0:
        print('\nBe careful with Max switching!!\n')

    # Variables

    Pg=model.addVars(G,lb=0,vtype=GRB.CONTINUOUS,name='Pg')
    Qg=model.addVars(G,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg')

    GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost')
    OF_MP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_MP')

    Phi_b=model.addVars(Bus,lb=0,vtype=GRB.CONTINUOUS,name='Phi_b')

    line_cont_radial =[]
    for l in line_cont_list:
        if l not in data['line_cont_notradial']:
            line_cont_radial += [l]

    Phi_l = model.addVars(line_cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Phi_l')

    TotalPhi=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalPhi')

    # Constraints

    eqTotBalance=model.addConstr( (   quicksum(Pdemand.loc[d]['Pd'] for d in DemandSet) <= quicksum( Pg[g] for g in G) ),
                                 name='EqTotBalance')


    model.addConstrs((Pg[g]<= Gen_data.loc[g]['Pmax']  for g in G ),name='Eq_Pgmax')
    model.addConstrs(( -Pg[g] <= -Gen_data.loc[g]['Pmin']  for g in G),name='Eq_Pgmin')
    model.addConstrs((Qg[g]<= Gen_data.loc[g]['Qmax']  for g in G ),name='Eq_Qgmax')
    model.addConstrs(( -Qg[g] <= -Gen_data.loc[g]['Qmin']  for g in G),name='Eq_Qgmin')

    if Pg_market_fix is not None:
        model.addConstrs(( Pg[g] == Pg_market_fix[g]  for g in G),name='Eq_PFix')


    eqOF=model.addConstr(GenCost==quicksum(Gen_data.loc[g]['b']*Pg[g] for g in G ) ,name='Eq_OF')


    model.addConstr( TotalPhi == quicksum(Phi_b[b] for b in Bus) + quicksum(Phi_l[l] for l in line_cont_list)
                    ,name='Eq_Phi')


    if FixedCost==True:

        OF_MP =  TotalPhi + quicksum(0.001*Gen_data.loc[g]['b']*Qg[g]*Qg[g] for g in G )

        model.addConstr(GenCost <= (1+Alpha)*quicksum(Gen_data.loc[g]['b']*Pg_market[g] for g in G) ,name='Eq_GC2')

    else:
        OF_MP = GenCost + TotalPhi + quicksum(0.001*Gen_data.loc[g]['b']*Qg[g]*Qg[g] for g in G )


    model.setObjective(OF_MP,GRB.MINIMIZE)
    model.update()

    NumBinVarsMP0 = model.NumBinVars
    NumConVarsMP0 = model.NumVars - model.NumBinVars
    NumConstrsMP0 = model.NumConstrs

    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time

    tot_time = 0
    tot_time += ex_time

    time_iteration = {}
    time_iter_detail = {'MP':{},'SP':{}}
    MP_time=0; SP_time=0
    MP_time+=ex_time


    # Initial central master solve

    PgMP={}; QgMP={}
    for g in G:
        PgMP[g]=Pg[g].x
        QgMP[g]=Qg[g].x


    # Topology initialization
    FinalTopology={'bus':{b:1 for b in Bus},
                   'g':{g:0 for g in G}, 'd':{d:0 for d in DemandSet},
                   'l_i':{(l,i,j):0 for l,i,j in Lines},
                  'Pgi':{}, 'Qgi':{},
                  'Pg':{},'Qg':{}}
    for g in G:
        FinalTopology['Pgi'][(g,'busbar1')]=PgMP[g]; FinalTopology['Pgi'][(g,'busbar2')]=0
        FinalTopology['Qgi'][(g,'busbar1')]=QgMP[g]; FinalTopology['Qgi'][(g,'busbar2')]=0
        FinalTopology['Pg'][g]=PgMP[g]; FinalTopology['Qg'][g]=QgMP[g]

    TopologyOld=deepcopy(FinalTopology) #topology old is passed to MP2

    topology_default0 = {'bus':{b:1 for b in Bus},
                   'g':{g:0 for g in G}, 'd':{d:0 for d in DemandSet},
                   'l_i':{(l,i,j):0 for l,i,j in Lines},
                   'Pgi':{},'Qgi':{}}


    Mu={}
    UB_k={}; LB_k={}
    UB=1e8; LB=0


    k=0         #OSP outer loop
    k1=0        #FSP inner loop
    GenCost_k1={}
    FSP_Obj_k1={}
    Topology_k={}
    UB_min=UB

    NumberofCuts=0


    all_sub_cont_list=[]
    for b in Bus:
        all_sub_cont_list+=[b,str(b+'-1'),str(b+'-2')]

    all_cont=all_sub_cont_list+line_cont_list #used for reporting LoadShedding_df
    linearization_cont=[0]+all_cont
    if Vol0_li is None:
        Vol0_li={(l,i,j,c):1.0 for l,i,j in Lines for c in linearization_cont}
    if delta0_li is None:
        delta0_li={(l,i,j,c):0.0 for l,i,j in Lines for c in linearization_cont}

    ShedCost_df = pd.DataFrame(0,columns=['ShedCost(c)'],index=all_cont) #not including c0
    PdShed_dic = pd.DataFrame(0,columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,busbar,all_cont]))


    line_cont_notradial =[0] # in case of MP1=ED
    for c in line_cont_list:
        if c in data['line_cont_notradial']:
            line_cont_notradial +=[c]


    # Create reusable FSP, OSP, and MP2 models.
    print('creating FSP, OSP, MP2 models.')
    FSP_models={}
    OSP_models={}
    MP2_models={}

    NumBinVars = {'MPi':{},'SP':{}}; NumConVars = {'MPi':{},'SP':{}}; NumConstrs = {'MPi':{},'SP':{}}

    NumBinVars['MP0'] = NumBinVarsMP0
    NumConVars['MP0'] = NumConVarsMP0
    NumConstrs['Mp0'] = NumConstrsMP0

    for c in line_cont_notradial:
        res = create_hmmp_line_fsp(data=data,cont_list=[c],
                                  Vol0_li=Vol0_li,delta0_li=delta0_li)
        FSP_models[c] = res['model']
        NumBinVars['SP'][c] = res['NumBinVars']
        NumConVars['SP'][c] = res['NumConVars']
        NumConstrs['SP'][c] = res['NumConstrs']

    for c in line_cont_list:
        res = create_hmmp_substation_osp(data=data,cont_list=[c],Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,Max_Sw_bus=Max_Sw_bus,
                                        Vol0_li=Vol0_li,delta0_li=delta0_li)
        OSP_models[c] = res['model']
        NumBinVars['SP'][c] = res['NumBinVars']
        NumConVars['SP'][c] = res['NumConVars']
        NumConstrs['SP'][c] = res['NumConstrs']

    for sub in Bus:
        res = create_hmmp_substation_mp(data=data,substation=sub,Max_Sw_bus=Max_Sw_bus,Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,
                                             FixedCost=FixedCost,Alpha=Alpha,Pg_market=Pg_market,Probabilistic=Probabilistic,
                                             Vol0_li=Vol0_li,delta0_li=delta0_li)

        MP2_models[sub] = res['model']
        NumBinVars['MPi'][sub] = res['NumBinVars']
        NumConVars['MPi'][sub] = res['NumConVars']
        NumConstrs['MPi'][sub] = res['NumConstrs']


        res = create_hmmp_substation_osp(data=data,cont_list=[sub,str(sub+'-1'),str(sub+'-2')]
                                                ,Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,Max_Sw_bus=Max_Sw_bus,
                                                FixedCost=FixedCost,Alpha=Alpha,Pg_market=Pg_market,Probabilistic=Probabilistic,
                                                Vol0_li=Vol0_li,delta0_li=delta0_li)
        OSP_models[sub] = res['model']
        NumBinVars['SP'][sub] = res['NumBinVars']
        NumConVars['SP'][sub] = res['NumConVars']
        NumConstrs['SP'][sub] = res['NumConstrs']


    solution_dict = {} #iter:sol


    current_time = time.localtime()
    formatted_time = time.strftime("%Y-%m-%d %H:%M:%S", current_time)
    print("Current time and date:", formatted_time)

    # Algorithm 1 outer loop
    while k < Max_iter:


        print('\n**Main teration %3i'%(k))

        # Central non-radial line-feasibility loop
        k1=0
        while k1 < Max_FSP_iter:

            print('FSP iteration %3i'%k1)

            FSP_Objc=0
            FSP_Obj=0

            for g in G: #TopologyMP=def, so we check the FSP for the default top
                topology_default0['Pgi'][(g,'busbar1')]=PgMP[g]; topology_default0['Pgi'][(g,'busbar2')]=0
                topology_default0['Qgi'][(g,'busbar1')]=QgMP[g]; topology_default0['Qgi'][(g,'busbar2')]=0


            for c in (line_cont_notradial):
                result = solve_hmmp_line_fsp(data=data,TopologyMP=topology_default0,model0=FSP_models[c],print_result=False) #Problem: should we check default topology for disptach, or the new one?!
                FSP_Objc=result['OF_FSP']
                FSP_Obj+=FSP_Objc
                Lambda=result['Mu']
                tot_time += result['time']
                SP_time+=result['time']

                if FSP_Objc > FSP_criteria:

                    model.addConstr(  ( FSP_Objc
                                +quicksum( Lambda['Pg'][g]*(Pg[g]-PgMP[g])  for g in G)
                                +quicksum( Lambda['Qg'][g]*(Qg[g]-QgMP[g])  for g in G)
                                        <= 0   )
                        ,name=str('FSP-'+str(c)+'-'+str(NumberofCuts)))
                    NumberofCuts+=1

                # Update the central master with optimality cuts.
            if FSP_Obj > FSP_criteria:
                k1+=1
                model.update()
                start_time = time.time()
                model.optimize()
                end_time = time.time()
                ex_time=end_time-start_time
                tot_time += ex_time
                MP_time+=ex_time

                for g in G:
                    PgMP[g]=Pg[g].x
                    QgMP[g]=Qg[g].x
                    if FinalTopology['g'][g]==0:
                        FinalTopology['Pgi'][(g,'busbar1')]=PgMP[g]; FinalTopology['Pgi'][(g,'busbar2')]=0
                        FinalTopology['Qgi'][(g,'busbar1')]=QgMP[g]; FinalTopology['Qgi'][(g,'busbar2')]=0
                    elif FinalTopology['g'][g]==1:
                        FinalTopology['Pgi'][(g,'busbar1')]=0; FinalTopology['Pgi'][(g,'busbar2')]=PgMP[g]
                        FinalTopology['Qgi'][(g,'busbar1')]=0; FinalTopology['Qgi'][(g,'busbar2')]=QgMP[g]
                    FinalTopology['Pg'][g]=PgMP[g]
                    FinalTopology['Qg'][g]=QgMP[g]

                TopologyOld=deepcopy(FinalTopology)


            else:
                print('\n\t\tFSP converged in iteration %i'%k1)
                break

            # End central FSP loop.


        # Parallel substation MP2 and OSP evaluation
        OF_OSP=0 ; OF_OSP_b={}
        Mu={}
        MP2_res = {}
        selected_split = []

        for sub in (Bus):
            print('\tMP2 substation %3s'%sub)

            result=solve_hmmp_substation_mp(data=data,substation=sub,
                                        MP2model=MP2_models[sub],OSPmodel=OSP_models[sub],FSPmodels=FSP_models,
                            Zfixed=topology_default0,split_list=[],#all default
                            PgFix=PgMP,QgFix=QgMP,
                                        line_cont_list=line_cont_list,Max_FSP_iter=Max_FSP_iter,FSP_criteria=FSP_criteria)

            Vol0_li.update(result['Vol0_li'])
            delta0_li.update(result['delta0_li'])
            MP2_res[sub] = result
            OF_OSP_b[sub]=result['OF_OSP']
            OF_OSP+=OF_OSP_b[sub]
            Mu[sub]=result['Mu']
            topology=result['TopologyMP']
            tot_time += result['time']
            MP_time += result['time']

            for cont in all_cont:
                if cont in result['ShedCost_df'].index:
                    ShedCost_df.loc[cont] = result['ShedCost_df'].loc[cont]
                    for d in DemandSet:
                        for i in busbar:
                            PdShed_dic.loc[(d,i,cont)] = result['PdShed_dic'].loc[(d,i,cont)]

            # Update the combined topology.
            FinalTopology['bus'][sub]=topology['bus'][sub]
            for g,i in G2B.select('*',sub):
                FinalTopology['g'][g]=topology['g'][g]
                if topology['g'][g]==0:
                    FinalTopology['Pgi'][(g,'busbar1')]=PgMP[g]; FinalTopology['Pgi'][(g,'busbar2')]=0
                    FinalTopology['Qgi'][(g,'busbar1')]=QgMP[g]; FinalTopology['Qgi'][(g,'busbar2')]=0
                elif topology['g'][g]==1:
                    FinalTopology['Pgi'][(g,'busbar1')]=0; FinalTopology['Pgi'][(g,'busbar2')]=PgMP[g]
                    FinalTopology['Qgi'][(g,'busbar1')]=0; FinalTopology['Qgi'][(g,'busbar2')]=QgMP[g]
            for d,i in D2B.select('*',sub):
                FinalTopology['d'][d]=topology['d'][d]
            for l,i,j in Lines.select('*',sub,'*'):
                FinalTopology['l_i'][l,i,j]=topology['l_i'][l,i,j]

        # Substation MP2 evaluations complete.


        for counter in range(1,Max_Sw_bus+1): # at counter c, only c+1 split is aloowed. so busbars are opened one by one

            split_candidates = []
            for b in Bus:
                if b not in selected_split: #selected splits are already opnned. the rest should be checked
                    if FinalTopology['bus'][b]==0:
                        split_candidates.append(b)


            if (len(split_candidates+selected_split) <= counter): #not violatied.
                print('Zmax not violated. Continue.')
                break

            else:
                print(f'Zmax violated or more than 1 split/iteration is fonud. {len(split_candidates)} buses are candidates to split:  ',split_candidates )

                best_delta=0
                for b in split_candidates:
                    if MP2_res[b]['delta_obj']>best_delta:
                        best_delta = MP2_res[b]['delta_obj']
                        best_sub = b
                print('bus %s is split.'%(best_sub))
                selected_split += [best_sub]

                flag_list = [b for b in split_candidates if b not in selected_split]
                print('flag list: ',flag_list)

                for sub in (flag_list):

                    result=solve_hmmp_substation_mp(data=data,substation=sub,
                                                MP2model=MP2_models[sub],OSPmodel=OSP_models[sub],FSPmodels=FSP_models,
                                    Zfixed=FinalTopology,split_list=selected_split, #all connected to default, split bus is final_top
                                    PgFix=PgMP,QgFix=QgMP,
                                                line_cont_list=line_cont_list,Max_FSP_iter=Max_FSP_iter,FSP_criteria=FSP_criteria,Max_Sw_bus=Max_Sw_bus)

                    Vol0_li.update(result['Vol0_li'])
                    delta0_li.update(result['delta0_li'])
                    MP2_res[sub] = result
                    OF_OSP_b[sub]=result['OF_OSP']
                    Mu[sub]=result['Mu']
                    topology=result['TopologyMP']
                    tot_time += result['time']
                    MP_time+=result['time']

                    for cont in all_cont:
                        if cont in result['ShedCost_df'].index:
                            ShedCost_df.loc[cont] = result['ShedCost_df'].loc[cont]
                            for d in DemandSet:
                                for i in busbar:
                                    PdShed_dic.loc[(d,i,cont)] = result['PdShed_dic'].loc[(d,i,cont)]

                    # Update the combined topology.
                    FinalTopology['bus'][sub]=topology['bus'][sub]
                    for g,i in G2B.select('*',sub):
                        FinalTopology['g'][g]=topology['g'][g]
                        if topology['g'][g]==0:
                            FinalTopology['Pgi'][(g,'busbar1')]=PgMP[g]; FinalTopology['Pgi'][(g,'busbar2')]=0
                            FinalTopology['Qgi'][(g,'busbar1')]=QgMP[g]; FinalTopology['Qgi'][(g,'busbar2')]=0
                        elif topology['g'][g]==1:
                            FinalTopology['Pgi'][(g,'busbar1')]=0; FinalTopology['Pgi'][(g,'busbar2')]=PgMP[g]
                            FinalTopology['Qgi'][(g,'busbar1')]=0; FinalTopology['Qgi'][(g,'busbar2')]=QgMP[g]
                    for d,i in D2B.select('*',sub):
                        FinalTopology['d'][d]=topology['d'][d]
                    for l,i,j in Lines.select('*',sub,'*'):
                        FinalTopology['l_i'][l,i,j]=topology['l_i'][l,i,j]

                OF_OSP = 0
                for sub in Bus:
                    OF_OSP+=OF_OSP_b[sub]


        # Busbar-split selection complete.


        # All substation MP2 and OSP models are solved.


        # Radial-line contingency OSP evaluation

        # Accumulate line and substation recourse costs.
        OF_OSP_c={}

        for c in line_cont_list:
            result=solve_hmmp_substation_osp(data=data,model0=OSP_models[c], cont_list=[c],TopologyMP=FinalTopology)

            Vol0_li.update(result['Vol0_li'])
            delta0_li.update(result['delta0_li'])

            OF_OSP_c[c]=result['OF_OSP']
            OF_OSP+=OF_OSP_c[c]          #OF_OSP has values before here
            tot_time += result['time']
            SP_time+=result['time']

            Mu_Pg={}; Mu_Qg={}
            for g in G:
                if FinalTopology['g'][g] == 0:
                    Mu_Pg[g]=result['Mu']['Pgi'][g,'busbar1']; Mu_Qg[g]=result['Mu']['Qgi'][g,'busbar1']
                elif FinalTopology['g'][g] == 1:
                    Mu_Pg[g]=result['Mu']['Pgi'][g,'busbar2']; Mu_Qg[g]=result['Mu']['Qgi'][g,'busbar2']

            Mu[c]={'Pg':Mu_Pg ,'Qg':Mu_Qg} #Mu has values before here


            for cont in all_cont:
                if cont in result['ShedCost_df'].index:
                    ShedCost_df.loc[cont] = result['ShedCost_df'].loc[cont]
                    for d in DemandSet:
                        for i in busbar:
                            PdShed_dic.loc[(d,i,cont)] = result['PdShed_dic'].loc[(d,i,cont)]


        obj = model.getObjective().getValue() #this is OF_MP.x

        LB=obj ; LB_k[k]=LB
        if FixedCost==True:
            UB=OF_OSP ; UB_k[k]=UB
        else:
            UB=GenCost.x+OF_OSP ; UB_k[k]=UB
        time_iteration[k] = tot_time

        time_iter_detail['MP'][k]=MP_time; MP_time=0
        time_iter_detail['SP'][k]=SP_time; SP_time=0


        k_min = min(UB_k, key=UB_k.get)
        UB_min = UB_k[k_min]


        solution_dict[k] = {'FinalTopology':FinalTopology,'Pg':PgMP,'Qg':QgMP,
                            'GenCost':GenCost.x,'TotalShedCost':OF_OSP,'ShedCost_df':ShedCost_df,'PdShed_dic':PdShed_dic , 'Cost_tot':UB,
                            'time':tot_time,'time_iteration':time_iteration,
                            'Vol0_li':deepcopy(Vol0_li),'delta0_li':deepcopy(delta0_li)}


        print('********UB_k=%.2f and LB=%.2f  =>  Gap=%.3f'%(UB,LB,UB-LB))
        print('UB_min: %.2f'%UB_min)


        if Pg_market_fix is not None:
            print('\nPg fixed, only one iteration! UB is %.2f'%(UB))
            time_iteration[k] = tot_time
            break


        if (UB-LB) >= OSP_criteria:

            for sub in Bus:  #multi heuristic cuts
                model.addConstr(  ( OF_OSP_b[sub]
                            +quicksum( Mu[sub]['Pg'][g]*(Pg[g]-PgMP[g])  for g in G)
                            +quicksum( Mu[sub]['Qg'][g]*(Qg[g]-QgMP[g])  for g in G)
                                <= Phi_b[sub]  )
                ,name=str('OSP-'+sub+'-'+str(NumberofCuts)))
                NumberofCuts+=1

            for c in line_cont_list:
                model.addConstr(  ( OF_OSP_c[c]
                                   +quicksum( Mu[c]['Pg'][g]*(Pg[g]-PgMP[g])  for g in G)
                                   +quicksum( Mu[c]['Qg'][g]*(Qg[g]-QgMP[g])  for g in G)
                                        <= Phi_l[c]  )
                        ,name=str('OSP-'+c+'-'+str(NumberofCuts)))
                NumberofCuts+=1


            # Update the central master with optimality cuts.

            k+=1
            model.update()
            start_time = time.time()
            model.optimize()
            end_time = time.time()
            ex_time=end_time-start_time
            tot_time += ex_time
            MP_time += ex_time


            # Update central dispatch and topology.
            for g in G:
                PgMP[g]=Pg[g].x
                QgMP[g]=Qg[g].x
            for g in G:
                if FinalTopology['g'][g]==0:
                    FinalTopology['Pgi'][(g,'busbar1')]=PgMP[g]
                    FinalTopology['Pgi'][(g,'busbar2')]=0
                    FinalTopology['Qgi'][(g,'busbar1')]=QgMP[g]
                    FinalTopology['Qgi'][(g,'busbar2')]=0
                elif FinalTopology['g'][g]==1:
                    FinalTopology['Pgi'][(g,'busbar1')]=0
                    FinalTopology['Pgi'][(g,'busbar2')]=PgMP[g]
                    FinalTopology['Qgi'][(g,'busbar1')]=0
                    FinalTopology['Qgi'][(g,'busbar2')]=QgMP[g]
                FinalTopology['Pg'][g]=PgMP[g]
                FinalTopology['Qg'][g]=QgMP[g]

            TopologyOld=deepcopy(FinalTopology)


        elif (UB-LB)<OSP_criteria:
            print('\nMain Gap is %.2f, MP2 loop break!'%(UB-LB))
            time_iteration[k] = tot_time
            break #while main iteration


    print('Best solution is UB_kmin=%.2f and LB=%.2f  =>  Gap=%.3f'%(UB_k[k_min],LB,UB_k[k_min]-LB))

    for b in Bus:
        if solution_dict[k_min]['FinalTopology']['bus'][b]==0:
            print('bus %s is open.'%b)


    # Results
    return { 'TopologyDict':solution_dict[k_min]['FinalTopology'],
            'Pg':solution_dict[k_min]['Pg'],
            'Qg':solution_dict[k_min]['Qg'],
            'UB_k':UB_k,
            'LB_k':LB_k,
            'Cost_tot':solution_dict[k_min]['Cost_tot'],
        'GenCost':solution_dict[k_min]['GenCost'],
        'TotalShedCost':solution_dict[k_min]['TotalShedCost'],
        'ShedCost_df':solution_dict[k_min]['ShedCost_df'],
        'PdShed_dic':solution_dict[k_min]['PdShed_dic'],
             'time':tot_time,
            'time_iteration':time_iteration,
            'time_iter_detail':time_iter_detail,
            'solution_dict':solution_dict,
'NumBinVars':NumBinVars,
'NumConVars':NumConVars,
'NumConstrs':NumConstrs,
'Vol0_li':solution_dict[k_min]['Vol0_li'],
'delta0_li':solution_dict[k_min]['delta0_li']
    }


