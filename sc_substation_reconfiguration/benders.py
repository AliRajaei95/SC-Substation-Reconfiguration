"""Benders-decomposition baselines from the Section IV-A case study.

This script includes two models: BD-C, which uses classical Benders cuts, and
BD-H, which uses heuristic cuts.

Select the cut strategy through the ``heuristic_cut`` input of
``BCC_Benders_Classic_AC_v62``: set it to ``False`` for BD-C or ``True`` for
BD-H.

The master problem selects topology and normal dispatch, while feasibility and
optimality subproblems assess line, busbar, and coupler contingencies.
"""


import pandas as pd
import numpy as np
import itertools
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
import time
from tqdm import tqdm

from .original_MIP_model import *  # Shared data preparation and AC-OPF utilities.


# Originally tested with Python 3.7 and Gurobi 9.0.


def create_FSP_line_v62(
    data,
    cont_list=None,
):
    """Build the line-contingency feasibility subproblem (FSP).

    The FSP holds topology and dispatch linking variables available for a later
    solve and minimizes feasibility slack. It does not permit load shedding or
    corrective redispatch.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        cont_list: Line-contingency identifiers represented in the model.

    Returns:
        A dictionary containing the unsolved Gurobi ``model`` and its binary,
        continuous, and constraint counts.
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

    
    # Feasibility subproblem model
    model=gp.Model('BCC Benders FSP-1')
    
    cont_list=cont_list.copy()   
        
   
    model.Params.OutputFlag=0

  
    # Variables
    Pgi=model.addVars(G,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pgi')
    Qgi=model.addVars(G,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qgi') 
    Pg=model.addVars(G,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pg')
    Qg=model.addVars(G,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg')


    Pdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pdi')
    Qdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Qdi')


    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')         
    Pflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_li')
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')         
    Qflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_li')
    
    Pflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_bus')  #From b1 to b2!
    Qflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_bus')  #From b1 to b2!

    
    delta_bi=model.addVars(Bus,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi') 
    delta_li=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_le') 
    V2_bi=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi') 
    V2_li=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_le') 


    z_bus=model.addVars(Bus,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_bus')    
    z_li=model.addVars(Lines,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_li')  
    z_g=model.addVars(G,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_g') 
    z_d=model.addVars(DemandSet,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_d') 

    OF_FSP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_FSP')             # obj func var
    
    sp_up=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sp_up')
    sp_dn=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sp_dn')
    sq_up=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sq_up')
    sq_dn=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sq_dn')
     

    # Generator constraints
    eq_Pg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar1']  
                             for g in G   ), name='eq_Pg1min')
    eq_Pg1max=model.addConstrs( (    Pgi[g,'busbar1'] <=(1-z_g[g])*Gen_data.loc[g]['Pmax']  
                             for g in G  ), name='eq_Pg1max')
    
    eq_Pg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar2']   
                             for g in G   ), name='eq_Pg2min')
    eq_Pg2max=model.addConstrs( (   Pgi[g,'busbar2']  <= z_g[g]*Gen_data.loc[g]['Pmax']   
                             for g in G   ), name='eq_Pg2max')
    
    model.addConstrs( (   Pg[g] == Pgi[g,'busbar1']+Pgi[g,'busbar2']   for g in G   ), name='eq_Pg')
    model.addConstrs( (   Qg[g] == Qgi[g,'busbar1']+Qgi[g,'busbar2']   for g in G   ), name='eq_Qg')

    
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
    -branch.loc[(l,i,j)]['b_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) #+ Ploss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')
    
    model.addConstrs( (  Qflow[l,i,j,c] == 
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j,c] - V2_li[l,j,i,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) 
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')
    
    # Branch losses


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


    # Busbar and coupler contingencies
    model.addConstrs((  Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])     
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmax')
    model.addConstrs((  -Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])     
                          for b in Bus for ll in cont_list if b!=ll), name='eq_busbarPflowmin')
    
    
    # Balance 
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
    
    
    # Feasibility objective
    model.addConstr(OF_FSP ==  quicksum( sp_up[b,m,c]+sp_dn[b,m,c]+sq_up[b,m,c]+sq_dn[b,m,c]    
                                        for b in Bus for m in busbar for c in cont_list)  ,name='Eq_OF_FSP')


    model.setObjective(OF_FSP,GRB.MINIMIZE)

    model.update()
    
        
    # Model statistics
    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs
    
        
    return {'model':model,
            'NumBinVars':NumBinVars,
            'NumConVars':NumConVars,
            'NumConstrs':NumConstrs
            }
    
    
def solve_FSP_line_v62(
    data,
    TopologyMP,
    model0,
    print_result=False,
):
    """Solve a prepared FSP for a topology proposed by the master problem.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        TopologyMP: Master-problem topology and dispatch values to impose on
            the subproblem.
        model0: Unsolved FSP model created by ``create_FSP_line_v62``.
        print_result: Print the feasibility objective when ``True``.

    Returns:
        A dictionary containing linking-constraint duals (``Mu``), solve time,
        and the FSP objective value (``OF_FSP``).
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
     

    model.addConstrs(  (model.getVarByName('Pg['+str(g)+']') == TopologyMP['Pg'][g]     for g in G )  ,name='eqBendersPg')
    
    model.addConstrs(  (model.getVarByName('Qg['+str(g)+']') == TopologyMP['Qg'][g]     for g in G  )  ,name='eqBendersQg')

    
    model.addConstrs(  (model.getVarByName('z_bus['+str(b)+']') == (TopologyMP['bus'][b])      for b in Bus )  ,name='eqBendersZ_bus')
        
    model.addConstrs(  (model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')  == (TopologyMP['l_i'][l,i,j])      for l,i,j in Lines )  ,name='eqBendersZ_li')
    
    model.addConstrs(  (model.getVarByName('z_g['+str(g)+']') == (TopologyMP['g'][g])      for g in G )  ,name='eqBendersZ_g')
    
    model.addConstrs(  (model.getVarByName('z_d['+str(d)+']') == (TopologyMP['d'][d])     for d in DemandSet )  ,name='eqBendersZ_d')
    
    
    model.update()
    # Solve the feasibility subproblem.
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time  #execution time
    
    
    status = model.Status
        
    
    if status == GRB.INFEASIBLE:
        print('\n\nFSP Optimization was stopped with infeasibility!')
        
        
        # Relax the bounds and try to make the model feasible
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

        # Report artificial variables introduced by the feasibility relaxation.
        print('\nSlack values:')
        slacks = model.getVars()[orignumvars:]
        for sv in slacks:
            if sv.X > 1e-9:
                print('%s = %g' % (sv.VarName, sv.X))
    
    
    # Extract dual multipliers for Benders cuts.
    MuPgi={}; MuQgi={}; Muz_bus={}; Muz_li={}; Muz_g={}; Muz_d={}; MuPg={}; MuQg={}
       
    for g in G:
        MuPg[g] = model.getConstrByName(str('eqBendersPg['+g+']')).pi
        MuQg[g] = model.getConstrByName(str('eqBendersQg['+g+']')).pi
    for b in Bus:
        Muz_bus[b]=model.getConstrByName(str('eqBendersZ_bus['+b+']')).pi
 
    for g in G:
        Muz_g[g]=model.getConstrByName(str('eqBendersZ_g['+g+']')).pi
    for d in DemandSet:
        Muz_d[d]=model.getConstrByName(str('eqBendersZ_d['+d+']')).pi
    for l,i,j in Lines:
        Muz_li[(l,i,j)]=model.getConstrByName(str('eqBendersZ_li['+l+','+i+','+j+']')).pi
    
    Mu={
        'Pg':MuPg,
        'Qg':MuQg,
        'bus':Muz_bus,
        'g':Muz_g,
        'd':Muz_d,
        'l_i':Muz_li
    }


    if print_result==True:
        
        print(model.getVarByName('OF_FSP').x)
        
    
    # Results
    return { 
        'Mu':Mu,
        'time':ex_time,
        'OF_FSP': model.getVarByName('OF_FSP').x         
    }
    
    
def create_OSP_substation_v62(
    data,
    cont_list=None,
    Up_redispatch=0,
    Dn_redispatch=1.0,
    Max_Sw_bus=0,
):
    """Build the contingency optimality subproblem (OSP).

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        cont_list: Busbar, coupler, or combined contingency identifiers to
            represent in the model.
        Up_redispatch: Fraction of generator capacity available for upward
            active-power redispatch.
        Dn_redispatch: Fraction of generator capacity available for downward
            active-power redispatch.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.

    Returns:
        A dictionary containing the unsolved Gurobi ``model`` and its binary,
        continuous, and constraint counts.
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

    
    # Optimality subproblem model
    model=gp.Model('BCC Benders Optimality-SP')
    
    cont_list=cont_list.copy()   
    model.Params.OutputFlag=0


    # Reserve cost
    res2prod_coeff=0.5
    Gen_data['c_res']=Gen_data['b']*res2prod_coeff
    
    
    # Variables
    Pgi=model.addVars(G,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pgi')  
    dPgi_up=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_up')
    dPgi_dn=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_dn')
    Qgi=model.addVars(G,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qgi')  
    dQgi_up=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_up')
    dQgi_dn=model.addVars(G,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_dn')
    Pg=model.addVars(G,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pg')
    Qg=model.addVars(G,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg')

    Pdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pdi')
    Pdi_Shed=model.addVars(DemandSet,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Pdi_Shed')
    Qdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Qdi')
    Qdi_Shed=model.addVars(DemandSet,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Qdi_Shed')

    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')         
    Pflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_li')
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')         
    Qflow_li=model.addVars(Lines,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_li')
    
    Pflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_bus')  #From b1 to b2!
    Qflow_bus=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_bus')  #From b1 to b2!
    
    delta_bi=model.addVars(Bus,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi') 
    delta_li=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_le') 
    V2_bi=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi') 
    V2_li=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_le') 


    z_bus=model.addVars(Bus,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_bus')    
    z_li=model.addVars(Lines,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_li')  
    z_g=model.addVars(G,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_g') 
    z_d=model.addVars(DemandSet,lb=0,ub=1,vtype=GRB.CONTINUOUS,name='z_d')  

    OF_OSP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_OSP')             # obj func var
    RDCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='RDCost')
    ShedCost=model.addVars(cont_list,lb=0,vtype=GRB.CONTINUOUS,name='ShedCost')
    TotalShedCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalShedCost')

 
    # Reliability switching 
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
    
    eq_MaxSw_bus=model.addConstr((  quicksum( (1-z_bus[b]) for b in Bus) <= Max_Sw_bus ) , name='eq_MaxSw_bus')        # Eq max line sw 
    
    
    # Symmetry breaking
    for b in Bus:
        lmin,b=L2B.select('*',b)[0]
        lmin,i,j=Lines.select(lmin,b,'*')[0]
        eq_symmetry1=model.addConstr((  z_li[lmin,i,j] ==0  ), name='eq_symmetry')

    eq_zbus_ref=model.addConstr( (  z_bus[Bus[0]]==1    ), name='eq_zbus_ref')

    
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
    
    
    model.addConstrs( (   Pg[g] == Pgi[g,'busbar1']+Pgi[g,'busbar2']   for g in G   ), name='eq_Pg')
    model.addConstrs( (   Qg[g] == Qgi[g,'busbar1']+Qgi[g,'busbar2']   for g in G   ), name='eq_Qg')

    
    Eq_dPgUp_res=model.addConstrs(( dPgi_up[g,i,c] <= Gen_data.loc[g]['Pmax']*Up_redispatch  
                          for g in G for i in busbar for c in cont_list ),name='Eq_dPgUp_res')
    Eq_dPgDn_res=model.addConstrs(( dPgi_dn[g,i,c] <= Gen_data.loc[g]['Pmax']*Dn_redispatch  
                            for g in G for i in busbar for c in cont_list ),name='Eq_dPgDn_res')
    model.addConstrs(( dPgi_up[g,i,c] == 0  
                          for g in G for i in busbar for c in cont_list if c==0 ),name='Eq_dPg0')
    model.addConstrs(( dPgi_dn[g,i,c] == 0  
                            for g in G for i in busbar for c in cont_list if c==0 ),name='Eq_dPgDn_res')
    
    Eq_dQgUp_res=model.addConstrs(( dQgi_up[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #*Up_redispatch  
                          for g in G for i in busbar for c in cont_list ),name='Eq_dQgUp_res')
    Eq_dQgDn_res=model.addConstrs(( dQgi_dn[g,i,c] <= Gen_data.loc[g]['Qmax']*1# *Dn_redispatch  
                            for g in G for i in busbar for c in cont_list ),name='Eq_dQgDn_res')
    model.addConstrs(( dQgi_up[g,i,c] == 0  
                          for g in G for i in busbar for c in cont_list if c==0 ),name='Eq_dQg0')
    model.addConstrs(( dQgi_dn[g,i,c] == 0  
                            for g in G for i in busbar for c in cont_list if c==0 ),name='Eq_dQgDn_res')

    
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
    model.addConstrs( (  Pdi_Shed[d,i,c]==0 
                              for d in DemandSet for i in busbar 
                               for c in cont_list
                      if c==0) , name='eq_PdShedlimit0')
    
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
    -branch.loc[(l,i,j)]['b_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) #+ Ploss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')
    
    model.addConstrs( (  Qflow[l,i,j,c] == 
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j,c] - V2_li[l,j,i,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) 
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')
    
    # Branch losses

 
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
    

    # Balance 
    
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
    
    model.addConstrs( (ShedCost[c] == quicksum(Pdemand.loc[d]['ShedCost']*Pdi_Shed[d,i,c] 
                                              for d in DemandSet for i in busbar)     for c in cont_list) ,name='Eq_ShedC')
    
    # Redispatch and load-shedding objective
    model.addConstr(TotalShedCost ==  quicksum(ShedCost[c]  for c in cont_list)   ,name='Eq_TotalShedCost')
    
    model.addConstr(OF_OSP == RDCost + TotalShedCost ,name='Eq_OF_OSP')


    model.setObjective(OF_OSP,GRB.MINIMIZE)

    model.update()
    
    
    # Model statistics
    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs
    
        
    return {'model':model,
            'NumBinVars':NumBinVars,
            'NumConVars':NumConVars,
            'NumConstrs':NumConstrs
            }
    
    
def solve_OSP_substation_v62(
    data,
    TopologyMP,
    model0,
    cont_list=None,
):
    """Solve a prepared OSP for a topology proposed by the master problem.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        TopologyMP: Master-problem topology and dispatch values to impose on
            the subproblem.
        model0: Unsolved OSP model created by
            ``create_OSP_substation_v62``.
        cont_list: Contingency identifiers represented by ``model0``.

    Returns:
        A dictionary containing linking-constraint duals, solve time, the OSP
        objective, contingency shedding costs, and demand-level shedding.
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
    model.Params.OutputFlag=0


    model.addConstrs(  (model.getVarByName('Pg['+str(g)+']') == TopologyMP['Pg'][g]     for g in G )  ,name='eqBendersPg')
    
    model.addConstrs(  (model.getVarByName('Qg['+str(g)+']') == TopologyMP['Qg'][g]     for g in G  )  ,name='eqBendersQg')
    
    model.addConstrs(  (model.getVarByName('z_bus['+str(b)+']') == (TopologyMP['bus'][b])      for b in Bus )  ,name='eqBendersZ_bus')
        
    model.addConstrs(  (model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')  == (TopologyMP['l_i'][l,i,j])      for l,i,j in Lines )  ,name='eqBendersZ_li')
    
    model.addConstrs(  (model.getVarByName('z_g['+str(g)+']') == (TopologyMP['g'][g])      for g in G )  ,name='eqBendersZ_g')
    
    model.addConstrs(  (model.getVarByName('z_d['+str(d)+']') == (TopologyMP['d'][d])     for d in DemandSet )  ,name='eqBendersZ_d')
    

    model.update()
    # Solve the optimality subproblem.
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time 
    

    # Extract dual multipliers for Benders cuts.
    MuPgi={}; MuQgi={}; Muz_bus={}; Muz_li={}; Muz_g={}; Muz_d={}; MuPg={}; MuQg={}
       
    for g in G:
        MuPg[g] = model.getConstrByName(str('eqBendersPg['+g+']')).pi
        MuQg[g] = model.getConstrByName(str('eqBendersQg['+g+']')).pi
    
    for b in Bus:
        Muz_bus[b]=model.getConstrByName(str('eqBendersZ_bus['+b+']')).pi  
    for g in G:
        Muz_g[g]=model.getConstrByName(str('eqBendersZ_g['+g+']')).pi
    for d in DemandSet:
        Muz_d[d]=model.getConstrByName(str('eqBendersZ_d['+d+']')).pi
    for l,i,j in Lines:
        Muz_li[(l,i,j)]=model.getConstrByName(str('eqBendersZ_li['+l+','+i+','+j+']')).pi
    
    Mu={
        'Pg':MuPg,
        'Qg':MuQg,
        'bus':Muz_bus,
        'g':Muz_g,
        'd':Muz_d,
        'l_i':Muz_li
    }
    
                
    # Collect contingency and demand-level load shedding.
    ShedCost_df = pd.DataFrame(columns=['ShedCost(c)'],index=cont_list, dtype=float)
    PdShed_dic = pd.DataFrame(columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,busbar,cont_list]), dtype=float)
    for c in cont_list:
        ShedCost_df.loc[c]=model.getVarByName('ShedCost['+str(c)+']').x 
        for d in DemandSet:
            for i in busbar:
                PdShed_dic.loc[d,i,c]=model.getVarByName('Pdi_Shed['+str(d)+','+str(i)+','+str(c)+']').x
    
    
    # Results
    return { 
        'Mu':Mu,
        'time':ex_time,
        'OF_OSP': model.getVarByName('OF_OSP').x, 
        'ShedCost_df':ShedCost_df,
        'PdShed_dic':PdShed_dic,
    }
    
    

def BCC_Benders_Classic_AC_v62(
    data,
    line_cont_list=[],
    Max_FSP_iter=10,
    FSP_criteria=0,
    FSP_single=True,
    Max_OSP_iter=10,
    OSP_criteria=0.1,
    substation_OSP=True,
    heuristic_cut=True,
    hops_cut=0,
    all_cont_OSP=False,
    Max_Sw_bus=0,
    Up_redispatch=0,
    Dn_redispatch=1.0,
    Zfixdict=None,
    PQgFix=None,
    PgFix=None,
    print_result=False,
):
    """Solve the BD-C or BD-H security-constrained reconfiguration model.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        line_cont_list: Line contingencies checked by the feasibility
            subproblem.
        Max_FSP_iter: Maximum number of FSP cut iterations per outer iteration.
        FSP_criteria: FSP objective tolerance used to declare feasibility.
        FSP_single: Solve line contingencies individually when ``True``.
        Max_OSP_iter: Maximum number of outer OSP iterations.
        OSP_criteria: Upper/lower-bound gap used as the OSP stopping tolerance.
        substation_OSP: Solve one OSP per substation when ``True``.
        heuristic_cut: Select BD-H heuristic cuts when ``True``; use BD-C
            classical Benders cuts when ``False``.
        hops_cut: Network-neighborhood depth used to localize heuristic cuts;
            zero uses only the affected substation.
        all_cont_OSP: Solve OSPs by individual contingency instead of by
            substation when ``True``.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.
        Up_redispatch: Fraction of generator capacity available for upward
            active-power redispatch.
        Dn_redispatch: Fraction of generator capacity available for downward
            active-power redispatch.
        Zfixdict: Optional dictionary of topology variables to fix.
        PQgFix: Optional dictionary of busbar-level active and reactive
            generation values to fix.
        PgFix: Optional dictionary of total active-generation values to fix.
        print_result: Print iteration and solution details when ``True``.

    Returns:
        A dictionary containing the selected topology, dispatch, convergence
        bounds, costs, load shedding, timings, and model-size statistics.
    """
    
    
    # Input data
    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
    Max_timelimit=data['Max_timelimit']=900 #900
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

    
    # Master problem
    model=gp.Model('Benders_MPnormal_FSP1_OSP2_v62')
    
    line_cont_list=line_cont_list.copy()
    
    
    model.Params.MIPGap=Max_MIPGap
    model.Params.timelimit=Max_timelimit
    model.Params.OutputFlag=0

    
    Max_Sw_bus=Max_Sw_bus  #Max number of bus splitting switching
   

    # Master-problem variables
    Pgi=model.addVars(G,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pgi')  
    Qgi=model.addVars(G,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qgi')  
    Pg=model.addVars(G,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pg')
    Qg=model.addVars(G,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg')
    

    Pdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pdi')
    Qdi=model.addVars(DemandSet,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Qdi')

    Pflow=model.addVars(Lines,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')         
    Pflow_li=model.addVars(Lines,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_li')
    Qflow=model.addVars(Lines,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')         
    Qflow_li=model.addVars(Lines,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_li')
    
    
    Pflow_bus=model.addVars(Bus,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_bus')  #From b1 to b2!
    Qflow_bus=model.addVars(Bus,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_bus')  #From b1 to b2!

    delta_bi=model.addVars(Bus,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi') 
    delta_li=model.addVars(Lines,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_le') 
    V2_bi=model.addVars(Bus,busbar,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi') 
    V2_li=model.addVars(Lines,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_le') 
    
   
    z_bus=model.addVars(Bus,vtype=GRB.BINARY,name='z_bus')    
    z_li=model.addVars(Lines,vtype=GRB.BINARY,name='z_line_i')  
    z_g=model.addVars(G,vtype=GRB.BINARY,name='z_g') 
    z_d=model.addVars(DemandSet,vtype=GRB.BINARY,name='z_d')

    
    GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost')

    OF_MP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_MP')             # obj func var
    
    
    # Contingency groups and recourse-cost variables
    all_sub_cont_list=[]
    for b in Bus:
        all_sub_cont_list+=[b,str(b+'-1'),str(b+'-2')]
            
    all_cont=all_sub_cont_list+line_cont_list #used for reporting LoadShedding_df
    

    if substation_OSP==True and all_cont_OSP==False:
        Phi_b=model.addVars(Bus,lb=0,vtype=GRB.CONTINUOUS,name='Phi_b')
        TotalPhi=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalPhi')
        
    elif substation_OSP==False and all_cont_OSP==True:
        Phi_c=model.addVars(all_cont,lb=0,vtype=GRB.CONTINUOUS,name='Phi_b')
        TotalPhi=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalPhi')  
    else:
        print('\n\nwrong OSP input!')


    # Initial topology
    if Zfixdict==None:
        for b in Bus:
            z_bus[b].Start=1 
        for l,i,j in Lines:
            z_li[l,i,j].Start=0
        for g in G:
            z_g[g].Start=0             
        for d in DemandSet:
            z_d[d].Start=0
        
        
    # Fix topology 
    if Zfixdict!=None:
        model.addConstrs((  z_bus[b] == Zfixdict['bus'][b]   for b in Bus), name='Fix_bus')
        model.addConstrs((  z_g[g] == Zfixdict['g'][g]   for g in G), name='Fix_g')
        model.addConstrs((  z_d[d] == Zfixdict['d'][d]   for d in DemandSet), name='Fix_d')
        model.addConstrs((  z_li[l,i,j] == Zfixdict['l_i'][l,i,j]   for l,i,j in Lines), name='Fix_li')
            
            
    # Fix dispatch
    if PQgFix!=None:
        model.addConstrs((  Pgi[g,i] == PQgFix['Pg'][g,i]   for g in G for i in busbar), name='Fix_Pg')
        model.addConstrs((  Qgi[g,i] == PQgFix['Qg'][g,i]   for g in G for i in busbar), name='Fix_Qg')
        
    if PgFix!=None:
        model.addConstrs( (   Pg[g]  ==  PgFix[g] 
                             for g in G ), name='eq_Pgfix')
        
        model.addConstrs( (   Pgi[g,'busbar1']  == (1-z_g[g])*PgFix[g] 
                             for g in G ), name='eq_Pg1')
        model.addConstrs( (  Pgi[g,'busbar2'] == z_g[g]*PgFix[g] 
                             for g in G  ), name='eq_Pg2')
        
        
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
    
    eq_MaxSw_bus=model.addConstr((  quicksum( (1-z_bus[b]) for b in Bus) <= Max_Sw_bus) , name='eq_MaxSw_bus') 
    

    for b in Bus:
        lmin,b=L2B.select('*',b)[0]
        lmin,i,j=Lines.select(lmin,b,'*')[0]
        eq_symmetry1=model.addConstr((  z_li[lmin,i,j] ==0  ), name='eq_symmetry')

    eq_zbus_ref=model.addConstr( (  z_bus[Bus[0]]==1    ), name='eq_zbus_ref')

        
    eq_Pg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar1']   
                             for g in G), name='eq_Pg1min')
    eq_Pg1max=model.addConstrs( (    Pgi[g,'busbar1'] <=(1-z_g[g])*Gen_data.loc[g]['Pmax']  
                             for g in G  ), name='eq_Pg1max')
    
    eq_Pg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar2']   
                             for g in G   ), name='eq_Pg2min')
    eq_Pg2max=model.addConstrs( (   Pgi[g,'busbar2']  <= z_g[g]*Gen_data.loc[g]['Pmax']   
                             for g in G   ), name='eq_Pg2max')
    

    model.addConstrs( (   Pg[g] == Pgi[g,'busbar1']+Pgi[g,'busbar2']   for g in G   ), name='eq_Pg')
    model.addConstrs( (   Qg[g] == Qgi[g,'busbar1']+Qgi[g,'busbar2']   for g in G   ), name='eq_Qg')

    
    # Reactive generation
    eq_Qg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar1']   
                             for g in G   ), name='eq_Qg1min')
    eq_Qg1max=model.addConstrs( (    Qgi[g,'busbar1'] <=(1-z_g[g])*Gen_data.loc[g]['Qmax']  
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

    
    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar1'] 
                             for l,i,j in Lines ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Pflow_li[l,i,j,'busbar1'] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar2'] 
                            for l,i,j in Lines) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Pflow_li[l,i,j,'busbar2'] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines) ,name='eq_flow2max')
    
    eq_flow=model.addConstrs((   Pflow[l,i,j] == Pflow_li[l,i,j,'busbar1']+Pflow_li[l,i,j,'busbar2']  
                                for l,i,j in Lines ), name='eq_Pflow')
    
    
    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar1'] 
                             for l,i,j in Lines) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Qflow_li[l,i,j,'busbar1'] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar2'] 
                            for l,i,j in Lines) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Qflow_li[l,i,j,'busbar2'] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines) ,name='eq_flow2max')
    
    eq_flow=model.addConstrs((   Qflow[l,i,j] == Qflow_li[l,i,j,'busbar1']+Qflow_li[l,i,j,'busbar2']  
                                for l,i,j in Lines), name='eq_Qflow')
    
    # Tight linear AC limits
    model.addConstrs( (  Pflow[l,i,j] + np.tan(np.pi/6)*Qflow[l,i,j] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines) , name='eqPij1')
    
    model.addConstrs( (  Pflow[l,i,j] - np.tan(np.pi/6)*Qflow[l,i,j] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines) , name='eqPij2')
    
    model.addConstrs( (  Pflow[l,i,j] + np.tan(np.pi/6)*Qflow[l,i,j] >= (beta-1)*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines) , name='eqPij3')
     
    model.addConstrs( (  Pflow[l,i,j] - np.tan(np.pi/6)*Qflow[l,i,j] >= (beta-1)*branch.loc[(l,i,j)]['limit'] 
                            for l,i,j in Lines) , name='eqPij4')

    
    # Linear AC branch-flow equations

    model.addConstrs( (  Pflow[l,i,j] == 
    0.5*branch.loc[(l,i,j)]['g_ij']*( V2_li[l,i,j] - V2_li[l,j,i] )
    -branch.loc[(l,i,j)]['b_ij']*(delta_li[l,i,j]-delta_li[l,j,i]) #+ Ploss[l,i,j]
                            for l,i,j in Lines ) , name='eqPij')
    
    model.addConstrs( (  Qflow[l,i,j] == 
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j] - V2_li[l,j,i])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j]-delta_li[l,j,i]) 
                            for l,i,j in Lines) , name='eqQij')
    
    # Branch losses


    # Voltage magnitude
    model.addConstrs( ( V2_bi[b,i] <= vmax**2  for b in Bus for i in busbar ) , name='eqvmaxb')
    model.addConstrs( ( V2_bi[b,i] >= vmin**2  for b in Bus for i in busbar ) , name='eqvminb')
    model.addConstrs( ( V2_li[l,i,j] <= vmax**2 for l,i,j in Lines ) , name='eqvmaxl')
    model.addConstrs( ( V2_li[l,i,j] >= vmin**2  for l,i,j in Lines ) , name='eqvminl')


    # Voltage angle
    
    eq_delta_bus1=model.addConstrs((   -BigM_b*(1-z_bus[b]) <= delta_bi[b,'busbar1']-delta_bi[b,'busbar2'] 
                               for b in Bus ), name='eq_delta_bus1')
    eq_delta_bus2=model.addConstrs((   delta_bi[b,'busbar1']-delta_bi[b,'busbar2'] <= BigM_b*(1-z_bus[b])
                               for b in Bus ), name='eq_delta_bus2')
    model.addConstrs(( -z_li[l,i,j]*Maxdelta<= delta_li[l,i,j] - delta_bi[i,'busbar1']   
                                for l,i,j in Lines), name='eq_delta_lb1Frmin')
    model.addConstrs((  delta_li[l,i,j] - delta_bi[i,'busbar1'] <= z_li[l,i,j]*Maxdelta   
                                for l,i,j in Lines), name='eq_delta_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*Maxdelta<= delta_li[l,i,j] - delta_bi[i,'busbar2']   
                                for l,i,j in Lines), name='eq_delta_lb2Frmin')
    model.addConstrs((  delta_li[l,i,j] - delta_bi[i,'busbar2'] <= (1-z_li[l,i,j])*Maxdelta   
                                for l,i,j in Lines), name='eq_delta_lb2Frmax')
    
    # Squared voltage magnitude
    eq_V2_bus1=model.addConstrs((   -MaxV2*(1-z_bus[b]) <= V2_bi[b,'busbar1']-V2_bi[b,'busbar2'] 
                               for b in Bus ), name='eq_V2_bus1')
    eq_V2_bus2=model.addConstrs((   V2_bi[b,'busbar1']-V2_bi[b,'busbar2'] <= MaxV2*(1-z_bus[b])
                               for b in Bus ), name='eq_V2_bus2')
    model.addConstrs(( -z_li[l,i,j]*MaxV2<= V2_li[l,i,j] - V2_bi[i,'busbar1']   
                                for l,i,j in Lines), name='eq_V2_lb1Frmin')
    model.addConstrs((  V2_li[l,i,j] - V2_bi[i,'busbar1'] <= z_li[l,i,j]*MaxV2   
                                for l,i,j in Lines), name='eq_V2_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*MaxV2<= V2_li[l,i,j] - V2_bi[i,'busbar2']   
                                for l,i,j in Lines), name='eq_V2_lb2Frmin')
    model.addConstrs((  V2_li[l,i,j] - V2_bi[i,'busbar2'] <= (1-z_li[l,i,j])*MaxV2   
                                for l,i,j in Lines), name='eq_V2_lb2Frmax')

    
    eq_delta_ref=model.addConstr((delta_bi[Bus[0],'busbar1']==0      ), name='ref_bus_angle' )
    
    
    model.addConstrs((  Pflow_bus[b]<=BigM_busbar*(z_bus[b]) for b in Bus), name='eq_busbarPflowmax')
    model.addConstrs((  -Pflow_bus[b]<=BigM_busbar*(z_bus[b]) for b in Bus), name='eq_busbarPflowmin')
    
    model.addConstrs((  Qflow_bus[b]<=BigM_busbar*(z_bus[b]) for b in Bus), name='eq_busbarQflowmax')
    model.addConstrs((  -Qflow_bus[b]<=BigM_busbar*(z_bus[b]) for b in Bus ), name='eq_busbarQflowmin')
    
    
    eq_balance_busbar1=model.addConstrs((                                        
                quicksum( Pgi[g,m]   for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m]   for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m] for l,i,j in Lines.select('*',b,'*'))
                +Pflow_bus[b]
                            for b in Bus for m in busbar if m=='busbar1' ), name='eq_balance1')
    eq_balance_busbar2=model.addConstrs((                                        
                quicksum( Pgi[g,m] for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m] for l,i,j in Lines.select('*',b,'*'))
                -Pflow_bus[b]
                            for b in Bus for m in busbar if m=='busbar2' ), name='eq_balance2')
    
    
    eq_Qbalance_busbar1=model.addConstrs((                                        
                quicksum( Qgi[g,m]   for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m]   for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m] for l,i,j in Lines.select('*',b,'*'))
                +Qflow_bus[b]
                            for b in Bus for m in busbar  if m=='busbar1' ), name='eq_balance1')
    eq_Qbalance_busbar2=model.addConstrs((                                        
                quicksum( Qgi[g,m]for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m] for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m] for l,i,j in Lines.select('*',b,'*'))
                -Qflow_bus[b]
                            for b in Bus for m in busbar if m=='busbar2'  ), name='eq_balance2')

    
    # Master objective
    model.addConstr(GenCost == quicksum(Gen_data.loc[g]['b']*Pg[g] for g in G) ,name='Eq_GC')
    
    
    if substation_OSP==True and all_cont_OSP==False:
        model.addConstr( TotalPhi == quicksum(Phi_b[b] for b in Bus)   ,name='Eq_Phi')
        OF_MP = GenCost + TotalPhi + quicksum(0.001*Gen_data.loc[g]['b']*Qg[g]*Qg[g] for g in G) 
    elif substation_OSP==False and all_cont_OSP==True:
        model.addConstr( TotalPhi == quicksum(Phi_c[c] for c in all_cont)   ,name='Eq_Phi')
        OF_MP = GenCost + TotalPhi + quicksum(0.001*Gen_data.loc[g]['b']*Qg[g]*Qg[g] for g in G)
        
        
    model.setObjective(OF_MP,GRB.MINIMIZE)

    model.update()
    NumBinVarsMP = model.NumBinVars
    NumConVarsMP = model.NumVars - model.NumBinVars
    NumConstrsMP = model.NumConstrs

    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time  #execution time
    
    tot_time = 0 
    tot_time += ex_time

    time_iteration = {}
    time_iter_detail = {'MP':{},'SP':{}}
    MP_time=0; SP_time=0
    MP_time+=ex_time
    
    
    # Build the Topology dictionary
    def UpdateTopology():
        PgMP={}; QgMP={} ;z_busMP={}; z_gMP={}; z_dMP={}; z_liMP={}
        
        for g in G:
                PgMP[g]=Pg[g].x
                QgMP[g]=Qg[g].x
        for b in Bus:
            z_busMP[b]=z_bus[b].x
        for g in G:
            z_gMP[g]=z_g[g].x
        for d in DemandSet:
            z_dMP[d]=z_d[d].x
        for l,i,j in Lines:
            z_liMP[(l,i,j)]=z_li[l,i,j].x
        
        return {
            'Pg':PgMP,
            'Qg':QgMP,
            'bus':z_busMP,
              'g':z_gMP,
              'd':z_dMP,
              'l_i':z_liMP}
    
    TopologyMP=UpdateTopology()
    
    
    # Benders iteration state
    NumberofCuts=0
    
    Mu={}
    UB_k={}; LB_k={}
    UB=1e6; LB=0 

    
    k1=0         #FSP inner loop
    k2=0         #OSP outer loop
    GenCost_k1={}
    FSP_Obj_k1={}
    Topology_k={}
    UB_min=UB
    
    
    ShedCost_df = pd.DataFrame(columns=['ShedCost(c)'],index=all_cont, dtype=float)
    PdShed_dic = pd.DataFrame(columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,busbar,all_cont]), dtype=float)


    # Create reusable FSP and OSP models.
    print('creating FSP, OSP models.')
    FSP_models={}
    OSP_models={}
    NumBinVars = {'SP':{}}; NumConVars = {'SP':{}}; NumConstrs = {'SP':{}}

    NumBinVars['MP'] = NumBinVarsMP
    NumConVars['MP'] = NumConVarsMP
    NumConstrs['Mp'] = NumConstrsMP

    for c in line_cont_list:
        res = create_FSP_line_v62(data,cont_list=[c])
        FSP_models[c] = res['model']
        NumBinVars['SP'][c] = res['NumBinVars'] 
        NumConVars['SP'][c] = res['NumConVars'] 
        NumConstrs['SP'][c] = res['NumConstrs']

    if substation_OSP==True and all_cont_OSP==False:
        for b in Bus:
            res = create_OSP_substation_v62(data,cont_list=[b,str(b+'-1'),str(b+'-2')]
                                                    ,Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,Max_Sw_bus=Max_Sw_bus)
            OSP_models[b] = res['model']
            NumBinVars['SP'][b] = res['NumBinVars'] 
            NumConVars['SP'][b] = res['NumConVars'] 
            NumConstrs['SP'][b] = res['NumConstrs']

            
    elif substation_OSP==False and all_cont_OSP==True:
        for c in all_cont:
            res =create_OSP_substation_v62(data,cont_list=[c]
                                                    ,Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,Max_Sw_bus=Max_Sw_bus)
            
            OSP_models[c] = res['model']
            NumBinVars['SP'][c] = res['NumBinVars'] 
            NumConVars['SP'][c] = res['NumConVars'] 
            NumConstrs['SP'][c] = res['NumConstrs']


    current_time = time.localtime()
    formatted_time = time.strftime("%Y-%m-%d %H:%M:%S", current_time)
    print("Current time and date:", formatted_time)
    
    # OSP outer loop
    for k2 in tqdm(range(Max_OSP_iter)):
        
        k1=0
        # FSP inner loop
        while k1 <= Max_FSP_iter:

            GenCost_k1[k1]=GenCost.x


            if FSP_single==True: #default

                ex_time=0

                FSP_Objc=0
                FSP_Obj=0
                

                for c in (line_cont_list):
                    result = solve_FSP_line_v62(data=data,TopologyMP=TopologyMP,model0=FSP_models[c],print_result=False)
                    FSP_Objc=result['OF_FSP']
                    FSP_Obj+=FSP_Objc
                    Mu=result['Mu']
                    tot_time += result['time']
                    SP_time+=result['time']

                    if FSP_Objc >= FSP_criteria:
                        # Add either a classical BD-C cut or heuristic BD-H cut.
                        if heuristic_cut==False:
                            model.addConstr(  ( FSP_Objc   
                                        +quicksum( Mu['Pg'][g]*(Pg[g]-TopologyMP['Pg'][g])  for g in G)
                                        +quicksum( Mu['Qg'][g]*(Qg[g]-TopologyMP['Qg'][g])  for g in G)
                                        +quicksum( Mu['bus'][b]*(z_bus[b]-TopologyMP['bus'][b])  for b in Bus)
                                        +quicksum( Mu['g'][g]*(z_g[g]-TopologyMP['g'][g])  for g in G)
                                        +quicksum( Mu['d'][d]*(z_d[d]-TopologyMP['d'][d])  for d in DemandSet)
                                        +quicksum( Mu['l_i'][(l,i,j)]*(z_li[l,i,j]-TopologyMP['l_i'][l,i,j])  for l,i,j in Lines)
                                                <= 0   )
                                ,name=str('FSP-'+str(NumberofCuts)))
                            NumberofCuts+=1
                        elif heuristic_cut==True:
                            model.addConstr(  ( FSP_Objc   
                                        +quicksum( Mu['Pg'][g]*(Pg[g]-TopologyMP['Pg'][g])  for g in G)
                                        +quicksum( Mu['Qg'][g]*(Qg[g]-TopologyMP['Qg'][g])  for g in G)
                                         +quicksum( Mu['bus'][b]*(z_bus[b]-TopologyMP['bus'][b])  for b in Bus)
                                        +quicksum( Mu['g'][g]*(z_g[g]-TopologyMP['g'][g])  for g in G)
                                        +quicksum( Mu['d'][d]*(z_d[d]-TopologyMP['d'][d])  for d in DemandSet)
                                        +quicksum( Mu['l_i'][(l,i,j)]*(z_li[l,i,j]-TopologyMP['l_i'][l,i,j])  for l,i,j in Lines)
                                                <= 0   )
                                ,name=str('FSP-'+str(NumberofCuts)))
                            NumberofCuts+=1

                FSP_Obj_k1[k1] = FSP_Obj

                if FSP_Obj > FSP_criteria:
                    k1+=1
                    model.update()
                    start_time = time.time()
                    model.optimize()
                    end_time = time.time()
                    ex_time=end_time-start_time  #execution time
                    
                    tot_time += ex_time
                    MP_time+=ex_time

                    TopologyMP=UpdateTopology()
                else:
                    break

        if FSP_Obj>FSP_criteria:
            print('FSP obj is: %.3f'%FSP_Obj)


        Topology_k[k2] = TopologyMP
        
        # Substation-grouped OSP evaluation
        if substation_OSP==True and all_cont_OSP==False: # substation cont_list
        
            OF_OSP=0 ; OF_OSP_b={}
            Mu={}

            for b in (Bus): #OSP========== 
                sub_cont_list = []
                sub_cont_list+=[b,str(b+'-1'),str(b+'-2')]
                
                
                result = solve_OSP_substation_v62(data=data,model0=OSP_models[b],TopologyMP=TopologyMP, cont_list=sub_cont_list)
                OF_OSP_b[b]=result['OF_OSP']
                OF_OSP+=OF_OSP_b[b]
                Mu[b]=result['Mu']
                tot_time += result['time']
                SP_time+=result['time']
                

                for cont in all_cont:
                    if cont in result['ShedCost_df'].index:
                        ShedCost_df.loc[cont] = result['ShedCost_df'].loc[cont]
                        for d in DemandSet:
                            for i in busbar:
                                PdShed_dic.loc[(d,i,cont)] = result['PdShed_dic'].loc[(d,i,cont)]
                

            obj = model.getObjective().getValue() #OF_MP.x
            
            LB=obj ; LB_k[k2]=LB    
            UB=GenCost.x+OF_OSP ; UB_k[k2]=UB
            time_iteration[k2] = tot_time

            time_iter_detail['MP'][k2]=MP_time; MP_time=0
            time_iter_detail['SP'][k2]=SP_time; SP_time=0
            
            
            k_min = min(UB_k, key=UB_k.get)
            UB_min = UB_k[k_min]
            

            if (UB_min-LB) >= OSP_criteria:
                
                # Add either global classical cuts or localized heuristic cuts.
                if heuristic_cut==False:
                    for sub in Bus: #multiple cut
                        model.addConstr(  ( OF_OSP_b[sub] 
                                +quicksum( Mu[sub]['Pg'][g]*(Pg[g]-TopologyMP['Pg'][g])  for g in G)
                                +quicksum( Mu[sub]['Qg'][g]*(Qg[g]-TopologyMP['Qg'][g])  for g in G)
                                +quicksum( Mu[sub]['bus'][b]*(z_bus[b]-TopologyMP['bus'][b])  for b in Bus )
                                +quicksum( Mu[sub]['g'][g]*(z_g[g]-TopologyMP['g'][g])  for g in G)
                                +quicksum( Mu[sub]['d'][d]*(z_d[d]-TopologyMP['d'][d])  for d in DemandSet )
                                +quicksum( Mu[sub]['l_i'][(l,i,j)]*(z_li[l,i,j]-TopologyMP['l_i'][l,i,j])  for l,i,j in Lines)
                                            <= Phi_b[sub]  )
                            ,name=str('OSP-'+sub+'-'+str(NumberofCuts)))
                        model.addConstr(  Phi_b[sub].x <= Phi_b[sub] ,name=str('OSP-Phi-'+sub+'-'+str(NumberofCuts)))
                        NumberofCuts+=2

                
                if heuristic_cut==True:
                    if hops_cut ==0:
                        for sub in Bus:  #multi heuristic cuts
                            model.addConstr(  ( OF_OSP_b[sub] 
                                +quicksum( Mu[sub]['Pg'][g]*(Pg[g]-TopologyMP['Pg'][g])  
                                          for g in G )
                                +quicksum( Mu[sub]['Qg'][g]*(Qg[g]-TopologyMP['Qg'][g])  
                                          for g in G )
                                +quicksum( Mu[sub]['bus'][b]*(z_bus[b]-TopologyMP['bus'][b])  for b in Bus if b==sub )
                                +quicksum( Mu[sub]['g'][g]*(z_g[g]-TopologyMP['g'][g]) for g,b in G2B.select('*',sub))
                                +quicksum( Mu[sub]['d'][d]*(z_d[d]-TopologyMP['d'][d])  for d,b in D2B.select('*',sub) )
                                +quicksum( Mu[sub]['l_i'][(l,i,j)]*(z_li[l,i,j]-TopologyMP['l_i'][l,i,j])  
                                          for l,i,j in Lines.select('*',sub,'*'))
                                            <= Phi_b[sub]  )
                            ,name=str('OSP-'+sub+'-'+str(NumberofCuts)))

                            model.addConstr(  Phi_b[sub].x <= Phi_b[sub] ,name=str('OSP-Phi-'+sub+'-'+str(NumberofCuts)))
                            NumberofCuts+=2
                    
                    elif hops_cut >= 1:
                        for sub in Bus:  #multi hops cuts
                            model.addConstr(  ( OF_OSP_b[sub] 
                                +quicksum( Mu[sub]['Pg'][g]*(Pg[g]-TopologyMP['Pg'][g])  
                                          for g in G )
                                +quicksum( Mu[sub]['Qg'][g]*(Qg[g]-TopologyMP['Qg'][g])  
                                          for g in G )
                       +quicksum( Mu[sub]['bus'][b]*(z_bus[b]-TopologyMP['bus'][b])  
                             for b in Bus if b in Find_nodes_hops(data,sub,hops_cut) )
                    +quicksum( Mu[sub]['g'][g]*(z_g[g]-TopologyMP['g'][g]) 
                              for b in Bus for g,i in G2B.select('*',b) if b in Find_nodes_hops(data,sub,hops_cut))
                    +quicksum( Mu[sub]['d'][d]*(z_d[d]-TopologyMP['d'][d])  
                               for b in Bus for d,b in D2B.select('*',b) if b in Find_nodes_hops(data,sub,hops_cut) )
                    +quicksum( Mu[sub]['l_i'][(l,i,j)]*(z_li[l,i,j]-TopologyMP['l_i'][l,i,j])  
                                  for l,i,j in Lines if l in Find_lines_hops(data,sub,hops_cut))
                                            <= Phi_b[sub]  )
                            ,name=str('OSP-'+sub+'-'+str(NumberofCuts)))

                            model.addConstr(  Phi_b[sub].x <= Phi_b[sub] ,name=str('OSP-Phi-'+sub+'-'+str(NumberofCuts)))
                            NumberofCuts+=2
                    
                    
                k2+=1
                model.update()
                start_time = time.time()
                model.optimize()
                end_time = time.time()
                ex_time=end_time-start_time  #execution time
                
                tot_time += ex_time
                MP_time += ex_time

                TopologyOld=TopologyMP.copy()
                TopologyMP=UpdateTopology()
                
                
            elif (UB_min-LB)<OSP_criteria:
                print('Gap is %.2f, OSP loop break!'%(UB_min-LB))
                time_iteration[k2] = tot_time
                break

                
        # Individual-contingency OSP evaluation
        elif substation_OSP==False and all_cont_OSP==True: # all cont_list
        
            OF_OSP=0 ; OF_OSP_c={}
            Mu={}


            for c in (all_cont): #OSP========== 
                
                    
                result = solve_OSP_substation_v62(data=data,cont_list=[c],TopologyMP=TopologyMP,model0=OSP_models[c])
                OF_OSP_c[c]=result['OF_OSP']
                OF_OSP+=OF_OSP_c[c]
                Mu[c]=result['Mu']
                tot_time += result['time']
                SP_time+=result['time']
                
                for cont in all_cont:
                    if cont in result['ShedCost_df'].index:
                        ShedCost_df.loc[cont] = result['ShedCost_df'].loc[cont]
                        for d in DemandSet:
                            for i in busbar:
                                PdShed_dic.loc[(d,i,cont)] = result['PdShed_dic'].loc[(d,i,cont)]
                

            obj = model.getObjective().getValue() #OF_MP.x
            LB=obj ; LB_k[k2]=LB    
            UB=GenCost.x+OF_OSP ; UB_k[k2]=UB
            time_iteration[k2] = tot_time
            time_iter_detail['MP'][k2]=MP_time; MP_time=0
            time_iter_detail['SP'][k2]=SP_time; SP_time=0
            
            k_min = min(UB_k, key=UB_k.get)
            UB_min = UB_k[k_min]
            

            if (UB_min-LB) >= OSP_criteria:
                
                for c in all_cont: #multiple cut
                        model.addConstr(  ( OF_OSP_c[c] 
                                +quicksum( Mu[c]['Pg'][g]*(Pg[g]-TopologyMP['Pg'][g])  for g in G)
                                +quicksum( Mu[c]['Qg'][g]*(Qg[g]-TopologyMP['Qg'][g])  for g in G )
                                +quicksum( Mu[c]['bus'][b]*(z_bus[b]-TopologyMP['bus'][b])  for b in Bus )
                                +quicksum( Mu[c]['g'][g]*(z_g[g]-TopologyMP['g'][g])  for g in G)
                                +quicksum( Mu[c]['d'][d]*(z_d[d]-TopologyMP['d'][d])  for d in DemandSet )
                                +quicksum( Mu[c]['l_i'][(l,i,j)]*(z_li[l,i,j]-TopologyMP['l_i'][l,i,j])  for l,i,j in Lines)
                                            <= Phi_c[c]  )
                            ,name=str('OSP-'+str(NumberofCuts)))
                        model.addConstr(  Phi_c[c].x <= Phi_c[c] ,name=str('OSP-Phi-'+str(NumberofCuts)))
                        NumberofCuts+=2

                
                k2+=1
                model.update()
                start_time = time.time()
                model.optimize()
                end_time = time.time()
                ex_time=end_time-start_time  #execution time
                
                tot_time += ex_time
                MP_time += ex_time

                TopologyOld=TopologyMP.copy()
                TopologyMP=UpdateTopology()
                
                
            elif (UB_min-LB)<OSP_criteria:
                print('Gap is %.2f, OSP loop break!'%(UB_min-LB))
                time_iteration[k2] = tot_time
                break
      
                
        else:
            print('\n\nwrong OSP input!')
            break #while OSP loop
                
            
    print('Best solution found in iteration: ',k_min, UB_min)


    # Fix pre-contingency dispatch
    Pgi_dict={}; Qgi_dict={}
    for g in G:
        for i in busbar:
            Pgi_dict[(g,i)]=Pgi[g,i].x
            Qgi_dict[(g,i)]=Qgi[g,i].x
            
    Pg_dict = {}; Qg_dict = {}
    for g in G:
        Pg_dict[g] = Pgi[g,'busbar1'].x*(1-z_g[g].x) + Pgi[g,'busbar2'].x*(z_g[g].x)
        Qg_dict[g] = Qgi[g,'busbar1'].x*(1-z_g[g].x) + Qgi[g,'busbar2'].x*(z_g[g].x)
    
    
    # Results
    return {
        'TopologyDict':TopologyMP,
           'UB_k':UB_k,
           'LB_k':LB_k,
        'Cost_tot':UB,
        'GenCost':GenCost.x,
        'TotalShedCost':OF_OSP,
                'Pg':Pg_dict,
        'Qg':Qg_dict,
        'Pgi':Pgi_dict,
        'Qgi':Qgi_dict,
        'ShedCost_df':ShedCost_df,
         'PdShed_dic': PdShed_dic,
        'time':tot_time,
        'time_iteration':time_iteration,
        'time_iter_detail':time_iter_detail,
        'NumBinVars':NumBinVars,
        'NumConVars':NumConVars,
        'NumConstrs':NumConstrs
    }

    


