"""Iterative topology-fixing heuristic baseline from Section IV-A.

1-Opt-H: an iterative one-step improvement heuristic
that starts from the initial topology T_0 and,
at each iteration, flips a single binary variable that reduces the objective,
until no improving move exists.
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


def create_iterative_heuristic_model(
    data,
    cont_list=None,
    Max_Sw_bus=0,
    Up_redispatch=0,
    Dn_redispatch=1.0,
    Zfixdict=None,
    Zinitial=None,
    PQgFix=None,
    PgFix=None,
    FixedCost=False,
    Alpha=0,
    Pg_market=None,
    Probabilistic=False,
    SolverTime=600,
    Threads=None,
    line_shedding=True,
    Vol0_li=None,
    delta0_li=None,
):
    """Build the reusable MIP for the iterative topology heuristic.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        cont_list: Contingency identifiers represented in the model.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.
        Up_redispatch: Fraction of generator capacity available for upward
            active-power redispatch.
        Dn_redispatch: Fraction of generator capacity available for downward
            active-power redispatch.
        Zfixdict: Optional dictionary of topology variables to fix.
        Zinitial: Optional initial topology used as a Gurobi warm start.
        PQgFix: Optional busbar-level active and reactive generation to fix.
        PgFix: Optional total active-generation values to fix.
        FixedCost: Enforce the market-generation cost limit when ``True``.
        Alpha: Relative margin applied to the market-generation cost limit.
        Pg_market: Reference market dispatch used by the fixed-cost constraint.
        Probabilistic: Weight contingency shedding costs by their probabilities
            when ``True``.
        SolverTime: Gurobi time limit in seconds.
        Threads: Optional number of Gurobi solver threads.
        line_shedding: Line-shedding configuration flag retained from the case
            study interface.
        Vol0_li: Optional branch-voltage linearization points; unit values are
            used when omitted.
        delta0_li: Optional branch-angle linearization points; zero values are
            used when omitted.

    Returns:
        The configured, unsolved Gurobi model.
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

    
    # Model configuration
    model=gp.Model('SBF_Preventive SC line contingency')
    model.Params.OutputFlag=1
    
    cont_list=cont_list.copy()
    # The contingency list includes the normal operating state (0).
    
    
    model.Params.MIPGap=Max_MIPGap
    model.Params.timelimit=SolverTime
    if Threads is not None:
        model.Params.Threads = Threads

    
    Max_Sw_bus=Max_Sw_bus  
  

    if Vol0_li is None:
        print('Vol0 is None')
        Vol0_li={}
        for l,i,j in Lines:
            for c in cont_list:
                Vol0_li[l,i,j,c]=1.0

    if delta0_li is None:
        delta0_li={}
        for l,i,j in Lines:
            for c in cont_list:
                delta0_li[l,i,j,c]=0.0
    
    
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
    V2_bi=model.addVars(Bus,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2_bi') 
    V2_li=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2_li') 
    
   
    z_bus=model.addVars(Bus,vtype=GRB.BINARY,name='z_bus')    
    z_li=model.addVars(Lines,vtype=GRB.BINARY,name='z_li')  
    z_g=model.addVars(G,vtype=GRB.BINARY,name='z_g') 
    z_d=model.addVars(DemandSet,vtype=GRB.BINARY,name='z_d') 

    OF=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF')             
    GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost')
    RDCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='RDCost')
    ShedCost=model.addVars(cont_list,lb=0,vtype=GRB.CONTINUOUS,name='ShedCost')
    Shed_c=model.addVars(cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Shed_c')

    TotalShedCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalShedCost')


    if Zfixdict==None:
        for b in Bus:
            z_bus[b].Start=1 
        for l,i,j in Lines:
            z_li[l,i,j].Start=0
        for g in G:
            z_g[g].Start=0             
        for d in DemandSet:
            z_d[d].Start=0
        
        
    # Fixed topology
    if Zfixdict!=None:
        model.addConstrs((  z_bus[b] == Zfixdict['bus'][b]   for b in Bus), name='Fix_bus')
        model.addConstrs((  z_g[g] == Zfixdict['g'][g]   for g in G), name='Fix_g')
        model.addConstrs((  z_d[d] == Zfixdict['d'][d]   for d in DemandSet), name='Fix_d')
        model.addConstrs((  z_li[l,i,j] == Zfixdict['l_i'][l,i,j]   for l,i,j in Lines), name='Fix_li')
            
            
    # Fixed dispatch
    if PQgFix!=None:
        model.addConstrs((  Pgi[g,i] == PQgFix['Pg'][g,i]   for g in G for i in busbar), name='Fix_Pg')
        model.addConstrs((  Qgi[g,i] == PQgFix['Qg'][g,i]   for g in G for i in busbar), name='Fix_Qg')
        
    if PgFix!=None:
        model.addConstrs( (   Pgi[g,'busbar1']  == (1-z_g[g])*PgFix[g] 
                             for g in G ), name='eq_Pg1')
        model.addConstrs( (  Pgi[g,'busbar2'] == z_g[g]*PgFix[g] 
                             for g in G  ), name='eq_Pg2')
        
        
    if Zinitial!=None:    
        for b in Bus:
            z_bus[b].Start=Zinitial['bus'][b] 
        for l,i,j in Lines:
            z_li[l,i,j].Start=Zinitial['l_i'][l,i,j]
        for g in G:
            z_g[g].Start=Zinitial['g'][g]             
        for d in DemandSet:
            z_d[d].Start=Zinitial['d'][d] 
        
        
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
    

    if Zfixdict==None:
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
    
    model.addConstrs( (   Pg[g] == Pgi[g,'busbar1']+Pgi[g,'busbar2']   for g in G   ), name='eq_Pg')
    model.addConstrs( (   Qg[g] == Qgi[g,'busbar1']+Qgi[g,'busbar2']   for g in G   ), name='eq_Qg')
    
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
    model.addConstrs(( dPgi_dn[g,i,c] == 0  
                            for g in G for i in busbar for c in cont_list if c==0 ),name='Eq_dPgDn_res')
    
    Eq_dQgUp_res=model.addConstrs(( dQgi_up[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #Up_redispatch  
                          for g in G for i in busbar for c in cont_list ),name='Eq_dQgUp_res')
    Eq_dQgDn_res=model.addConstrs(( dQgi_dn[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #Dn_redispatch  
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
    model.addConstrs( (  Ploss[l,i,j,c] == Ploss_vol[l,i,j,c] + Ploss_delta[l,i,j,c]
                            for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqPloss_total')
    model.addConstrs( (  Qloss[l,i,j,c] == Qloss_vol[l,i,j,c] + Qloss_delta[l,i,j,c]
                            for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqQloss_total')

    model.addConstrs( (  0 <= Ploss[l,i,j,c] + epsilon[l,i,j,c]
                            for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqPloss_pos')

    model.addConstrs( (  Ploss_vol[l,i,j,c] ==
                branch.loc[(l,i,j)]['g_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                -0.5*branch.loc[(l,i,j)]['g_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                            for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqPloss_vol')
    model.addConstrs( (  Qloss_vol[l,i,j,c] ==
                -branch.loc[(l,i,j)]['b_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])
                +0.5*branch.loc[(l,i,j)]['b_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2
                            for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqQloss_vol')

    model.addConstrs( (  Ploss_delta[l,i,j,c] ==
                branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                -0.5*branch.loc[(l,i,j)]['g_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                            for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqPloss_delta')
    model.addConstrs( (  Qloss_delta[l,i,j,c] ==
                -branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])*(delta_li[l,i,j,c]-delta_li[l,j,i,c])
                +0.5*branch.loc[(l,i,j)]['b_ij']*(delta0_li[l,i,j,c]-delta0_li[l,j,i,c])**2
                            for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqQloss_delta')


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
    

    # Active- and reactive-power balance
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
    
    model.addConstr(GenCost == quicksum(Gen_data.loc[g]['b']*Pg[g] for g in G) ,name='Eq_GC')
    
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


    if FixedCost==True:

        model.addConstr(OF ==  RDCost + TotalShedCost  ,name='Eq_OF')
    
        model.addConstr(GenCost <= (1+Alpha)*quicksum(Gen_data.loc[g]['b']*Pg_market[g] for g in G) ,name='Eq_GC2')

    else:
        model.addConstr(OF == GenCost + RDCost + TotalShedCost  ,name='Eq_OF')


    model.setObjective(OF,GRB.MINIMIZE)

    model.update()


    # Configured model
    return model
    
    
def solve_iterative_heuristic_iteration(
    data,
    TopologyFix,
    model0,
    cont_list,
    K_Opt=1,
    print_result=False,
):
    """Solve one topology-refinement iteration.

    Previously accepted topology decisions are fixed, while at most ``K_Opt``
    additional decisions may change during this solve.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        TopologyFix: Topology decisions fixed by earlier iterations.
        model0: Reusable MIP created by ``create_iterative_heuristic_model``.
        cont_list: Contingency identifiers represented by ``model0``.
        K_Opt: Maximum number of additional topology decisions to optimize.
        print_result: Print detailed solution information when ``True``.

    Returns:
        A dictionary containing solve time, topology decisions, updated fixed
        decisions, dispatch, optimality gap, costs, load shedding, and updated
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

    
    # Copy the reusable model.
    model=model0.copy()
    model.Params.OutputFlag=1


    # Fix previously accepted topology decisions.
    counter = 0
    for b in Bus:
        if b in TopologyFix['bus'].keys():
            model.addConstr(  (model.getVarByName('z_bus['+str(b)+']') == TopologyFix['bus'][b] )  ,name='eqBendersZ_bus')
            model.getVarByName('z_bus['+str(b)+']').Start = TopologyFix['bus'][b]
            counter+=1
        else:
            model.getVarByName('z_bus['+str(b)+']').Start = 1


    for l,i,j in Lines:
        if (l,i,j) in TopologyFix['l_i'].keys():
            model.addConstr(  (model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')  == TopologyFix['l_i'][l,i,j]   )  ,name='eqBendersZ_li')
            model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']').Start = TopologyFix['l_i'][l,i,j]
            counter+=1
        else:
            model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']').Start = 0


    for g in G:
        if g in TopologyFix['g'].keys():
            model.addConstr(  (model.getVarByName('z_g['+str(g)+']') == TopologyFix['g'][g])    ,name='eqBendersZ_g')
            model.getVarByName('z_g['+str(g)+']').Start = TopologyFix['g'][g] 
            counter+=1
        else:
            model.getVarByName('z_g['+str(g)+']').Start = 0

    for d in DemandSet:
        if d in TopologyFix['d'].keys():
            model.addConstr(  (model.getVarByName('z_d['+str(d)+']') == TopologyFix['d'][d])    ,name='eqBendersZ_d')
            model.getVarByName('z_d['+str(d)+']').Start = TopologyFix['d'][d]
            counter+=1
        else:
            model.getVarByName('z_d['+str(d)+']').Start = 0

    if print_result==True:
        print('counter is: ',counter)


    # Limit the number of additional topology changes.
    model.addConstr((  quicksum( (1- model.getVarByName('z_bus['+str(b)+']') ) for b in Bus)
                     + quicksum( model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')   for l,i,j in Lines  ) 
                      + quicksum( model.getVarByName('z_g['+str(g)+']') for g in G)
                       +quicksum( model.getVarByName('z_d['+str(d)+']') for d in DemandSet)
                         <= counter + K_Opt) , name='eq_MaxSw_bus') 


    model.update()
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time 


    # Solve the refinement problem.
    
    status = model.Status

    
    if status == GRB.INFEASIBLE:
        print('\n\nOptimization was stopped with infeasibility!')
        

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
    
    
    # Report results.

    
    ShedCost_df = pd.DataFrame(columns=['ShedCost(c)'],index=cont_list, dtype=float)
    PdShed_dic = pd.DataFrame(columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,busbar,cont_list]), dtype=float)
    for c in cont_list:
        ShedCost_df.loc[c]=model.getVarByName('ShedCost['+str(c)+']').x 
        for d in DemandSet:
            for i in busbar:
                PdShed_dic.loc[d,i,c]=model.getVarByName('Pdi_Shed['+str(d)+','+str(i)+','+str(c)+']').x 


    if print_result==True:
        
        ex_time_df=pd.DataFrame(data=ex_time,index=['Time'],columns=['Time'])
        print('Execution Time', ex_time_df)

        
        OF_df = pd.DataFrame(columns=['OF'],data=[model.getVarByName('OF').x])
        GenCost_df = pd.DataFrame(columns=['GenCost'],data=[model.getVarByName('GenCost').x])
        RDCost_df= pd.DataFrame(columns=['RD_Cost'],data=[model.getVarByName('RDCost').x])
        TotalShedCost_df= pd.DataFrame(columns=['TotalShedCost'],data=[model.getVarByName('TotalShedCost').x])

        
        print('\n\n ======================= Results ======================= \n ')
        print('# Objective:',OF_df )
        print('# Gen Cost: ',GenCost_df )
        print('# RD Cost :',RDCost_df)
        

    # Build the topology dictionary.
    
    z_bus_dict={}; z_g_dict={}; z_d_dict={}; z_li_dict={}
    
    
    for b in Bus:
        z_bus_dict[b]=model.getVarByName('z_bus['+str(b)+']').x
    for g in G:
        z_g_dict[g]=model.getVarByName('z_g['+str(g)+']') .x
    for d in DemandSet:
        z_d_dict[d]=model.getVarByName('z_d['+str(d)+']') .x
    for l,i,j in Lines:
        z_li_dict[(l,i,j)]=model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']').x
    
        
    TopologyDict={'bus':z_bus_dict,
          'g':z_g_dict,
          'd':z_d_dict,
          'l_i':z_li_dict,
         }
    
    TopologyFix_updated = {
        'bus':{},
        'g':{},
        'd':{},
        'l_i':{}
        }
    
    for b in Bus:
        if z_bus_dict[b]==0:
            TopologyFix_updated['bus'][b]=0
    for g in G:
        if z_g_dict[g]==1:
            TopologyFix_updated['g'][g]=1
    for d in DemandSet:
        if z_d_dict[d]==1:
            TopologyFix_updated['d'][d]=1
    for l,i,j in Lines:
        if z_li_dict[(l,i,j)]==1:
            TopologyFix_updated['l_i'][(l,i,j)] =  z_li_dict[(l,i,j)]


    # Extract pre-contingency dispatch.
            
    Pg_dict = {}; Qg_dict = {}
    for g in G:
        Pg_dict[g] = model.getVarByName('Pg['+str(g)+']').x
        Qg_dict[g] =model.getVarByName('Qg['+str(g)+']').x

    Vol0_li = {}
    delta0_li = {}
    for l,i,j in Lines:
        for c in cont_list:
            suffix = '['+str(l)+','+str(i)+','+str(j)+','+str(c)+']'
            Vol0_li[l,i,j,c] = np.sqrt(model.getVarByName('V2_li'+suffix).x)
            delta0_li[l,i,j,c] = model.getVarByName('delta_li'+suffix).x

    
    # Results
    return {
        'time':ex_time,
        'TopologyDict':TopologyDict,
        'TopologyFix_updated':TopologyFix_updated,
        'Pg':Pg_dict,
        'Qg':Qg_dict,
        'MIPGap':model.MIPGap*100,
        'Cost_tot':model.getVarByName('OF').x,
        'GenCost':model.getVarByName('GenCost').x,
        'TotalShedCost':model.getVarByName('TotalShedCost').x,
        'ShedCost_df':ShedCost_df,
        'PdShed_dic':PdShed_dic,
        'Vol0_li':Vol0_li,
        'delta0_li':delta0_li
        }
    
    
def SC_SR_1OptH(
    data,
    cont_list=None,
    K_Opt=1,
    Max_iteration=10000,
    Max_Sw_bus=0,
    Up_redispatch=0,
    Dn_redispatch=1.0,
    Zfixdict=None,
    Zinitial=None,
    PQgFix=None,
    PgFix=None,
    FixedCost=False,
    Alpha=0,
    Pg_market=None,
    Probabilistic=False,
    SolverTime=600,
    Threads=None,
    line_shedding=True,
    Vol0_li=None,
    delta0_li=None,
    print_result=False,
):
    """Run the iterative topology-fixing heuristic to convergence.

    Args:
        data: Network and model-parameter dictionary returned by
            ``read_data_AC``.
        cont_list: Security contingencies; normal operation is added internally.
        K_Opt: Maximum number of new topology decisions optimized per iteration.
        Max_iteration: Maximum number of topology-refinement iterations.
        Max_Sw_bus: Maximum number of busbar-splitting actions allowed.
        Up_redispatch: Fraction of generator capacity available for upward
            active-power redispatch.
        Dn_redispatch: Fraction of generator capacity available for downward
            active-power redispatch.
        Zfixdict: Optional dictionary of topology variables to fix.
        Zinitial: Optional initial topology used as a Gurobi warm start.
        PQgFix: Optional busbar-level active and reactive generation to fix.
        PgFix: Optional total active-generation values to fix.
        FixedCost: Enforce the market-generation cost limit when ``True``.
        Alpha: Relative margin applied to the market-generation cost limit.
        Pg_market: Reference market dispatch used by the fixed-cost constraint.
        Probabilistic: Weight contingency shedding costs by their probabilities
            when ``True``.
        SolverTime: Gurobi time limit per solve in seconds.
        Threads: Optional number of Gurobi solver threads.
        line_shedding: Line-shedding configuration flag retained from the case
            study interface.
        Vol0_li: Optional branch-voltage linearization points.
        delta0_li: Optional branch-angle linearization points.
        print_result: Print detailed solution information when ``True``.

    Returns:
        A dictionary containing iteration histories, the final cost and
        dispatch, load shedding, and total solution time.
    """
    

    # Include normal operation in the modeled states.
    cont_list = cont_list + [0]

    # Build one reusable model for all refinement iterations.
    model0 = create_iterative_heuristic_model(data=data,cont_list=cont_list,
                                        Max_Sw_bus=Max_Sw_bus,
                                    Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,
                                           Zfixdict=Zfixdict,#fixed topology
                                    Zinitial=Zinitial,
                                   PQgFix=PQgFix, PgFix=PgFix,
                                   FixedCost=FixedCost,Alpha=Alpha,Pg_market=Pg_market,
                                   Probabilistic=Probabilistic,
                                    SolverTime=SolverTime,Threads=Threads,
                                    line_shedding=line_shedding,
                                     Vol0_li=Vol0_li,
                                     delta0_li=delta0_li )
    

    # Initially, no topology decisions are fixed.
    TopologyFix = {
        'bus':{},
        'g':{},
        'd':{},
        'l_i':{}
    }

    UB_k={}; UB=1e6

    GenCost_k={}
    Topology_k={}
    Topology_updated_k={}
    UB_min=UB

    time_iteration = {}
    tot_time=0; 

   
    TopologyFix_updated = TopologyFix.copy()

    # Iterative topology-refinement loop
    for iter in tqdm(range(Max_iteration)):


        res = solve_iterative_heuristic_iteration(data=data,TopologyFix=TopologyFix_updated,model0=model0,cont_list=cont_list,K_Opt=K_Opt,print_result=print_result)


        TopologyFix_updated = res['TopologyFix_updated']
        UB = res['Cost_tot']
        tot_time += res['time']

        Topology_k[iter] = res['TopologyDict']
        Topology_updated_k[iter] = res['TopologyFix_updated']
        UB_k[iter] = UB
        time_iteration[iter] = tot_time
        GenCost_k[iter] = res['GenCost']


        # Stop when the objective no longer improves materially.
        if UB < UB_min*0.999:
            UB_min = UB
            print('iter %i, UB: %.2f'%(iter,UB))
        else:
            print('Terminate! ===> UB is %.2f'%UB)
            break

        if tot_time > 2*24*3600:
            print('time is up at %i, terminate!'%tot_time)


    # Results
    return {
        'time_iteration':time_iteration,
         'TopologyDict_k':Topology_k,
           'UB_k':UB_k,
        'Cost_tot':UB,
        'GenCost':GenCost_k[iter],
        'GenCost_k':GenCost_k,
        'TotalShedCost':res['TotalShedCost'],
        'Pg':res['Pg'],
        'Qg':res['Qg'],
        'ShedCost_df':res['ShedCost_df'],
        'PdShed_dic':res['PdShed_dic'],
        'Vol0_li':res['Vol0_li'],
        'delta0_li':res['delta0_li'],
        'time':tot_time,
    }
