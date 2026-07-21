"""Fixed-topology large-system baseline from the Section IV-A case study.

This specialized implementation evaluates the PEGASE 1354-bus system without
busbar splitting, reducing the model size for the large-network comparison.
"""

# %% [markdown]
# # Busbar Coupler Contingency (BCC) Project
#  
# **The proposed approach with linear AC**
# 
# **MP1: normal dispatch**
# - **FSP1: feasibility of line contingency**
# **MP2: find each substation configuration (parallel)**
# **OSP1: substation-contingency with load shedding.**
# - Problem: load shedding can help when active overloading is happening. However, if the generators have negative Qg, load shedding cannot help and we may encounter infeasibility problems. That's why we allow redisptaching Qg in contingencies.
# 
# - OSP-2:we allow line shedding here!
# 
# 
# 
# 
# 

# %% [markdown]
# # Especially designed for no busbar splitting, PEGASE 1354-bus
# # No MP2, we fix the dispatch and default topology, we evaluate with OSP

# %%
# Libraries
#%reset 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import itertools
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
import networkx as nx
import time
from copy import deepcopy
import random
import sys
from tqdm import tqdm
import pickle

from .original_MIP_model import *  # Shared data preparation and AC-OPF utilities.



# from BCC_Plot_Topology_v1 import plot_topology_v1


# tested with Python 3.7.0 & Gurobi 10.0.3

# %%


# %% [markdown]
# # BCC Feasibility-SP: non-radial line contingencies

# %%
def create_FSP_line_v65(data,cont_list=None):
    
    """This FSP-1 should only include line contingencies, without any load shedding and re-disptach.
    However, other contingency formulations have not been removed.
    Here we only use this function for normal operation or non-radial line contingencies."""
    
    
    #======  data
    

    Bus=data['Bus']    # for b in Bus
    branch=data['branch']
    Lines=data['Lines']
  
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']    
    D2B=data['D2B']                 
    Gen_data=data['Gen_data']
    G=data['G']            
    G2B=data['G2B']   
    vmin = data['vmin']
    vmax = data['vmax']
    beta = data['beta']
    
    
    #======
    
    model=gp.Model('BCC Benders FSP-1')
    
    cont_list=cont_list.copy()   
    #cont_list+=[0]             
        
   
    model.Params.OutputFlag=0

    # branch['limit'] = 0.95*branch['limit']

    
#=============================== Variables ==================================================

    Pg=model.addVars(G,lb=0,vtype=GRB.CONTINUOUS,name='Pg')
    Qg=model.addVars(G,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg') 



    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')         
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')         
    
    
    delta=model.addVars(Bus,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta') 
    V2=model.addVars(Bus,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2') 


    OF_FSP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_FSP')             # obj func var
    
    sp_up=model.addVars(Bus,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sp_up')
    sp_dn=model.addVars(Bus,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sp_dn')
    sq_up=model.addVars(Bus,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sq_up')
    sq_dn=model.addVars(Bus,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='sq_dn')
     

#=============================== Equations ==================================================


    

    
    #== Gen
    

    eqPgmax=model.addConstrs((Pg[g]<= Gen_data.loc[g]['Pmax']  for g in G ),name='Eq_Pgmax')
    eqPgmin=model.addConstrs(( -Pg[g] <= -Gen_data.loc[g]['Pmin']  for g in G),name='Eq_Pgmin')
    
    
    
    #Qg
    eqPgmax=model.addConstrs((Qg[g]<= Gen_data.loc[g]['Qmax']  for g in G ),name='Eq_Pgmax')
    eqPgmin=model.addConstrs(( -Qg[g] <= -Gen_data.loc[g]['Qmin']  for g in G),name='Eq_Pgmin')


    
    ###=== PF equations ===================================================

##===== line Contingency 
    eqPflow1_cont=model.addConstrs( ( Pflow[l,i,j,c]==0   
                            for l,i,j in Lines 
                           for c in cont_list               
                           if c==l ) , name='Eq_Pflow_cont')
    
    eqQflow1_cont=model.addConstrs( ( Qflow[l,i,j,c]==0   
                            for l,i,j in Lines 
                           for c in cont_list               
                           if c==l ) , name='Eq_Qflow_cont')
    


    # tight AC limits
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij1')
    
    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij2')
    
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij3')
     
    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit'] 
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij4')

    
    
    # AC Pij-V equations! 

    model.addConstrs( (  Pflow[l,i,j,c] == 
    0.5*branch.loc[(l,i,j)]['g_ij']*( V2[i,c] - V2[j,c] )
    -branch.loc[(l,i,j)]['b_ij']*(delta[i,c]-delta[j,c] ) 
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')
    
    model.addConstrs( (  Qflow[l,i,j,c] == 
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2[i,c] - V2[j,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta[i,c]-delta[j,c]) 
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')
    
 


    # voltage
    model.addConstrs( ( V2[b,c] <= vmax**2 for b in Bus for c in cont_list) , name='eqvmaxb')
    model.addConstrs( ( V2[b,c] >= vmin**2  for b in Bus for c in cont_list) , name='eqvminb')

    
    eq_delta_ref=model.addConstrs((delta[Bus[0],ll]==0     for ll in cont_list ), name='ref_bus_angle' ) 
    
    
    # Balance 
    eq_balance_busbar1=model.addConstrs((                                        
                quicksum( Pg[g]   for g,i in G2B.select('*',b)  )
                -quicksum( Pdemand.loc[d]['Pd']  for d,i in D2B.select('*',b)   ) + sp_up[b,ll] - sp_dn[b,ll] ==
                quicksum(Pflow[l,i,j,ll] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus for ll in cont_list ), name='eq_balance1')
    
    
    
    eq_Qbalance_busbar1=model.addConstrs((                                        
                quicksum( Qg[g]   for g,b in G2B.select('*',b)  )
                -quicksum( Pdemand.loc[d]['Qd']   for d,b in D2B.select('*',b)   ) + sq_up[b,ll] - sq_dn[b,ll] ==
                quicksum(Qflow[l,i,j,ll] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus for ll in cont_list  ), name='eq_balance1')
    
    
    model.addConstr(OF_FSP ==  quicksum( sp_up[b,c]+sp_dn[b,c]+sq_up[b,c]+sq_dn[b,c]    
                                        for b in Bus  for c in cont_list)  ,name='Eq_OF_FSP')


    model.setObjective(OF_FSP,GRB.MINIMIZE)

    model.update()
    
        
    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs
    
        
    return {'model':model,
            'NumBinVars':NumBinVars,
            'NumConVars':NumConVars,
            'NumConstrs':NumConstrs
            }
    
    

# %%
def solve_FSP_line_v65(data,TopologyMP,model0,cont_list,
                                    print_result=False):
    
    """This FSP-1 should only include line contingencies, without any load shedding and re-disptach.
    However, other contingency formulations have not been removed.
    Here we only use this function for normal operation or non-radial line contingencies."""
    
    
    #======  data
    

#     Max_timelimit=data['Max_timelimit']=600 #900
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    Lines=data['Lines']
  
    DemandSet=data['Demandset']    
    G=data['G']            

    #======
    
    model=model0.copy()
     

#=============================== Equations ==================================================


    model.addConstrs(  (model.getVarByName('Pg['+str(g)+']') == TopologyMP['Pg'][g]     for g in G )  ,name='eqBendersPg')
    model.addConstrs(  (model.getVarByName('Qg['+str(g)+']') == TopologyMP['Qg'][g]     for g in G )  ,name='eqBendersQg')
    
    model.update()
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time  #execution time
    
    

    #====== 
    status = model.Status
        
    
    if status == GRB.INFEASIBLE:
        print('\n\nFSP Optimization was stopped with infeasibility!')
        
#         display(TopologyMP)
        

        # Relax the bounds and try to make the model feasible
        print('\n\nThe model is infeasible; relaxing the bounds\n\n')
        orignumvars = model.NumVars
        # relaxing only variable bounds
        model.feasRelaxS(0, False, True, False)
        # for relaxing variable bounds and constraint bounds use
        # model.feasRelaxS(0, False, True, True)

        model.optimize()

        status = model.Status
        if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
                print('The relaxed model cannot be solved \
                       because it is infeasible or unbounded')
        #sys.exit(1)
        if status != GRB.OPTIMAL:
            print('Optimization was stopped with status %d' % status)
#             sys.exit(1)

        # print the values of the artificial variables of the relaxation
        print('\nSlack values:')
        slacks = model.getVars()[orignumvars:]
        for sv in slacks:
            if sv.X > 1e-9:
                print('%s = %g' % (sv.VarName, sv.X))
    
    
    
    #=================================== Obtaining dual values! ============================================     
    
    MuPg={}; MuQg={}; 
       
    for g in G:
        MuPg[g]=model.getConstrByName(str('eqBendersPg['+g+']')).pi
        MuQg[g]=model.getConstrByName(str('eqBendersQg['+g+']')).pi

    
    Mu={
        'Pg':MuPg,
        'Qg':MuQg
    }
    
    # for l,i,j in Lines:
    #     for c in cont_list:
    #         print(l,i,j,c,' : %.2f'%model.getVarByName('Pflow['+str(l)+','+str(i)+','+str(j)+','+str(c)+']').x     )


    if print_result==True:
        
        print(model.getVarByName('OF_FSP').x)
 
        
    return { 
        'Mu':Mu,
        'time':ex_time,
        #'z_lineZc':z_lineZc_df,
        'OF_FSP': model.getVarByName('OF_FSP').x         
    }
    
    
    

# %%


# %% [markdown]
# # OSP substation

# %%
def create_OSP_substation_v65(data,substation,cont_list=None,
                                    Up_redispatch=0,Dn_redispatch=1.0,
                           Max_Sw_bus=0):
    
    """OSP for a fixed topology and a given contingency set.
    OSP cost comes from load shedding.
    con_list should be [sub,sub-1,sub-2]
    
    """
    
    #======  data
    
    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
#     Max_timelimit=data['Max_timelimit']=600 #900
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

    
    
    #======
    
    model=gp.Model('OSP-v65')
    
    cont_list=cont_list.copy() 

    model.Params.OutputFlag=0

    Bus_sub = [substation]
    Bus_non_sub =[b for b in Bus if b!= substation]

    G_sub = [g for g,x in G2B.select('*',substation)]
    G_non_sub = [g for g in G if g not in G_sub ]

    D_sub = [d for d,x in D2B.select('*',substation)]
    D_non_sub = [d for d in DemandSet if d not in D_sub]

    Lines_sub = [(l,i,j) for l,i,j in Lines.select('*',substation,'*') ]
    Lines_non_sub = [(l,i,j) for l,i,j in Lines if (l,i,j) not in Lines_sub ]



    
    Pgi=model.addVars(G_sub,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pgi') 
    dPgi_up=model.addVars(G_sub,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_up')
    dPgi_dn=model.addVars(G_sub,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPgi_dn')
    Qgi=model.addVars(G_sub,busbar,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qgi')  
    dQgi_up=model.addVars(G_sub,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_up')
    dQgi_dn=model.addVars(G_sub,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQgi_dn')

    Pg = model.addVars(G_non_sub,lb=0,vtype=GRB.CONTINUOUS,name='Pg')
    dPg_up=model.addVars(G_non_sub,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPg_up')
    dPg_dn=model.addVars(G_non_sub,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dPg_dn')
    Qg=model.addVars(G_non_sub,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg')  
    dQg_up=model.addVars(G_non_sub,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQg_up')
    dQg_dn=model.addVars(G_non_sub,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='dQg_dn')


    Pdi=model.addVars(D_sub,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Pdi')
    Pdi_Shed=model.addVars(D_sub,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Pdi_Shed')
    Qdi=model.addVars(D_sub,busbar,lb=0,vtype=GRB.CONTINUOUS,name='Qdi')
    Qdi_Shed=model.addVars(D_sub,busbar,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Qdi_Shed')

    Pd_Shed=model.addVars(D_non_sub,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Pd_Shed')
    Qd_Shed=model.addVars(D_non_sub,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='Qd_Shed')

    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')         
    Pflow_li=model.addVars(Lines_sub,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_li')
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')         
    Qflow_li=model.addVars(Lines_sub,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_li')
    
    Pflow_bus=model.addVars(Bus_sub,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow_bus')  #From b1 to b2!
    Qflow_bus=model.addVars(Bus_sub,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow_bus')  #From b1 to b2!

    delta_bi=model.addVars(Bus_sub,busbar,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_bi') 
    delta_li=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='delta_li') 
    V2_bi=model.addVars(Bus_sub,busbar,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2_bi') 
    V2_li=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='V2_li')  


    z_bus=model.addVars(Bus_sub,vtype=GRB.CONTINUOUS,name='z_bus')    
    z_li=model.addVars(Lines_sub,vtype=GRB.CONTINUOUS,name='z_li')  
    z_g=model.addVars(G_sub,vtype=GRB.CONTINUOUS,name='z_g') 
    z_d=model.addVars(D_sub,vtype=GRB.CONTINUOUS,name='z_d')
 

    OF_OSP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_OSP')             # obj func var
    #GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost')
    RDCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='RDCost')
    ShedCost=model.addVars(cont_list,lb=0,vtype=GRB.CONTINUOUS,name='ShedCost')
    TotalShedCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalShedCost')

 
        

#=============================== Equations ==================================================
    # Binary Equations ==================================================

    # Reliability switching 
    # eq_2line_busbar2=model.addConstrs(( 2*(1-z_bus[b]) <= 
    #                                quicksum(z_li[l,i,j] for l,i,j in Lines.select('*',b,'*')) 
    #                         for b in Bus_sub ), name='eq_2line_busbar2' ) 
    # eq_2line_busbar1=model.addConstrs(( 2*(1-z_bus[b]) <= 
    #                                quicksum(1-z_li[l,i,j] for l,i,j in Lines.select('*',b,'*')) 
    #                         for b in Bus_sub ), name='eq_2line_busbar1' )
    
    
    # eq_2line_busbar3=model.addConstrs(( z_bus[b] == 1 
    #                         for b in Bus_sub
    #                          if NumberL2B[b]<=3  ), name='eq_2line_busbar3' )
    
    # # max Switching ================
    # eq_MaxSw_bus=model.addConstr((  quicksum( (1-z_bus[b]) for b in Bus_sub) <= Max_Sw_bus) , name='eq_MaxSw_bus')        
    
    # symmetry 
    for b in Bus_sub:
        lmin,b=L2B.select('*',b)[0]
        lmin,i,j=Lines.select(lmin,b,'*')[0]
        eq_symmetry1=model.addConstr((  z_li[lmin,i,j] ==0  ), name='eq_symmetry1')

    # eq_zbus_ref=model.addConstr( (  z_bus[Bus_sub[0]]==1    ), name='eq_zbus_ref')
    model.addConstr( (  z_bus[substation]==1    ), name='eq_zbus_ref')

    
    #== Gen
    eq_Pg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar1']+dPgi_up[g,'busbar1',c]-dPgi_dn[g,'busbar1',c]   
                             for g in G_sub for c in cont_list  ), name='eq_Pg1min')
    eq_Pg1max=model.addConstrs( (    Pgi[g,'busbar1']+dPgi_up[g,'busbar1',c]-dPgi_dn[g,'busbar1',c] <=(1-z_g[g])*Gen_data.loc[g]['Pmax']  
                             for g in G_sub for c in cont_list ), name='eq_Pg1max')
    
    eq_Pg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Pmin'] <= Pgi[g,'busbar2'] +dPgi_up[g,'busbar2',c]-dPgi_dn[g,'busbar2',c]  
                             for g in G_sub for c in cont_list  ), name='eq_Pg2min')
    eq_Pg2max=model.addConstrs( (   Pgi[g,'busbar2']+dPgi_up[g,'busbar2',c]-dPgi_dn[g,'busbar2',c]  <= z_g[g]*Gen_data.loc[g]['Pmax']   
                             for g in G_sub for c in cont_list  ), name='eq_Pg2max')
    
    #Qg
    eq_Qg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c]   
                             for g in G_sub for c in cont_list  ), name='eq_Qg1min')
    eq_Qg1max=model.addConstrs( (    Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c] <=(1-z_g[g])*Gen_data.loc[g]['Qmax']  
                             for g in G_sub for c in cont_list ), name='eq_Qg1max')
    
    eq_Qg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar2'] +dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]  
                             for g in G_sub for c in cont_list  ), name='eq_Qg2min')
    eq_Qg2max=model.addConstrs( (   Qgi[g,'busbar2']+dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]  <= z_g[g]*Gen_data.loc[g]['Qmax']   
                             for g in G_sub for c in cont_list  ), name='eq_Pg2max')

    # =======   Corrective redispatch/ reserve equations =====================================
    
    Eq_dPgUp_res=model.addConstrs(( dPgi_up[g,i,c] <= Gen_data.loc[g]['Pmax']*Up_redispatch  
                          for g in G_sub for i in busbar for c in cont_list ),name='Eq_dPgUp_res')
    Eq_dPgDn_res=model.addConstrs(( dPgi_dn[g,i,c] <= Gen_data.loc[g]['Pmax']*Dn_redispatch  
                            for g in G_sub for i in busbar for c in cont_list ),name='Eq_dPgDn_res')
    # model.addConstrs(( dPgi_up[g,i,c] == 0  
    #                       for g in G_sub for i in busbar for c in cont_list if c==0 ),name='Eq_dPg0')
    # model.addConstrs(( dPgi_dn[g,i,c] == 0  
    #                         for g in G_sub for i in busbar for c in cont_list if c==0 ),name='Eq_dPgDn_res')
    
    Eq_dQgUp_res=model.addConstrs(( dQgi_up[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #*Up_redispatch  
                          for g in G_sub for i in busbar for c in cont_list ),name='Eq_dQgUp_res')
    Eq_dQgDn_res=model.addConstrs(( dQgi_dn[g,i,c] <= Gen_data.loc[g]['Qmax']*1 #Dn_redispatch  
                            for g in G_sub for i in busbar for c in cont_list ),name='Eq_dQgDn_res')
    # model.addConstrs(( dQgi_up[g,i,c] == 0  
    #                       for g in G_sub for i in busbar for c in cont_list if c==0 ),name='Eq_dQg0')
    # model.addConstrs(( dQgi_dn[g,i,c] == 0  
    #                         for g in G_sub for i in busbar for c in cont_list if c==0 ),name='Eq_dQgDn_res')

    

    # non sub 
    model.addConstrs( (    -Pg[g ]-dPg_up[g ,c]+dPg_dn[g ,c] <= -Gen_data.loc[g]['Pmin']      
                             for g in G_non_sub for c in cont_list  ), name='eq_Pgmin')
    model.addConstrs( (     Pg[g]+dPg_up[g,c]-dPg_dn[g ,c] <= Gen_data.loc[g]['Pmax']  
                             for g in G_non_sub for c in cont_list ), name='eq_Pgmax')
    
    
    model.addConstrs( (    -Qg[g ]-dQg_up[g ,c]+dQg_dn[g ,c]   <= -Gen_data.loc[g]['Qmin']   
                             for g in G_non_sub for c in cont_list  ), name='eq_Qg1min')
    model.addConstrs( (    Qg[g ]+dQg_up[g ,c]-dQg_dn[g ,c] <= Gen_data.loc[g]['Qmax']  
                             for g in G_non_sub for c in cont_list ), name='eq_Qg1max')
    
    
    model.addConstrs(( dPg_up[g,c] <= Gen_data.loc[g]['Pmax']*Up_redispatch  
                          for g in G_non_sub for c in cont_list ),name='Eq_dPgUp_res')
    model.addConstrs(( dPg_dn[g,c] <= Gen_data.loc[g]['Pmax']*Dn_redispatch  
                            for g in G_non_sub for c in cont_list ),name='Eq_dPgDn_res')
    # model.addConstrs(( dPg_up[g,c] == 0  
    #                       for g in G_non_sub for c in cont_list if c==0 ),name='Eq_dPg0')
    # model.addConstrs(( dPg_dn[g,c] == 0  
    #                         for g in G_non_sub for c in cont_list if c==0 ),name='Eq_dPgDn_res')
    
    model.addConstrs(( dQg_up[g,c] <= Gen_data.loc[g]['Qmax']*1 #*Up_redispatch  
                          for g in G_non_sub for c in cont_list ),name='Eq_dQgUp_res')
    model.addConstrs(( dQg_dn[g,c] <= Gen_data.loc[g]['Qmax']*1 #Dn_redispatch  
                            for g in G_non_sub for c in cont_list ),name='Eq_dQgDn_res')
    # model.addConstrs(( dQg_up[g,c] == 0  
    #                       for g in G_non_sub for c in cont_list if c==0 ),name='Eq_dQg0')
    # model.addConstrs(( dQg_dn[g,c] == 0  
    #                         for g in G_non_sub for c in cont_list if c==0 ),name='Eq_dQgDn_res')




    #Pd
    
    eq_Pd1=model.addConstrs( ( Pdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Pd']   
                            for d in D_sub ), name='eq_Pd1')
    eq_Pd2=model.addConstrs( ( Pdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Pd']   
                            for d in D_sub ), name='eq_Pd2')
    eq_Qd1=model.addConstrs( ( Qdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Qd']   
                            for d in D_sub ), name='eq_Qd1')
    eq_Qd2=model.addConstrs( ( Qdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Qd']   
                            for d in D_sub ), name='eq_Qd2')
    
    # Load sheddding
    model.addConstrs( (  Pdi_Shed[d,'busbar1',c]==Pdi[d,'busbar1'] 
                              for b in Bus_sub for d,b in D2B.select('*',b) 
                               for c in cont_list if c==str(b+'-1') ) , name='eq_Pd1Shed')
    model.addConstrs( (  Pdi_Shed[d,'busbar2',c]==Pdi[d,'busbar2'] 
                              for b in Bus_sub for d,b in D2B.select('*',b) 
                               for c in cont_list if c==str(b+'-2') ) , name='eq_Pd2Shed')
    model.addConstrs( (  Pdi_Shed[d,i,c]<=Pdi[d,i] 
                              for d in D_sub for i in busbar 
                               for c in cont_list) , name='eq_PdShedlimit')
    # model.addConstrs( (  Pdi_Shed[d,i,c]==0 
    #                           for d in D_sub for i in busbar 
    #                            for c in cont_list
    #                   if c==0) , name='eq_PdShedlimit0')
    
    model.addConstrs( (  Qdi_Shed[d,i,c]== (Pdemand.loc[d]['Qd'])/(Pdemand.loc[d]['Pd'])*Pdi_Shed[d,i,c]
                              for d in D_sub for i in busbar 
                               for c in cont_list) , name='eq_QdShed')

    # non_sub shedding
    model.addConstrs( (  Pd_Shed[d,c]<= Pdemand.loc[d]['Pd'] 
                                for d in D_non_sub 
                                for c in cont_list) , name='eq_PdShedlimit')
    # model.addConstrs( (  Pd_Shed[d,c]==0 
    #                             for d in D_non_sub 
    #                             for c in cont_list
    #                     if c==0) , name='eq_PdShedlimit0')

    model.addConstrs( (  Qd_Shed[d,c]== (Pdemand.loc[d]['Qd'])/(Pdemand.loc[d]['Pd'])*Pd_Shed[d,c]
                                for d in D_non_sub  
                                for c in cont_list) , name='eq_QdShed')
        
    
    ###=== PF equations ===================================================

##===== line Contingency 
    eqPflow1_cont=model.addConstrs( ( Pflow[l,i,j,c]==0   
                            for l,i,j in Lines 
                           for c in cont_list               
                           if c==l ) , name='Eq_Pflow_cont')
    eqPflow_ei_cont=model.addConstrs( ( Pflow_li[l,i,j,m,c]==0   
                            for l,i,j in Lines_sub 
                            for m in busbar
                           for c in cont_list                     
                           if c==l ) , name='Eq_Pflow_ei_cont')
    eqQflow1_cont=model.addConstrs( ( Qflow[l,i,j,c]==0   
                            for l,i,j in Lines 
                           for c in cont_list               
                           if c==l ) , name='Eq_Qflow_cont')
    eqQflow_ei_cont=model.addConstrs( ( Qflow_li[l,i,j,m,c]==0   
                            for l,i,j in Lines_sub 
                            for m in busbar
                           for c in cont_list                     
                           if c==l ) , name='Eq_Qflow_ei_cont')

    ##========== flow limits   
    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar1',c] 
                             for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Pflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Pflow_li[l,i,j,'busbar2',c] 
                            for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow2min')
    
    eq_flow2max=model.addConstrs( ( Pflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow2max')
    
    eq_flow=model.addConstrs((   Pflow[l,i,j,c] == Pflow_li[l,i,j,'busbar1',c]+Pflow_li[l,i,j,'busbar2',c]  
                                for l,i,j in Lines_sub for c in cont_list), name='eq_Pflow')
    
    
    eq_flow1min=model.addConstrs( ( -(1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar1',c] 
                             for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow1min')
    eq_flow1max=model.addConstrs( ( Qflow_li[l,i,j,'busbar1',c] <= (1-z_li[l,i,j])*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow1max')
    eq_flow2min=model.addConstrs( ( -z_li[l,i,j]*branch.loc[(l,i,j)]['limit'] <= Qflow_li[l,i,j,'busbar2',c] 
                            for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow2min')
    eq_flow2max=model.addConstrs( ( Qflow_li[l,i,j,'busbar2',c] <= z_li[l,i,j]*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines_sub  for c in cont_list ) ,name='eq_flow2max')
    
    eq_flow=model.addConstrs((   Qflow[l,i,j,c] == Qflow_li[l,i,j,'busbar1',c]+Qflow_li[l,i,j,'busbar2',c]  
                                for l,i,j in Lines_sub for c in cont_list), name='eq_Qflow')


    
    # tight AC limits
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij1')
    
    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij2')
    
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij3')
     
    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit'] 
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij4')

    
    
    # AC Pij-V equations! 

    model.addConstrs( (  Pflow[l,i,j,c] == 
    0.5*branch.loc[(l,i,j)]['g_ij']*( V2_li[l,i,j,c] - V2_li[l,j,i,c] )
    -branch.loc[(l,i,j)]['b_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) #+ Ploss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')
    
    # if c!=str(i+'-1') if c!=str(i+'-2') if c!=str(j+'-1') if c!=str(j+'-2')
    
    model.addConstrs( (  Qflow[l,i,j,c] == 
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j,c] - V2_li[l,j,i,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) 
#     -V2_li[l,i,j,c]*(branch.loc[(l,i,j)]['b']/2) #+ Qloss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c  ) , name='eqQij')
    

    
    for b in Bus_non_sub:
        lmin,b,j1=Lines.select('*',b,'*')[0]
        model.addConstrs( (  delta_li[lmin,b,j1,c] == delta_li[l,i,j,c]
                           for l,i,j in Lines.select('*',b,'*') for c in cont_list  ) , name='eqDelta')
        
        model.addConstrs( (  V2_li[lmin,b,j1,c] == V2_li[l,i,j,c]
                           for l,i,j in Lines.select('*',b,'*') for c in cont_list  ) , name='eqV2')



    
    


    # voltage    
    model.addConstrs( ( V2_bi[b,i,c] <= vmax**2  for b in Bus_sub for i in busbar for c in cont_list) , name='eqvmaxb')
    model.addConstrs( ( V2_bi[b,i,c] >= vmin**2  for b in Bus_sub for i in busbar for c in cont_list) , name='eqvminb')
    model.addConstrs( ( V2_li[l,i,j,c] <= vmax**2 for l,i,j in Lines for c in cont_list) , name='eqvmaxl')
    model.addConstrs( ( V2_li[l,i,j,c] >= vmin**2  for l,i,j in Lines for c in cont_list) , name='eqvminl')
    
    eq_delta_bus1=model.addConstrs((   -BigM_b*(1-z_bus[b]) <= delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll] 
                               for b in Bus_sub for ll in cont_list 
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_delta_bus1')
    eq_delta_bus2=model.addConstrs((   delta_bi[b,'busbar1',ll]-delta_bi[b,'busbar2',ll] <= BigM_b*(1-z_bus[b])
                               for b in Bus_sub for ll in cont_list 
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_delta_bus2')

    model.addConstrs(( -z_li[l,i,j]*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar1',c]   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar1',c] <= z_li[l,i,j]*Maxdelta   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*Maxdelta<= delta_li[l,i,j,c] - delta_bi[i,'busbar2',c]   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmin')
    model.addConstrs((  delta_li[l,i,j,c] - delta_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*Maxdelta   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_delta_lb2Frmax')
    
    #V2
    eq_V2_bus1=model.addConstrs((   -MaxV2*(1-z_bus[b]) <= V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll] 
                               for b in Bus_sub for ll in cont_list 
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll), name='eq_V2_bus1')
    eq_V2_bus2=model.addConstrs((   V2_bi[b,'busbar1',ll]-V2_bi[b,'busbar2',ll] <= MaxV2*(1-z_bus[b])
                               for b in Bus_sub for ll in cont_list 
                                    if b!=ll if str(b+'-1')!=ll  if str(b+'-2')!=ll ), name='eq_V2_bus2')
    
    model.addConstrs(( -z_li[l,i,j]*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar1',c]   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar1',c] <= z_li[l,i,j]*MaxV2   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb1Frmax')
    model.addConstrs(( -(1-z_li[l,i,j])*MaxV2<= V2_li[l,i,j,c] - V2_bi[i,'busbar2',c]   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmin')
    model.addConstrs((  V2_li[l,i,j,c] - V2_bi[i,'busbar2',c] <= (1-z_li[l,i,j])*MaxV2   
                                for l,i,j in Lines_sub
                                   for c in cont_list
                                    if c!=l), name='eq_V2_lb2Frmax')

    
    # eq_delta_ref=model.addConstrs((delta_bi[Bus_sub[0],'busbar1',ll]==0     for ll in cont_list ), name='ref_bus_angle' )
    eq_delta_ref=model.addConstrs((delta_li[Lines[0][0],Lines[0][1],Lines[0][2] ,ll]==0     for ll in cont_list ), name='ref_bus_angle' ) 
 


    #busbar & coupler contingency
    model.addConstrs((  Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])     
                          for b in Bus_sub for ll in cont_list if b!=ll), name='eq_busbarPflowmax')
    model.addConstrs((  -Pflow_bus[b,ll]<=BigM_busbar*(z_bus[b])     
                          for b in Bus_sub for ll in cont_list if b!=ll), name='eq_busbarPflowmin')
    
    model.addConstrs((  Pflow_bus[b,ll]==0     
                          for b in Bus_sub for ll in cont_list if ll==b ), name='eq_CouplerCont')
    model.addConstrs((  Pflow_bus[b,ll]==0     
                          for b in Bus_sub for ll in cont_list if ll==str(b+'-1') ), name='eq_bus1Cont')
    model.addConstrs((  Pflow_bus[b,ll]==0     
                          for b in Bus_sub for ll in cont_list if ll==str(b+'-2') ), name='eq_bus2Cont')
    
    model.addConstrs( ( Pflow_li[l,i,j,'busbar1',c]==0
                            for b in Bus_sub for l,i,j in Lines.select('*',b,'*') 
                           for c in cont_list                     
                           if c==str(b+'-1') ) , name='Eq_Pflow_ei_bus1contFr')   
    model.addConstrs( ( Pflow_li[l,i,j,'busbar2',c]==0
                            for b in Bus_sub for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list                     
                           if c==str(b+'-2') ) , name='Eq_Pflow_ei_bus2contFr')
    
    model.addConstrs((  Qflow_bus[b,ll]<=BigM_busbar*(z_bus[b])     
                          for b in Bus_sub for ll in cont_list if b!=ll), name='eq_busbarQflowmax')
    model.addConstrs((  -Qflow_bus[b,ll]<=BigM_busbar*(z_bus[b])     
                          for b in Bus_sub for ll in cont_list if b!=ll), name='eq_busbarQflowmin')
    
    model.addConstrs((  Qflow_bus[b,ll]==0     
                          for b in Bus_sub for ll in cont_list if ll==b ), name='eq_CouplerCont')
    model.addConstrs((  Qflow_bus[b,ll]==0     
                          for b in Bus_sub for ll in cont_list if ll==str(b+'-1') ), name='eq_bus1Cont')
    model.addConstrs((  Qflow_bus[b,ll]==0     
                          for b in Bus_sub for ll in cont_list if ll==str(b+'-2') ), name='eq_bus2Cont')
    
    model.addConstrs( ( Qflow_li[l,i,j,'busbar1',c]==0
                            for b in Bus_sub for l,i,j in Lines.select('*',b,'*') 
                           for c in cont_list                     
                           if c==str(b+'-1') ) , name='Eq_Qflow_ei_bus1contFr')   
    model.addConstrs( ( Qflow_li[l,i,j,'busbar2',c]==0
                            for b in Bus_sub for l,i,j in Lines.select('*',b,'*')
                           for c in cont_list                     
                           if c==str(b+'-2') ) , name='Eq_Qflow_ei_bus2contFr')
    


    ## ==== Balance ================ 
    
    eq_balance_busbar1=model.addConstrs((                                        
                quicksum( Pgi[g,m] +dPgi_up[g,m,ll]-dPgi_dn[g,m,ll]  for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m] - Pdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Pflow_bus[b,ll]
                            for b in Bus_sub for m in busbar for ll in cont_list 
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_balance_busbar2=model.addConstrs((                                        
                quicksum( Pgi[g,m]+dPgi_up[g,m,ll]-dPgi_dn[g,m,ll] for g,b in G2B.select('*',b)  )
                -quicksum( Pdi[d,m]- Pdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Pflow_bus[b,ll]
                            for b in Bus_sub for m in busbar for ll in cont_list 
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')
    
    
    eq_Qbalance_busbar1=model.addConstrs((                                        
                quicksum( Qgi[g,m] +dQgi_up[g,m,ll]-dQgi_dn[g,m,ll]  for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m] - Qdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                +Qflow_bus[b,ll]
                            for b in Bus_sub for m in busbar for ll in cont_list 
                        if m=='busbar1'  if ll!=str(b+'-1') ), name='eq_balance1')
    eq_Qbalance_busbar2=model.addConstrs((                                        
                quicksum( Qgi[g,m]+dQgi_up[g,m,ll]-dQgi_dn[g,m,ll] for g,b in G2B.select('*',b)  )
                -quicksum( Qdi[d,m]- Qdi_Shed[d,m,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow_li[l,i,j,m,ll] for l,i,j in Lines.select('*',b,'*'))
                -Qflow_bus[b,ll]
                            for b in Bus_sub for m in busbar for ll in cont_list 
                            if m=='busbar2'  if ll!=str(b+'-2') ), name='eq_balance2')
    



    #balance for non-sub
    model.addConstrs((                                        
                quicksum( Pg[g] +dPg_up[g,ll]-dPg_dn[g,ll]  for g,b in G2B.select('*',b)  )
                - quicksum( Pdemand.loc[d]['Pd'] - Pd_Shed[d,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Pflow[l,i,j,ll] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus_non_sub for ll in cont_list ), name='eq_balance')
    
    model.addConstrs((                                        
                quicksum( Qg[g] +dQg_up[g,ll]-dQg_dn[g,ll]  for g,b in G2B.select('*',b)  )
                - quicksum( Pdemand.loc[d]['Qd'] - Qd_Shed[d,ll]  for d,b in D2B.select('*',b)   ) ==
                quicksum(Qflow[l,i,j,ll] for l,i,j in Lines.select('*',b,'*'))
                            for b in Bus_non_sub for ll in cont_list ), name='eq_balance')

    
    
    
    #cost
    model.addConstr(RDCost == quicksum( Gen_data.loc[g]['c_res']*(dPgi_up[g,i,c]+0.0001*dPgi_dn[g,i,c]) 
                                       for g in G_sub for i in busbar for c in cont_list)
                                + quicksum( Gen_data.loc[g]['c_res']*(dPg_up[g,c]+0.0001*dPg_dn[g,c]) 
                                       for g in G_non_sub for c in cont_list)
                              ,name='Eq_RD')
    
    model.addConstrs( (ShedCost[c] == quicksum(Pdemand.loc[d]['ShedCost']*Pdi_Shed[d,i,c] 
                                              for d in D_sub for i in busbar)
                                        + quicksum(Pdemand.loc[d]['ShedCost']*Pd_Shed[d,c] 
                                              for d in D_non_sub )     for c in cont_list) ,name='Eq_ShedC')
    
    model.addConstr(TotalShedCost ==  quicksum(ShedCost[c]  for c in cont_list)   ,name='Eq_TotalShedCost')
    
    model.addConstr(OF_OSP == RDCost + TotalShedCost ,name='Eq_OF_OSP')


    model.setObjective(OF_OSP,GRB.MINIMIZE)

    model.update()

    # model.write('OSP_'+str(substation)+'.lp')
    
    
    
    
    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs
    
        
    return {'model':model,
            'NumBinVars':NumBinVars,
            'NumConVars':NumConVars,
            'NumConstrs':NumConstrs
            }
    
    
    

# %%


# %%
def solve_OSP_substation_v65(data,substation,TopologyMP,model0,cont_list=None):
    
    """OSP for a fixed topology and a given contingency set.
    OSP cost comes from load shedding."""
    
    #======  data
    
#     Max_timelimit=data['Max_timelimit']=600 #900
    Bus=data['Bus']    # for b in Bus
    busbar=data['busbar']
    Lines=data['Lines']
    
    DemandSet=data['Demandset']    
    D2B=data['D2B']                 
    G=data['G']            
    G2B=data['G2B']   
    

    
    #======
    
    model=model0.copy()

    Bus_sub = [substation]
    Bus_non_sub =[b for b in Bus if b!= substation]

    G_sub = [g for g,x in G2B.select('*',substation)]
    G_non_sub = [g for g in G if g not in G_sub ]

    D_sub = [d for d,x in D2B.select('*',substation)]
    D_non_sub = [d for d in DemandSet if d not in D_sub]

    Lines_sub = [(l,i,j) for l,i,j in Lines.select('*',substation,'*') ]
    Lines_non_sub = [(l,i,j) for l,i,j in Lines if (l,i,j) not in Lines_sub ]


    

#=============================== Equations ==================================================


    model.addConstrs(  (model.getVarByName('Pgi['+str(g)+','+str(i)+']') == TopologyMP['Pgi'][g,i]     for g in G_sub for i in busbar )  ,name='eqBendersPgi')

    model.addConstrs(  (model.getVarByName('Qgi['+str(g)+','+str(i)+']') == TopologyMP['Qgi'][g,i]     for g in G_sub for i in busbar )  ,name='eqBendersQgi')


    model.addConstrs( (   model.getVarByName('Pg['+str(g)+']')   == TopologyMP['Pg'][g]     for g in G_non_sub ), name='eqBendersPg')
    model.addConstrs( (   model.getVarByName('Qg['+str(g)+']')   == TopologyMP['Qg'][g]     for g in G_non_sub ), name='eqBendersQg')
        
    model.addConstrs(  (model.getVarByName('z_bus['+str(b)+']') == (TopologyMP['bus'][b])      for b in Bus_sub )  ,name='eqBendersZ_bus')
        
    model.addConstrs(  (model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')  == (TopologyMP['l_i'][l,i,j])      for l,i,j in Lines_sub )  ,name='eqBendersZ_li')
    
    model.addConstrs(  (model.getVarByName('z_g['+str(g)+']') == (TopologyMP['g'][g])      for g in G_sub )  ,name='eqBendersZ_g')
    
    model.addConstrs(  (model.getVarByName('z_d['+str(d)+']') == (TopologyMP['d'][d])     for d in D_sub )  ,name='eqBendersZ_d')
    

    model.update()
    #model.write('BCC_Benders_v4_OSP.lp')
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time 
    

#     status = model.Status

    
#     if status == GRB.INFEASIBLE:
#         print('\n\nMP2 was stopped with infeasibility!')
#         print('cont list: ',cont_list)


#         # Relax the bounds and try to make the model feasible
#         print('\n\nThe model is infeasible; relaxing the bounds\n\n')
#         orignumvars = model.NumVars
#         # relaxing only variable bounds
# #         model.feasRelaxS(0, False, True, False)
#         # for relaxing variable bounds and constraint bounds use
# #         model.feasRelaxS(0, False, True, True)
#         model.feasRelaxS(0, False, False, True) #relaxing constraints
#         model.optimize()
#         status = model.Status
#         if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
#                 print('The relaxed model cannot be solved \
#                        because it is infeasible or unbounded')
#         if status != GRB.OPTIMAL:
#             print('Optimization was stopped with status %d' % status)

#         # print the values of the artificial variables of the relaxation
#         print('\nSlack values:')
#         slacks = model.getVars()[orignumvars:]
#         for sv in slacks:
#             if sv.X > 1e-9:
#                 print('%s = %g' % (sv.VarName, sv.X))


    #=================================== Obtaining dual values! ============================================ 
    
    
    MuPgi={}; MuQgi={}; Muz_bus={}; Muz_li={}; Muz_g={}; Muz_d={}
       
    for g in G_sub:
        for m in busbar:
            MuPgi[(g,m)]=model.getConstrByName(str('eqBendersPgi['+g+','+m+']')).pi
            MuQgi[(g,m)]=model.getConstrByName(str('eqBendersQgi['+g+','+m+']')).pi
    for b in Bus_sub:
        Muz_bus[b]=model.getConstrByName(str('eqBendersZ_bus['+b+']')).pi  
    for g in G_sub:
        Muz_g[g]=model.getConstrByName(str('eqBendersZ_g['+g+']')).pi
    for d in D_sub:
        Muz_d[d]=model.getConstrByName(str('eqBendersZ_d['+d+']')).pi
    for l,i,j in Lines_sub:
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
    for g in G_sub:
        if model.getVarByName('z_g['+str(g)+']').x == 0:
            Mu_Pg[g]=Mu['Pgi'][g,'busbar1']
            Mu_Qg[g]=Mu['Qgi'][g,'busbar1']
        elif model.getVarByName('z_g['+str(g)+']').x == 1:
            Mu_Pg[g]=Mu['Pgi'][g,'busbar2']
            Mu_Qg[g]=Mu['Qgi'][g,'busbar2']
    for g in G_non_sub:
        Mu_Pg[g] = model.getConstrByName(str('eqBendersPg['+g+']')).pi
        Mu_Qg[g] = model.getConstrByName(str('eqBendersQg['+g+']')).pi

                
    Mu['Pg']=Mu_Pg
    Mu['Qg']=Mu_Qg
    
    
    
    
    ShedCost_df = pd.DataFrame(0,columns=['ShedCost(c)'],index=cont_list, dtype=float)
    PdShed_dc = pd.DataFrame(0,columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,cont_list]), dtype=float)
    for c in cont_list:
        ShedCost_df.loc[c]=model.getVarByName('ShedCost['+str(c)+']').x 
        for d in D_sub:
            PdShed_dc.loc[d,c]=model.getVarByName('Pdi_Shed['+str(d)+',busbar1,'+str(c)+']').x+model.getVarByName('Pdi_Shed['+str(d)+',busbar2,'+str(c)+']').x
        for d in D_non_sub:
            PdShed_dc.loc[d,c]=model.getVarByName('Pd_Shed['+str(d)+','+str(c)+']').x





        
    
    return { 
        'Mu':Mu,
        'time':ex_time,
        #'z_lineZc':z_lineZc_df,
        'OF_OSP': model.getVarByName('OF_OSP').x, 
        'ShedCost_df':ShedCost_df,
                        'PdShed_dc':PdShed_dc
    }
    
    
    

# %%


# %%


# %% [markdown]
# # Main function

# %%
def BCC_1354_full_AC(data,line_cont_list=[],
                            Max_iter=10,
                              Max_FSP_iter=10,FSP_criteria=0,
                               Pg_market_fix=None,
                            Max_Sw_bus=0,Up_redispatch=0,Dn_redispatch=1.0):
    
    #======  data
    
    Sbase=data['Sbase']
    Max_MIPGap=data['Max_MIPGap']
#     Max_timelimit=data['Max_timelimit']=600 #900
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
    

    #======
    
    model=gp.Model('BCC_main_SCOPF')
    model.Params.OutputFlag=0

    if Max_Sw_bus!=0:
        print('\n no switching here!!\n')

#=============================== Variables ==================================================

    Pg=model.addVars(G,lb=0,vtype=GRB.CONTINUOUS,name='Pg')
    Qg=model.addVars(G,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg')
    
    GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost')
    OF_MP=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF_MP') 
    
    Phi_b=model.addVars(Bus,lb=0,vtype=GRB.CONTINUOUS,name='Phi_b')
    
    line_cont_radial =[]
    line_cont_list0 = line_cont_list+[0]
    for l in line_cont_list:
        if l not in data['line_cont_notradial']:
            line_cont_radial += [l]
            
    Phi_l = model.addVars(line_cont_list0,lb=0,vtype=GRB.CONTINUOUS,name='Phi_l')
    
    TotalPhi=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='TotalPhi')
    
#=============================== Equations ==================================================

    eqTotBalance=model.addConstr( (   quicksum(Pdemand.loc[d]['Pd'] for d in DemandSet) <= quicksum( Pg[g] for g in G) ),
                                 name='EqTotBalance')  

    
    model.addConstrs((Pg[g]<= Gen_data.loc[g]['Pmax']  for g in G ),name='Eq_Pgmax')
    model.addConstrs(( -Pg[g] <= -Gen_data.loc[g]['Pmin']  for g in G),name='Eq_Pgmin')

    model.addConstrs((Qg[g]<= Gen_data.loc[g]['Qmax']  for g in G ),name='Eq_Qgmax')
    model.addConstrs(( -Qg[g] <= -Gen_data.loc[g]['Qmin']  for g in G),name='Eq_Qgmin')

    if Pg_market_fix is not None:
        model.addConstrs(( Pg[g] == Pg_market_fix[g]  for g in G),name='Eq_PFix')
    

    
    eqOF=model.addConstr(GenCost==quicksum(Gen_data.loc[g]['b']*Pg[g] for g in G ) ,name='Eq_OF')
    
    
    model.addConstr( TotalPhi == quicksum(Phi_b[b] for b in Bus) + quicksum(Phi_l[l] for l in line_cont_list0)  
                    ,name='Eq_Phi')
#     model.addConstr(OF_MP == GenCost + TotalPhi   ,name='Eq_OF_MP')
    
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

    # =======================================================================================================================
    
    # MP-1 is solved
        
    PgMP={}  ; QgMP={}
    for g in G:
        PgMP[g]=Pg[g].x
        QgMP[g]=Qg[g].x
        
        
    # topology initialization
    FinalTopology={'bus':{b:1 for b in Bus},
                   'g':{g:0 for g in G}, 'd':{d:0 for d in DemandSet},
                   'l_i':{(l,i,j):0 for l,i,j in Lines}, 
                  'Pgi':{}, 'Qgi':{},
                  'Pg':{} ,'Qg':{}}
    for g in G:
        FinalTopology['Pgi'][(g,'busbar1')]=PgMP[g]; FinalTopology['Pgi'][(g,'busbar2')]=0
        FinalTopology['Qgi'][(g,'busbar1')]=QgMP[g]; FinalTopology['Qgi'][(g,'busbar2')]=0
        FinalTopology['Pg'][g]=PgMP[g]; FinalTopology['Qg'][g]=QgMP[g]
        
    
    
    topology_default0 = {'bus':{b:1 for b in Bus},
                   'g':{g:1 for g in G}, 'd':{d:0 for d in DemandSet},
                   'l_i':{(l,i,j):0 for l,i,j in Lines},
                   'Pgi':{},'Qgi':{},
                   'Pg':{},'Qg':{}}
    
    
    
    
    Mu={}
    UB_k={}; LB_k={}
    UB=1e8; LB=0 


    k=0         #OSP outer loop
    k1=0        #FSP inner loop
    UB_min=UB
    
    NumberofCuts=0
    
    
    all_sub_cont_list=[]
    for b in Bus:
        all_sub_cont_list+=[b,str(b+'-1'),str(b+'-2')]
            
    all_cont=all_sub_cont_list+line_cont_list #used for reporting LoadShedding_df
    
    ShedCost_df = pd.DataFrame(0,columns=['ShedCost(c)'],index=all_cont, dtype=float) #not including c0
    PdShed_dc = pd.DataFrame(0,columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,all_cont]), dtype=float)


    
    #======================================
    
    line_cont_notradial =[0] # in case of MP1=ED
    for c in line_cont_list:
        if c in data['line_cont_notradial']:
            line_cont_notradial +=[c]



    # create FSP and OSP models
    print('creating FSP, OSP, MP2 models.')
    FSP_models={}
    




    for c in line_cont_notradial:
        FSP_models[c] = create_FSP_line_v65(data=data,cont_list=[c])['model']
    
        
    
            
    
    solution_dict = {} #iter:sol

        
    current_time = time.localtime()
    formatted_time = time.strftime("%Y-%m-%d %H:%M:%S", current_time)
    print("Current time and date:", formatted_time)
    
    # MP-2 substation outer loop
    while k < Max_iter:
        
        
        print('\n**Main teration %3i'%(k))
        
        # =================================== FSP inner loop for non-radial lines
        k1=0
        for k1 in tqdm(range(Max_FSP_iter)):
            
            # print('FSP iteration %3i'%k1)
        
            FSP_Objc=0
            FSP_Obj=0

            for g in G: #TopologyMP=def, so we check the FSP for the default top
                topology_default0['Pgi'][(g,'busbar1')]=PgMP[g]; topology_default0['Pgi'][(g,'busbar2')]=0
                topology_default0['Pg'][g]=PgMP[g]
                topology_default0['Qgi'][(g,'busbar1')]=QgMP[g]; topology_default0['Qgi'][(g,'busbar2')]=0
                topology_default0['Qg'][g]=QgMP[g]


            for c in (line_cont_notradial):
#                 print('\t\tcont %s'%(c),end='')
                # result = BCC_AC_FSP_line_v61(data,FinalTopology,cont_list=[c],print_result=False)


                result = solve_FSP_line_v65(data=data,TopologyMP=topology_default0,model0=FSP_models[c],
                                            cont_list=[c],print_result=False) #Problem: should we check default topology for disptach, or the new one?!
                #TopologyMP=FinalTopology
                FSP_Objc=result['OF_FSP']
                FSP_Obj+=FSP_Objc
                Lambda=result['Mu']
                tot_time += result['time']
            
                if FSP_Objc > FSP_criteria:

                    model.addConstr(  ( FSP_Objc   
                                +quicksum( Lambda['Pg'][g]*(Pg[g]-PgMP[g])  for g in G)
                                +quicksum( Lambda['Qg'][g]*(Qg[g]-QgMP[g])  for g in G)       
                                        <= 0   )
                        ,name=str('FSP-'+str(c)+'-'+str(NumberofCuts)))
                    NumberofCuts+=1

                #=========== MP update
            if FSP_Obj > FSP_criteria:
                k1+=1
                model.update()
                start_time = time.time()
                model.optimize()
                end_time = time.time()
                ex_time=end_time-start_time  
                tot_time += ex_time

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

                    
     
                    
            else:
                # print('\n\t\tFSP converged in iteration %i'%k1)
                break
                
            # end FSP1 loop =========================================================================================

    
        
        # =============================================== MP2 with OSP loop
        OF_OSP=0 ; OF_OSP_b={}
        Mu={}

        for sub in tqdm(Bus):

            res = create_OSP_substation_v65(data=data,substation=sub, cont_list=[sub,str(sub+'-1'),str(sub+'-2')]
                                                ,Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,Max_Sw_bus=Max_Sw_bus)
            OSP_model = res['model']
            


            
            result = solve_OSP_substation_v65(data=data,substation=sub, cont_list=[sub,str(sub+'-1'),str(sub+'-2')],TopologyMP=FinalTopology, model0=OSP_model)
            
            OF_OSP_b[sub]=result['OF_OSP']
            OF_OSP+=OF_OSP_b[sub]
            Mu[sub]=result['Mu']
            tot_time += result['time']

            
            for cont in all_cont:
                if cont in result['ShedCost_df'].index:
                    ShedCost_df.loc[cont] = result['ShedCost_df'].loc[cont]
                    for d in DemandSet:
                            PdShed_dc.loc[(d,cont)] = result['PdShed_dc'].loc[(d,cont)]
        
        #============================== MP2 loop done

            
        
        
        
        # ================================ OSP for radial line cont
        
        #OF_OSP has a value here
        OF_OSP_c={}
        
        for c in tqdm(line_cont_list0):
            # print('\tOSP line cont %3s'%c)

            res = create_OSP_substation_v65(data=data,cont_list=[c],substation=Bus[0],Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,Max_Sw_bus=Max_Sw_bus)
            OSP_model = res['model']

            result=solve_OSP_substation_v65(data=data,model0=OSP_model, cont_list=[c],substation=Bus[0],TopologyMP=FinalTopology)
            Mu[c]=result['Mu']

            OF_OSP_c[c]=result['OF_OSP']
            OF_OSP+=OF_OSP_c[c]          #OF_OSP has values before here


            for cont in all_cont:
                if cont in result['ShedCost_df'].index:
                    ShedCost_df.loc[cont] = result['ShedCost_df'].loc[cont]
                    for d in DemandSet:
                            PdShed_dc.loc[(d,cont)] = result['PdShed_dc'].loc[(d,cont)]
        
        
        
        obj = model.getObjective().getValue() #this is OF_MP.x
            
        LB=obj ; LB_k[k]=LB    
        UB=GenCost.x+OF_OSP ; UB_k[k]=UB
        time_iteration[k] = tot_time

        k_min = min(UB_k, key=UB_k.get)
        UB_min = UB_k[k_min]

        solution_dict[k] = {'FinalTopology':FinalTopology,'Pg':PgMP,'Qg':QgMP,
                            'GenCost':GenCost.x,'TotalShedCost':OF_OSP,'ShedCost_df':ShedCost_df,
                            'ShedCost_df':ShedCost_df,'PdShed_dc':PdShed_dc,
                            'Cost_tot':UB,
                            'time':tot_time,'time_iteration':time_iteration}


        
        print('********UB_k=%.2f and LB=%.2f  =>  Gap=%.3f'%(UB,LB,UB-LB))
        print('UB_min: %.2f'%UB_min)

        if Pg_market_fix is not None:
            print('\nPg fixed, only one iteration! UB is %.2f'%(UB))
            time_iteration[k] = tot_time
            break


        
        print('\nMain Gap is %.2f, loop break because we have one iteration!'%(UB-LB))
        time_iteration[k] = tot_time
        break #while main iteration
    


    
    
    # plt.figure(figsize=(5,3))
    # plt.plot(UB_k.keys(),UB_k.values(),marker='o')
    # plt.plot(LB_k.keys(),LB_k.values(),marker='o')
    # plt.legend(['UB','LB'],loc='best')
    # #plt.xticks(range(k+1))
    # plt.title('UB & LB')
    # plt.show()
    
    
    
    
            
#     display(FinalTopology)
            
            
            
    print('Best solution is UB_kmin=%.2f and LB=%.2f  =>  Gap=%.3f'%(UB_k[k_min],LB,UB_k[k_min]-LB))
    
    for b in Bus:
        if solution_dict[k_min]['FinalTopology']['bus'][b]==0:
            print('bus %s is open.'%b)
    
    
    
    return { 'TopologyDict':solution_dict[k_min]['FinalTopology'],
            'Pg':solution_dict[k_min]['Pg'],
            # 'Qg':solution_dict[k_min]['Qg'],
            'UB_k':UB_k,
            'LB_k':LB_k,
            'Cost_tot':solution_dict[k_min]['Cost_tot'],
        'GenCost':solution_dict[k_min]['GenCost'],
        'TotalShedCost':solution_dict[k_min]['TotalShedCost'],
        'ShedCost_df':solution_dict[k_min]['ShedCost_df'],
                        'PdShed_dc':solution_dict[k_min]['PdShed_dc'],
             'time':tot_time,
            'time_iteration':time_iteration,
    }
            
        
        
        
        
        
    
    
    
    
    

# %%




# %%
# sys=14
# lm=1
# dem=1.0
# sw=1
# CaseStudies={}


# data=read_data_AC(File='IEEE_14_bus_Data_PGLib_ACOPF.xlsx',DemFactor=1.0,LineLimit=1,print_data=False)

# line_cont = [] #data['line_cont_notradial']
# lc=len(line_cont)



    
# res = BCC_1354_full_AC(data,line_cont_list=line_cont,Max_Sw_bus=0,
#                                    Max_iter=1,Max_FSP_iter=20,FSP_criteria=0)


# # CaseStudies['full'] = res

# # File_name = f'CaseStudy2_AC_full_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

# # with open(File_name, 'wb') as f:
# #     pickle.dump(CaseStudies, f)


# %%


# %%
# AllCaseStudies = {}

# sys=1354
# dem=1.0
# sw=0
# lm=1
# lc=0

# Models=['full','Proposed'] 

# for model in Models:
#     File_name = f'CaseStudy2_{model}_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'
#     with open(File_name, 'rb') as f:
#         AllCaseStudies[model] = pickle.load(f)[model]
        
# File_name = f'CaseStudy2_Allmodels_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

# with open(File_name, 'wb') as f:
#     pickle.dump(AllCaseStudies, f)

# %%



