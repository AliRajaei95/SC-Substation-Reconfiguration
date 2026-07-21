"""Iterative topology-fixing heuristic baseline from Section IV-A.

The method repeatedly solves the MIP while fixing and refining candidate
substation configurations under the selected contingency set.
"""

# %% [markdown]
# # Busbar Coupler Contingency (BCC) Project
#  
# # BCC_full_model + ACOPF equations (linear)
# 
# 
# 
# 

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
from tqdm import tqdm
import sys

from .original_MIP_model import *  # Shared data preparation and AC-OPF utilities.

# # Iterative_Heuristic
# %%

def Create_v65_AC_MIP(data,cont_list=None,
                                        Max_Sw_bus=0,
                                    Up_redispatch=0,Dn_redispatch=1.0,
                                           Zfixdict=None,#fixed topology
                                    Zinitial=None,
                                   PQgFix=None, PgFix=None,
                                   FixedCost=False,Alpha=0,Pg_market=None,
                                   Probabilistic=False,
                                    SolverTime=600,Threads=None,
                                    line_shedding=True,
                                     Vol0_li=None
                                          ):
    
    
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
    
    model=gp.Model('SBF_Preventive SC line contingency')
    model.Params.OutputFlag=1
    
    cont_list=cont_list.copy()
     #full_cont_list has all the states including the 0 state of normal
    
#     print("======== contingency list ======")
#     print('        ',cont_list,'\n\n')
    
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
                Vol0_li[l,i,j,c]=1
    
    
  
    
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
    
#     Ploss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Ploss') 
#     epsilon=model.addVars(Lines,cont_list,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='epsilon') 
#     Qloss=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qloss')

    
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


#### ========== initial value


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
            
            
        
        
    #Fix the dispatch
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
        
        
         

# Binary Equations ==================================================

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
    
    # max Switching ================
    eq_MaxSw_bus=model.addConstr((  quicksum( (1-z_bus[b]) for b in Bus) <= Max_Sw_bus) , name='eq_MaxSw_bus') 
    

#     symmetry 
    if Zfixdict==None:
        for b in Bus:
            lmin,b=L2B.select('*',b)[0]
            lmin,i,j=Lines.select(lmin,b,'*')[0]
            eq_symmetry1=model.addConstr((  z_li[lmin,i,j] ==0  ), name='eq_symmetry')
#         print(lmin,i,j)

    eq_zbus_ref=model.addConstr( (  z_bus[Bus[0]]==1    ), name='eq_zbus_ref')

        
    
    #== Gen
    
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
    
    #Qg
    eq_Qg1min=model.addConstrs( (   (1-z_g[g])*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c]   
                             for g in G for c in cont_list  ), name='eq_Qg1min')
    eq_Qg1max=model.addConstrs( (    Qgi[g,'busbar1']+dQgi_up[g,'busbar1',c]-dQgi_dn[g,'busbar1',c] <=(1-z_g[g])*Gen_data.loc[g]['Qmax']  
                             for g in G for c in cont_list ), name='eq_Qg1max')
    
    eq_Qg2min=model.addConstrs( (   z_g[g]*Gen_data.loc[g]['Qmin'] <= Qgi[g,'busbar2'] +dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]  
                             for g in G for c in cont_list  ), name='eq_Qg2min')
    eq_Qg2max=model.addConstrs( (   Qgi[g,'busbar2']+dQgi_up[g,'busbar2',c]-dQgi_dn[g,'busbar2',c]  <= z_g[g]*Gen_data.loc[g]['Qmax']   
                             for g in G for c in cont_list  ), name='eq_Pg2max')
    

    
    
    
    # =======   Corrective redispatch/ reserve equations =====================================
    
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

    

    #Pd
    eq_Pd1=model.addConstrs( ( Pdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Pd']   
                            for d in DemandSet ), name='eq_Pd1')
    eq_Pd2=model.addConstrs( ( Pdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Pd']   
                            for d in DemandSet ), name='eq_Pd2')
    eq_Qd1=model.addConstrs( ( Qdi[d,'busbar1'] == (1-z_d[d])*Pdemand.loc[d]['Qd']   
                            for d in DemandSet ), name='eq_Qd1')
    eq_Qd2=model.addConstrs( ( Qdi[d,'busbar2'] == z_d[d]*Pdemand.loc[d]['Qd']   
                            for d in DemandSet ), name='eq_Qd2')
    
    # Load sheddding
    model.addConstrs( (  Pdi_Shed[d,'busbar1',c]==Pdi[d,'busbar1'] 
                              for b in Bus for d,b in D2B.select('*',b) 
                               for c in cont_list if c==str(b+'-1') ) , name='eq_Pd1Shed')
    model.addConstrs( (  Pdi_Shed[d,'busbar2',c]==Pdi[d,'busbar2'] 
                              for b in Bus for d,b in D2B.select('*',b) 
                               for c in cont_list if c==str(b+'-2') ) , name='eq_Pd2Shed')
    model.addConstrs( (  Pdi_Shed[d,i,c]<=Pdi[d,i] 
                              for d in DemandSet for i in busbar 
                               for c in cont_list) , name='eq_PdShedlimit')
    # model.addConstrs( (  Pdi_Shed[d,i,c]==0 
    #                           for d in DemandSet for i in busbar 
    #                            for c in cont_list
    #                   if c==0) , name='eq_PdShedlimit0')
    
    model.addConstrs( (  Qdi_Shed[d,i,c]== (Pdemand.loc[d]['Qd'])/(Pdemand.loc[d]['Pd'])*Pdi_Shed[d,i,c]
                              for d in DemandSet for i in busbar 
                               for c in cont_list) , name='eq_QdShed')
    
    # if line_shedding==False:
    #     model.addConstrs( (  Pdi_Shed[d,i,c]==0 
    #                               for d in DemandSet for i in busbar 
    #                                for c in cont_list
    #                           if c in data['radial_lines'] ) , name='eq_PdShed_line0')


###=== PF equations ===================================================

##===== line Contingency 
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

    ##========== flow limits   
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
    
    model.addConstrs( (  Qflow[l,i,j,c] == 
    -0.5*branch.loc[(l,i,j)]['b_ij']*(V2_li[l,i,j,c] - V2_li[l,j,i,c])
    -branch.loc[(l,i,j)]['g_ij']*(delta_li[l,i,j,c]-delta_li[l,j,i,c]) 
#     -V2_li[l,i,j,c]*(branch.loc[(l,i,j)]['b']/2) #+ Qloss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')
    
    #Loss
#     model.addConstrs( (  Ploss[l,i,j,c] == 
#                 branch.loc[(l,i,j)]['g_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])\
#                 -0.5*branch.loc[(l,i,j)]['g_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2  
#                        +epsilon[l,i,j,c]
#                             for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqPloss')
#     model.addConstrs( (  Qloss[l,i,j,c] == 
#                 -branch.loc[(l,i,j)]['b_ij']*((Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])/(Vol0_li[l,i,j,c]+Vol0_li[l,j,i,c]))*(V2_li[l,i,j,c]-V2_li[l,j,i,c])\
#                 +0.5*branch.loc[(l,i,j)]['b_ij']*(Vol0_li[l,i,j,c]-Vol0_li[l,j,i,c])**2  
#                             for l,i,j in Lines for c in cont_list if l!=c  ) , name='eqQloss')


    # voltage
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
    
    #V2
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


    #busbar & coupler contingency
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
    

## ==== Balance ================ 
    
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



    # cost!
    
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
#     model.write('BCC_v60_AC_full.lp')


    
    
    return model
    
    

# %%
def Solve_MIP_v65(data,TopologyFix,model0,cont_list,K_Opt=1,print_result=False):
    
    """Solve for K_opt

    fix the ones that need to be fixed and count them. 
    the rest can be changed with max k
    
    
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
    
    model=model0.copy()
    model.Params.OutputFlag=1


#=============================== Equations ==================================================

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


    # long cut
    model.addConstr((  quicksum( (1- model.getVarByName('z_bus['+str(b)+']') ) for b in Bus)
                     + quicksum( model.getVarByName('z_li['+str(l)+','+str(i)+','+str(j)+']')   for l,i,j in Lines  ) 
                      + quicksum( model.getVarByName('z_g['+str(g)+']') for g in G)
                       +quicksum( model.getVarByName('z_d['+str(d)+']') for d in DemandSet)
                         <= counter + K_Opt) , name='eq_MaxSw_bus') 



    


    model.update()
    #model.write('BCC_Benders_v4_OSP.lp')
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time 


    
    
    #===================== infeasibility debug
    
    status = model.Status

    
    if status == GRB.INFEASIBLE:
        print('\n\nOptimization was stopped with infeasibility!')
        

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
    
    
 

    # =======================================================================================================================
    # =======================================================================================================================
    # ======================================== Report Result ================================= 

    
    
    
    
    
    
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
        # print('# Total Shedding Cost: %.3f'%(TotalShedCost.x))
        # if TotalShedCost.x!=0:
        #     print('Average Shedding over %i contingency : %.3f MW' %( len(cont_list)-1,  100*(ShedCost_df.values.sum())/(len(cont_list)-1)  ) )
            # print('Average: ',ShedCost_df.values/ShedCost_df.value_counts)
#         display('# Shedding cost',ShedCost_df)
        

                    
                
    # ====================================================================================
    
    # Build the Topology dictionary
    
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

    # print(TopologyFix_updated)



    
 
            
    # Fix pre-contingency dispatch
            
    Pg_dict = {}; Qg_dict = {}
    for g in G:
        Pg_dict[g] = model.getVarByName('Pg['+str(g)+']').x
        Qg_dict[g] =model.getVarByName('Qg['+str(g)+']').x
        # print(g,Pg_dict[g])
    
    
    
    

    
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
        'PdShed_dic':PdShed_dic
        }
    
    

# %% main iterative heuristic function

def BCC_v65_AC_iterative_heuristic(data,cont_list=None,
                                   K_Opt=1, #sw change in each iteration
                                   Max_iteration=10000,
                                        Max_Sw_bus=0,
                                    Up_redispatch=0,Dn_redispatch=1.0,
                                           Zfixdict=None,#fixed topology
                                    Zinitial=None,
                                   PQgFix=None, PgFix=None,
                                   FixedCost=False,Alpha=0,Pg_market=None,
                                   Probabilistic=False,
                                    SolverTime=600,Threads=None,
                                    line_shedding=True,
                                     Vol0_li=None,
                                           print_result=False):
    

    cont_list = cont_list + [0]

    model0 = Create_v65_AC_MIP(data=data,cont_list=cont_list,
                                        Max_Sw_bus=Max_Sw_bus,
                                    Up_redispatch=Up_redispatch,Dn_redispatch=Dn_redispatch,
                                           Zfixdict=Zfixdict,#fixed topology
                                    Zinitial=Zinitial,
                                   PQgFix=PQgFix, PgFix=PgFix,
                                   FixedCost=FixedCost,Alpha=Alpha,Pg_market=Pg_market,
                                   Probabilistic=Probabilistic,
                                    SolverTime=SolverTime,Threads=Threads,
                                    line_shedding=line_shedding,
                                     Vol0_li=Vol0_li )
    



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

    for iter in tqdm(range(Max_iteration)):


        res = Solve_MIP_v65(data=data,TopologyFix=TopologyFix_updated,model0=model0,cont_list=cont_list,K_Opt=K_Opt,print_result=print_result)


        TopologyFix_updated = res['TopologyFix_updated']
        UB = res['Cost_tot']
        tot_time += res['time']

        Topology_k[iter] = res['TopologyDict']
        Topology_updated_k[iter] = res['TopologyFix_updated']
        UB_k[iter] = UB
        time_iteration[iter] = tot_time
        GenCost_k[iter] = res['GenCost']


        if UB < UB_min*0.999:
            UB_min = UB
            print('iter %i, UB: %.2f'%(iter,UB))
        else:
            print('Terminate! ===> UB is %.2f'%UB)
            break

        if tot_time > 2*24*3600:
            print('time is up at %i, terminate!'%tot_time)



    




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
        'time':tot_time,
    }
