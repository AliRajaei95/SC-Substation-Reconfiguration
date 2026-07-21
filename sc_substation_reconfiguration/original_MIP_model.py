"""Original monolithic MIP baseline from the Section IV-A case study.

This module also provides the input-data reader and shared linear AC-OPF
utilities used by the other comparison models.
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
import sys

# from BCC_Plot_Topology_v1 import plot_topology_v1


# tested with Python 3.7.0 & Gurobi 9.0

# %% [markdown]
# # Import system data 
# 

# %%
def read_data_AC(File='IEEE_14_bus_Data.xlsx',print_data=False,LineLimit=1.0,DemFactor=1.0,remove_line=[],
                 line_prob=1,coupler_prob=1,busbar_prob=0.1):
    # read the data file and convert it into a dictionary of dataframes 
     
    
    # print('Function Read data ...')
    # print('Demand factor is %.2f'%DemFactor)
    
    Sbase=100
    
    data={}
    
    
    Bus=pd.read_excel(File,sheet_name='Bus',skiprows=0,index_col=[0],usecols='A')
    Bus=list(Bus.index)                   # for b in Bus
    
    busbar=['busbar1','busbar2']     # set of busbars at each substation, # for i in busbar
    End=['Fr','To']                  #which end of line, for e in End

    


##============================ Branch bi-directional ========================

    branch=pd.read_excel(File,sheet_name='Branch',skiprows=1,index_col=[0,1,2],usecols='A:G')
    if remove_line!=[]:
        print('\n removing lines: ',remove_line)
        branch.drop(remove_line,axis=0,level=0,inplace=True)
        
    line=pd.read_excel(File,sheet_name='Branch',skiprows=1,index_col=0,usecols='A')
    if remove_line!=[]:
        line.drop(remove_line,axis=0,inplace=True)
    line=list(line.index)    #line is not bi-directional
        
    for l,i,j in branch.index:
        branch.loc[(l,j,i)]=branch.loc[(l,i,j)] #bi-directional

        


    branch['limit']=(LineLimit*branch['limit'])/Sbase


    branch['z2']= branch['r']**2+branch['x']**2
    branch['g_ij']= branch['r']/branch['z2']
    branch['b_ij']= -branch['x']/branch['z2']


    br_list=list(branch.index)
    Lines=gp.tuplelist(br_list)   #Lines is a gurobipy tupilelist that shows CONECCTIVITY of the network




    L2B=[]   #set all lines connected to bus b


    for l,i,j in Lines:
        L2B=L2B+[(l,i)]
    L2B=gp.tuplelist(L2B)     

    NumberL2B=dict.fromkeys(Bus)    #number of connected lines to each bus (substation)
    for i in NumberL2B.keys():
        n=0
        NumberL2B[i]=n
        for l in L2B.select('*',i):
            n+=1
            NumberL2B[i]=n
    



    # Demand Set
    Pdemand=pd.read_excel(File,sheet_name='DemandSet',skiprows=1,index_col=2,usecols='A:F')
    Pdemand.drop(columns=['Unnamed: 0','Unnamed: 1'],axis=1,inplace=True)


    Pdemand.rename(columns={1:'Pd'})
    Pdemand.rename(columns={2:'Qd'})
    Pdemand.fillna(0,inplace=True)
    Pdemand['Pd']=(DemFactor*Pdemand['Pd'])/Sbase          # pu
    Pdemand['Qd']=(DemFactor*Pdemand['Qd'])/Sbase          # pu

    

    DemandSet=list(Pdemand.index)                               #for d in DemandSet

    D2B_df=pd.read_excel(File,sheet_name='D2B',skiprows=0,index_col=[1,2])
    d2b_list=list(D2B_df.index)
    D2B=gp.tuplelist(d2b_list)



#==================================== Generation data =============
    # Gen_data 
    Gen_data=pd.read_excel(File,sheet_name='Gen',skiprows=2,index_col=0,usecols='A:G')
    Gen_data['Pmax']=Gen_data['Pmax']/Sbase
    Gen_data['Pmin']=Gen_data['Pmin']/Sbase
    Gen_data['Qmax']=Gen_data['Qmax']/Sbase
    Gen_data['Qmin']=Gen_data['Qmin']/Sbase
    
    # Reserve cost
    res2prod_coeff=0.5
    Gen_data['c_res']=Gen_data['b']*res2prod_coeff

    #Gen_data['b']=Gen_data['b']*Sbase
    G=list(Gen_data.index) #Set Generation

  


    # G2B is a GP object showing Gen-2-bus
    G2B_df=pd.read_excel(File,sheet_name='G2B',skiprows=0,index_col=[1,2])
    g2b_list=list(G2B_df.index)
    G2B=gp.tuplelist(g2b_list)

   
    
    line_cont=[l for l in line]
    coupler_cont=[b for b in Bus]
    busbar_cont=[str(b+'-1') for b in Bus]
    busbar_cont+=[str(b+'-2') for b in Bus]
    all_cont = line_cont + coupler_cont + busbar_cont

    Prob_cont = {0:1.0}
    # line_prob=1.0
    # coupler_prob=1.0
    # busbar_prob=0.1
    
    for c in all_cont:
        if c in line_cont:
            Prob_cont[c]=line_prob
        elif c in coupler_cont:
            Prob_cont[c]=coupler_prob
        elif c in busbar_cont:
            Prob_cont[c]=busbar_prob

    
    sub_cont_list = []
    for b in Bus:
        sub_cont_list += [[b,str(b+'-1'),str(b+'-2')]]
        
        
    # radial lines
    radial_lines=Find_radial_lines(Lines,Bus)
    line_cont_notradial=[l for l in line if l not in radial_lines]
    
    
    
    graph = Build_graph(Bus,Lines)

    
    
    if print_data==True:
        display(Pdemand)
        display(branch)
        display(Gen_data)
        
    
    data['Sbase']=Sbase
    data['Max_MIPGap']=0.01
    data['Max_timelimit']=600 #900
    data['Bus']=Bus    # for b in Bus
    data['busbar']=busbar
    data['branch']=branch
    data['Lines']=Lines
    data['line']=line
    data['L2B']=L2B
    data['NumberL2B']=NumberL2B
    data['Pdemand']=Pdemand
    data['Demandset']=DemandSet     # for d in Demand_set
    data['D2B']=D2B                  # for d,b in D2B
    data['Gen_data']=Gen_data
    data['G']=G             #for g in Gset   
    data['G2B']=G2B           
    data['vmin']=0.9
    data['vmax']=1.1
    data['beta']=0
    
    data['BigM_l']=2*branch['limit']
    data['BigM_b']=(2*3.14)/3
    data['Maxdelta']=3.4
    data['MaxV2']=2
    data['BigM_busbar']=10
    
    
    
    data['line_cont']=line_cont
    data['coupler_cont']=coupler_cont
    data['busbar_cont']=busbar_cont
    data['sub_cont_list']=sub_cont_list
    data['all_cont'] = all_cont
    data['line_cont_notradial']=line_cont_notradial
    data['radial_lines'] = radial_lines
    data['Prob_cont'] = Prob_cont
    
    data['graph'] = graph
    
    
    
    
    return data 

#====================================================================== 
def Find_radial_lines(Lines,nodes):
    
    node_list=nodes.copy()
    
    G=nx.Graph()
    
    edge_list=[]
    for l,i,j in Lines:
        edge_list=edge_list+[(i,j)]
    
    G.add_nodes_from(node_list)
    G.add_edges_from(edge_list)
    
    radial_lines=[]
    
    for l,i,j in Lines:
        G.remove_edge(i,j)
        if nx.is_connected(G)==False:
            if l not in radial_lines:
                radial_lines+=[l]
        G.add_edges_from([ (i,j)])
    
    return radial_lines



def Build_graph(nodes,Lines):

#     nodes=data['Bus']
#     Lines=data['Lines']
#     line=data['line']

    node_list=nodes.copy()

    graph=nx.Graph()

    graph.add_nodes_from(node_list)
    for l,i,j in Lines:
        graph.add_edge(i, j)

    return graph



# find n-hops lines
def Find_lines_hops(data,sub,hops=1):

    nodes=data['Bus']
    Lines=data['Lines']
    
    graph=data['graph']

        
    lines_hop_list=[]
    
    for l,i,j in Lines:
        length = (nx.shortest_path_length(graph, source=sub, target=i)
                  +nx.shortest_path_length(graph, source=sub, target=j))/2 
                

        if length - 0.5 <= hops:
            lines_hop_list+=[l]

    return lines_hop_list

# find n-hops lines
def Find_nodes_hops(data,sub,hops=1):

    Bus=data['Bus']
    Lines=data['Lines']
    
    graph=data['graph']

        
    nodes_hop_list=[]
    
    for b in Bus:
        length = nx.shortest_path_length(graph, source=sub, target=b)
                   
                

        if length <= hops:
            nodes_hop_list+=[b]

    return nodes_hop_list

# %%
def data_cb(model, where): 
    if where == gp.GRB.Callback.MIP:
        cur_obj = model.cbGet(gp.GRB.Callback.MIP_OBJBST)
        cur_bd = model.cbGet(gp.GRB.Callback.MIP_OBJBND)

        # Did objective value or best bound change?
        if model._obj != cur_obj or model._bd != cur_bd:
            model._obj = cur_obj
            model._bd = cur_bd
            model._data.append([time.time() - model._start, cur_obj, cur_bd])

# %%
def ED_market(data,print_result=False):


    #======  data

    Sbase=data['Sbase']
    Bus=data['Bus']    # for b in Bus
    branch=data['branch']
    Lines=data['Lines']
    L2B=data['L2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']     # for d in Demand_set
    D2B=data['D2B']                  # for d,b in D2B
    Gen_data=data['Gen_data']
    G=data['G']             #for g in Gset
    G2B=data['G2B']               #for g,b in G2B
    vmin = data['vmin']
    vmax = data['vmax']



    #======



#====
    
    

    
    model=gp.Model('AC_SC_OPF')
    model.Params.OutputFlag=0

    
    Pg=model.addVars(G,lb=0,vtype=GRB.CONTINUOUS,name='Pg')  
    
    

    GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost') 
    OF=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF')                          

    
    
    # equations ===============
    
    eqTotBalance=model.addConstr( (   quicksum(Pdemand.loc[d]['Pd'] for d in DemandSet) <= quicksum( Pg[g] for g in G) ),
                                 name='EqTotBalance')  

    
    model.addConstrs((Pg[g]<= Gen_data.loc[g]['Pmax']  for g in G ),name='Eq_Pgmax')
    model.addConstrs(( -Pg[g] <= -Gen_data.loc[g]['Pmin']  for g in G),name='Eq_Pgmin')



    model.addConstr(GenCost==quicksum(Gen_data.loc[g]['b']*Pg[g] for g in G ) ,name='Eq_GC')

    model.addConstr(OF==GenCost ,name='Eq_OF')

  

    model.setObjective(OF,GRB.MINIMIZE)

    model.update()


    start_time = time.time()
    model.optimize()
    # model.write('AC_SCOPF.lp')

    end_time = time.time()
    ex_time=end_time-start_time  #execution time
    

    # # =======================================================================================================================
    # # ======================================== Report Result =================================



    Pg_dict={}

    for ind,g in enumerate(G):
        Pg_dict[g]=Pg[g].x

    if print_result==True:


        print('# GenCost %.4f'%GenCost.x)
        Pgtot=0
        for g in G:
            print('Pg ',g,Pg[g].x)




    return {
        'model':model,
        'Pg':Pg_dict,
        'GenCost':OF.x,
        }


# %% [markdown]
# # market AC-SC-OPF

# %%
def AC_SC_OPF_lp(data,cont_list=[],
                    print_result=False):
    
    """"""


    #======  data

    Sbase=data['Sbase']
    Bus=data['Bus']    # for b in Bus
    branch=data['branch']
    Lines=data['Lines']
    L2B=data['L2B']
    Pdemand=data['Pdemand']
    DemandSet=data['Demandset']     # for d in Demand_set
    D2B=data['D2B']                  # for d,b in D2B
    Gen_data=data['Gen_data']
    G=data['G']             #for g in Gset
    G2B=data['G2B']               #for g,b in G2B
    vmin = data['vmin']
    vmax = data['vmax']



    #======



#====
    
    print('market running.')
    cont_list=cont_list.copy()
    cont_list+=[0]
    

    
    model=gp.Model('AC_SC_OPF')
    model.Params.OutputFlag=0

    
    Pg=model.addVars(G,lb=0,vtype=GRB.CONTINUOUS,name='Pg')  
    Qg=model.addVars(G,lb=-GRB.INFINITY ,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qg')  
    
    
    V2=model.addVars(Bus,cont_list,lb=0,vtype=GRB.CONTINUOUS,name='V2')  
    theta=model.addVars(Bus,cont_list,lb=-GRB.INFINITY ,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='theta')  
    

    

    Pflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Pflow')   
    Qflow=model.addVars(Lines,cont_list,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='Qflow')
    

    GenCost=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='GenCost') 
    OF=model.addVar(lb=0,vtype=GRB.CONTINUOUS,name='OF')                          

    
    

                


    
    # equations ===============
    
    model.addConstrs( (  
        quicksum( Pg[g] for g,i in G2B.select('*',b) )
            -  quicksum(Pdemand.loc[d]['Pd'] for d,i in D2B.select('*',b)) == 
            quicksum(Pflow[l,i,j,c] for l,i,j in Lines.select('*',b,'*')  if l!=c  )  
        for b in Bus
        for c in cont_list),name='eq_Pbalance') 
    
    model.addConstrs( (  
        quicksum( Qg[g] for g,i in G2B.select('*',b) )
            -  quicksum(Pdemand.loc[d]['Qd'] for d,i in D2B.select('*',b)) == 
            quicksum(Qflow[l,i,j,c] for l,i,j in Lines.select('*',b,'*')  if l!=c  )  
        for b in Bus
        for c in cont_list),name='eq_Qbalance') 


    model.addConstrs( ( Pflow[l,i,j,c]==0
                            for l,i,j in Lines for c in cont_list if l==c   ) , name='eqp_c')
    
    model.addConstrs( ( Qflow[l,i,j,c]==0
                        for l,i,j in Lines for c in cont_list if l==c   ) , name='eqq_c')
    


    model.addConstrs( (  Pflow[l,i,j,c] == 
    0.5*branch.loc[(l,i,j)]['g_ij']*( V2[i,c] - V2[j,c] )
    -branch.loc[(l,i,j)]['b_ij']*(theta[i,c]-theta[j,c]) #+ Ploss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c   ) , name='eqPij')
    
    model.addConstrs( (  Qflow[l,i,j,c] == 
    -0.5*branch.loc[(l,i,j)]['b_ij']*( V2[i,c]-V2[j,c])
    -branch.loc[(l,i,j)]['g_ij']*(theta[i,c]-theta[j,c]) 
     #+ Qloss[l,i,j,c]
                            for l,i,j in Lines
                           for c in cont_list if l!=c    ) , name='eqQij')

    
    


    model.addConstrs( ( V2[b,c] <= vmax**2 for b in Bus for c in cont_list) , name='eqvmax')
    model.addConstrs( ( V2[b,c] >= vmin**2  for b in Bus for c in cont_list) , name='eqvmin')
    
    
    beta=0
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij1')
    
    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] <= (1-beta)*branch.loc[(l,i,j)]['limit']
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij2')
    
    model.addConstrs( (  Pflow[l,i,j,c] + np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit']  
                            for l,i,j in Lines for c in cont_list    ) , name='eqPij3')
     
    model.addConstrs( (  Pflow[l,i,j,c] - np.tan(np.pi/6)*Qflow[l,i,j,c] >= (beta-1)*branch.loc[(l,i,j)]['limit'] 
                            for l,i,j in Lines for c in cont_list   ) , name='eqPij4')
   
    
   

    model.addConstrs((Pg[g]<= Gen_data.loc[g]['Pmax']  for g in G ),name='Eq_Pgmax')
    model.addConstrs(( -Pg[g] <= -Gen_data.loc[g]['Pmin']  for g in G),name='Eq_Pgmin')
    model.addConstrs((Qg[g]<= Gen_data.loc[g]['Qmax']  for g in G ),name='Eq_Qgmax')
    model.addConstrs(( -Qg[g] <= -Gen_data.loc[g]['Qmin']  for g in G),name='Eq_Qgmin')



    model.addConstrs((theta[Bus[0],c]==0       for c in cont_list), name='ref_bus_angle' ) 


    model.addConstr(GenCost==quicksum(Gen_data.loc[g]['b']*Pg[g] for g in G ) ,name='Eq_GC')

    model.addConstr(OF==GenCost ,name='Eq_OF')

  

    model.setObjective(OF,GRB.MINIMIZE)

    model.update()

    

    start_time = time.time()
    model.optimize()
    # model.write('AC_SCOPF.lp')

    end_time = time.time()
    ex_time=end_time-start_time  #execution time
    
    
    
    status = model.Status

    
    if status == GRB.INFEASIBLE:
        print('\n\nAC SCOPF Optimization was stopped with infeasibility!')
        print('cont list: ',cont_list)
        
        

        # Relax the bounds and try to make the model feasible
        print('\n\nThe model is infeasible; relaxing the bounds\n\n')
        orignumvars = model.NumVars
        # relaxing only variable bounds
#         model.feasRelaxS(0, False, True, False)
        # for relaxing variable bounds and constraint bounds use
#         model.feasRelaxS(0, False, True, True)
        model.feasRelaxS(0, False, False, True) #relaxing constraints

        model.optimize()

        status = model.Status
        if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
                print('The relaxed model cannot be solved \
                       because it is infeasible or unbounded')
        #sys.exit(1)
        if status != GRB.OPTIMAL:
            print('Optimization was stopped with status %d' % status)
            sys.exit(1)

        # print the values of the artificial variables of the relaxation
        print('\nSlack values:')
        slacks = model.getVars()[orignumvars:]
        for sv in slacks:
            if sv.X > 1e-9:
                print('%s = %g' % (sv.VarName, sv.X))


    # # =======================================================================================================================
    # # ======================================== Report Result =================================



    Pg_dict={}; Qg_dict={}

    for ind,g in enumerate(G):
        Pg_dict[g]=Pg[g].x
        Qg_dict[g]=Qg[g].x

        

    if print_result==True:

    


        print('# GenCost %.4f'%GenCost.x)
        Pgtot=0
        Qgtot=0
        for g in G:
            print('Pg ',g,Pg[g].x)
            Pgtot+=Pg[g].x
            Qgtot+=Qg[g].x
#             print('Qg ',g,Qg[g].x)
        # print(Pgtot,Qgtot)
        



    return {
        'model':model,
        'Pg':Pg_dict,
        'Qg':Qg_dict,
        'GenCost':OF.x,
        #'OF':OF.x,
        # 'time':ex_time,
        # 'Pgt_array':Pgt_array
        # #'Pflowdict':Pflowdict,
        # #'LMPdict':LMPdict
        }

# %%


# %% [markdown]
# # full BCC model

# %%
def BCC_v60_AC_full(data,cont_list=None,
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
    cont_list+=[0]             #full_cont_list has all the states including the 0 state of normal
    
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


    print("Number of variables:", model.NumVars)
    print("Number of constraints:", model.NumConstrs)

    NumBinVars = model.NumBinVars
    NumConVars = model.NumVars - model.NumBinVars
    NumConstrs = model.NumConstrs
    
    
    
    start_time = time.time()

    model._obj = None
    model._bd = None
    model._data = []
    model._start = time.time()
    model.optimize(callback=data_cb)
#     model.optimize()
    end_time = time.time()
    ex_time=end_time-start_time  #execution time
    
    time_data = model._data
    
    
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

    
    
    
    
    
    
    ShedCost_df = pd.DataFrame(columns=['ShedCost(c)'],index=cont_list)
    for c in cont_list:
        ShedCost_df.loc[c]=Shed_c[c].x

    ShedCost_df = pd.DataFrame(columns=['ShedCost(c)'],index=cont_list, dtype=float)
    PdShed_dic = pd.DataFrame(columns=['shed'],index=pd.MultiIndex.from_product([DemandSet,busbar,cont_list]), dtype=float)
    for c in cont_list:
        ShedCost_df.loc[c]=Shed_c[c].x
        for d in DemandSet:
            for i in busbar:
                PdShed_dic.loc[d,i,c]=Pdi_Shed[d,i,c].x
        
    
    
    if print_result==True:
        
        ex_time_df=pd.DataFrame(data=ex_time,index=['Time'],columns=['Time'])
        print('Execution Time', ex_time_df)

        ind=[]
        for l,i,j in Lines:
            for c in cont_list:
                ind=ind+[(l,i,j,c)]
            
        Pflow_df = pd.DataFrame(columns=['flow'], index=pd.MultiIndex.from_tuples(ind)) 
        for l,i,j in Lines:
            for c in cont_list:
                Pflow_df.loc[l,i,j,c]=Pflow[l,i,j,c].x
                
        
        OF_df = pd.DataFrame(columns=['OF'],data=[OF.x])
        GenCost_df = pd.DataFrame(columns=['GenCost'],data=[GenCost.x])
        RDCost_df= pd.DataFrame(columns=['RD_Cost'],data=[RDCost.x])
        TotalShedCost_df= pd.DataFrame(columns=['TotalShedCost'],data=[TotalShedCost.x])
#         ShedCost_df = pd.DataFrame(columns=['ShedCost(c)'],index=cont_list)
#         for c in cont_list:
#             ShedCost_df.loc[c]=ShedCost[c].x
        
        

                
        nb=0
        for b in Bus:    
                if z_bus[b].x==0:
                    print('bus ',b,' splitted.')
                    nb=nb+1            
        print(nb,'bus open!') 
  
        
        
        print('\n\n ======================= Results ======================= \n ')
        print('# Objective: %.3f' %(OF.x) )
        print('# Gen Cost: %.3f'%(GenCost.x) )
        print('# RD Cost : %.3f'% (RDCost.x))
        print('# Total Shedding Cost: %.3f'%(TotalShedCost.x))
        if TotalShedCost.x!=0:
            print('Average Shedding over %i contingency : %.3f MW' %( len(cont_list)-1,  100*(ShedCost_df.values.sum())/(len(cont_list)-1)  ) )
            # print('Average: ',ShedCost_df.values/ShedCost_df.value_counts)
#         display('# Shedding cost',ShedCost_df)
        

            
        # print('\nLoad shedding')
        # for c in cont_list:
        #     for d in DemandSet:
        #         for i in busbar:
        #             if Pdi_Shed[d,i,c].x> 1e-6 :
        #                 print('contingency %s , load %s-%s shed : %.3f'%(c,d,i,Pdi_Shed[d,i,c].x))
        
    
                    
                
    # ====================================================================================
    
    # Build the Topology dictionary
    
    z_bus_dict={}; z_g_dict={}; z_d_dict={}; z_li_dict={}; z_line_dict={}
    
    
    for b in Bus:
        z_bus_dict[b]=z_bus[b].x
    for g in G:
        z_g_dict[g]=z_g[g].x
    for d in DemandSet:
        z_d_dict[d]=z_d[d].x
    for l,i,j in Lines:
        z_li_dict[(l,i,j)]=z_li[l,i,j].x
    
        
    TopologyDict={'bus':z_bus_dict,
          'g':z_g_dict,
          'd':z_d_dict,
          'l_i':z_li_dict,
         }
    
 
            
    # Fix pre-contingency dispatch
    Pgi_dict={}; Qgi_dict={}
    for g in G:
        for i in busbar:
            Pgi_dict[(g,i)]=Pgi[g,i].x
            Qgi_dict[(g,i)]=Qgi[g,i].x
            
    Pg_dict = {}; Qg_dict = {}
    for g in G:
        Pg_dict[g] = Pg[g].x
        Qg_dict[g] = Qg[g].x
        print(g,Pg[g].x)
    
    
    
    

#     vol_s = {'s_vmax':{},'s_vmin':{} }
    
#     for b in Bus:
#         for c in cont_list:
#             vol_s['s_vmax'][b,c] = s_vmax[b,c].x
#             vol_s['s_vmin'][b,c] = s_vmin[b,c].x
    
    
   
    
    return {
        'time':ex_time,
        'time_iteration':time_data,
        'TopologyDict':TopologyDict,
        'Pg':Pg_dict,
        'Qg':Qg_dict,
        'Pgi':Pgi_dict,
        'Qgi':Qgi_dict,
        'MIPGap':model.MIPGap*100,
        'Cost_tot':OF.x,
        'GenCost':GenCost.x,
        'TotalShedCost':TotalShedCost.x,
        'ShedCost_df':ShedCost_df,
        'PdShed_dic':PdShed_dic,
        'NumBinVars':NumBinVars,
        'NumConVars':NumConVars,
        'NumConstrs':NumConstrs 
#         'vol_s':vol_s
        }
    
    
    
    

# %%


# %%
# NumBinVars = model.NumBinVars
# NumConVars = model.NumVars - model.NumBinVars
# NumConstrs = model.NumConstrs


# 'NumBinVars':NumBinVars,
# 'NumConVars':NumConVars,
# 'NumConstrs':NumConstrs

# %% [markdown]
# # full model

# %%
# data=read_data_AC(File='IEEE_14_bus_Data_PGLib_ACOPF.xlsx',print_data=False)

# line_cont = data['line_cont_notradial']
# coupler_cont =data['coupler_cont']
# busbar_cont = data['busbar_cont']

# # Pg_market=AC_SC_OPF_lp(data=data,cont_list=line_cont,print_result=True)['Pg']



# res_full=BCC_v60_AC_full(data=data,
#                         cont_list=[] #line_cont+coupler_cont+busbar_cont #busbar_cont+coupler_cont+line_cont
#                         # ,FixedCost=True,Pg_market=Pg_market,Alpha=0.1
#                         # ,Probabilistic=True
#                         ,Up_redispatch=0,Dn_redispatch=1.0,Zfixdict=None,SolverTime=60,Threads=2,
#                                                    Max_Sw_bus=1,print_result=True)
# # # # # # for g in data['G']:
# # # # # #     for i in data['busbar']:
# # # # # #         print(g,i,res_full['Pg'][g,i])
# # # # # #         print(g,i,res_full['Qg'][g,i])


# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%



