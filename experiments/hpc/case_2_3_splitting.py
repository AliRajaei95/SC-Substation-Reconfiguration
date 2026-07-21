# %% [markdown]
# # Case studies 2
# # Proposed solution methods
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
import pickle

# from BCC_Plot_Topology_v1 import plot_topology_v1


# tested with Python 3.10.0 & Gurobi 9.5

# %%
from sc_substation_reconfiguration.original_MIP_model import *
from sc_substation_reconfiguration.hmmp import *
from sc_substation_reconfiguration.hmmp_optimality import *





# %%
sys=14
lm=1
dem=1.2
SolverTime=18000


sw_list=[3]
Prob_list=[0,0.01,0.1,1]

prob=0
data=read_data_AC(File=f'IEEE_{sys}_bus_Data_PGLib_ACOPF.xlsx',DemFactor=dem,LineLimit=lm,busbar_prob=prob,print_data=False)

coupler_cont =data['coupler_cont']; busbar_cont = data['busbar_cont']
line_cont = data['line_cont_notradial']
lc = len(line_cont)

# Pg_market=AC_SC_OPF_lp(data=data,cont_list=line_cont,print_result=True)['Pg']

# Pg_market=ED_market(data=data,print_result=True)['Pg']







# %% [markdown]
# # PgFix

# # #Proposed

# for sw in sw_list:
#     CaseStudies = {}

#     res = BCC_v61_AC_MultiBenders_Main(data,line_cont_list=line_cont,Max_Sw_bus=sw,
#                                     Pg_market_fix=Pg_market,
#                                     Max_iter=10,Max_FSP_iter=0,FSP_criteria=0,OSP_criteria=1)


#     CaseStudies['Proposed'] = res


#     File_name = f'Results/Case2.3/CaseStudy2.3_EM_Proposed_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

#     with open(File_name, 'wb') as f:
#         pickle.dump(CaseStudies, f)




#     # %%
#     CaseStudies = {}

#     res = BCC_v66_AC_MultiBenders_OptimalSplit(data,line_cont_list=line_cont,
#                                     Pg_market_fix=Pg_market,
#                                     Max_iter=10,Max_FSP_iter=0,FSP_criteria=0,OSP_criteria=1,
#                                     Max_Sw_bus=sw)


#     CaseStudies['Proposed_new'] = res


#     File_name = f'Results/Case2.3/CaseStudy2.3_EM_Proposed_new_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

#     with open(File_name, 'wb') as f:
#         pickle.dump(CaseStudies, f)



#     #CO
#     CaseStudies = {}

#     res = BCC_v60_AC_full(data=data,
#                             cont_list=line_cont +busbar_cont+coupler_cont,
#                             Zfixdict=None,
#                             PgFix=Pg_market,
#                             SolverTime=SolverTime,
#                             Max_Sw_bus=sw,print_result=True)

#     CaseStudies['CO'] = res

#     File_name = f'Results/Case2.3/CaseStudy2.3_EM_CO_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

#     with open(File_name, 'wb') as f:
#         pickle.dump(CaseStudies, f)


# #combine all
# AllCaseStudies = {sw:{} for sw in sw_list}


# Cases=['Proposed','Proposed_new','CO']

# market='EM' #PgFix

# for sw in sw_list:
#     for case in Cases:
#         File_name = f'Results/Case2.3/CaseStudy2.3_{market}_{case}_14bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'
#         with open(File_name, 'rb') as f:
#             AllCaseStudies[sw][case] = pickle.load(f)[case]

        
# File_name = f'Results/Case2.3/CaseStudy2.3_All_{market}_{dem}dem_{lm}line_{lc}lc.pkl'

# with open(File_name, 'wb') as f:
#     pickle.dump(AllCaseStudies, f)






# ===============================================================================================================================================================


# %% [markdown]
# # Probabilistic
# 

# %%
#Prob 


for sw in sw_list:
    for prob in Prob_list:

        print('sw is:',sw)
        print('busbar prob is:',prob)

        data=read_data_AC(File=f'IEEE_{sys}_bus_Data_PGLib_ACOPF.xlsx',DemFactor=dem,LineLimit=lm,busbar_prob=prob,print_data=False)


        # #Proposed
        CaseStudies = {}

        res = BCC_v61_AC_MultiBenders_Main(data,line_cont_list=line_cont,Max_Sw_bus=sw,
                                        Probabilistic=True,
                                        Max_iter=10,Max_FSP_iter=10,FSP_criteria=0,OSP_criteria=1)


        CaseStudies['Proposed'] = res


        File_name = f'Results/Case2.3/CaseStudy2.3_{prob}Prob_Proposed_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

        with open(File_name, 'wb') as f:
            pickle.dump(CaseStudies, f)




        # %%
        CaseStudies = {}

        res = BCC_v66_AC_MultiBenders_OptimalSplit(data,line_cont_list=line_cont,
                                        Probabilistic=True,
                                        Max_iter=10,Max_FSP_iter=10,FSP_criteria=0,OSP_criteria=1,
                                        Max_Sw_bus=sw)


        CaseStudies['Proposed_new'] = res


        File_name = f'Results/Case2.3/CaseStudy2.3_{prob}Prob_Proposed_new_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

        with open(File_name, 'wb') as f:
            pickle.dump(CaseStudies, f)



        #CO
        CaseStudies = {}

        res = BCC_v60_AC_full(data=data,
                                cont_list=line_cont +busbar_cont+coupler_cont,
                                Zfixdict=None,
                                Probabilistic=True,
                                SolverTime=SolverTime,
                                Max_Sw_bus=sw,print_result=True)

        CaseStudies['CO'] = res

        File_name = f'Results/Case2.3/CaseStudy2.3_{prob}Prob_CO_{sys}bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

        with open(File_name, 'wb') as f:
            pickle.dump(CaseStudies, f)




# AllCaseStudies = {(sw,prob):{} for sw in sw_list for prob in Prob_list}


# Cases=['Proposed','Proposed_new','CO']


# for sw in sw_list:
#     for prob in Prob_list:
#         for case in Cases:
#             # if case == 'CO':
#             File_name = f'Results/CaseStudy2.3_{prob}Prob_{case}_14bus_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'
#             with open(File_name, 'rb') as f:
#                 AllCaseStudies[(sw,prob)][case] = pickle.load(f)[case]
            

        
# File_name = f'Results/CaseStudy2.3_All_Prob_{dem}dem_{lm}line_{lc}lc.pkl'

# with open(File_name, 'wb') as f:
#     pickle.dump(AllCaseStudies, f)


