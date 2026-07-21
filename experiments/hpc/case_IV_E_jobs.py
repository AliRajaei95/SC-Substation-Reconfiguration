
import os
import subprocess



sw_list=[0,1,2]
Prob_list=[0,0.01,0.1,1]
dem=1.0
lm=1
lc=19


# Directory to save the generated scripts
output_dir = "jobs"
os.makedirs(output_dir, exist_ok=True)

for sw in sw_list:
    for prob in Prob_list:
        code_filename = f"{output_dir}/case2.3_{sw}sw_{prob}prob_14cong.py"
        with open(code_filename, "w") as code_file:
            code_file.write(f"""import pickle

from sc_substation_reconfiguration.original_MIP_model import *
from sc_substation_reconfiguration.hmmp import *
                                                    
data_file = './IEEE_14_bus_Data_CongestionEq.xlsx'

# sys=14
lm=1
dem=1.0
SolverTime=72000



print('sw is:',{sw})
print('busbar prob is:',{prob})

data=read_data_AC(File=data_file,DemFactor={dem},LineLimit=lm,busbar_prob={prob},print_data=False)

coupler_cont=data['coupler_cont']; busbar_cont = data['busbar_cont']
line_cont = data['line_cont_notradial']
lc = len(line_cont)


# #Proposed

res = SC_SR_HMMP(data,line_cont_list=line_cont,Max_Sw_bus={sw},
                                Probabilistic=True,
                                Max_iter=10,Max_FSP_iter=0,FSP_criteria=0,OSP_criteria=1)


File_name = f'Results/Case2.3/CaseStudy2.3_{prob}Prob_Proposed_14busCongEq_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

with open(File_name, 'wb') as f:
    pickle.dump(res, f)




#CO

res = SC_SR_OrgMIP(data=data,
                        cont_list=line_cont +busbar_cont+coupler_cont,
                        Zfixdict=None,
                        Probabilistic=True,
                        SolverTime=SolverTime,
                        Max_Sw_bus={sw},print_result=True)

File_name = f'Results/Case2.3/CaseStudy2.3_{prob}Prob_CO_14busCongEq_{sw}sw_{dem}dem_{lm}line_{lc}lc.pkl'

with open(File_name, 'wb') as f:
    pickle.dump(res, f)
""")
        # Generate the job script (job-i.sh)
        job_filename = f"{output_dir}/{sw}sw_{prob}p_c2.3.sh"
        with open(job_filename, "w") as job_file:
            job_file.write(f"""#!/bin/bash
                       
#SBATCH --job-name="{sw}sw_{prob}p_c2.3"
#SBATCH --output="{output_dir}/{sw}sw_{prob}p_c2.3.out"
#SBATCH --time=21:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=2GB
#SBATCH --account=research-eemcs-ese

module load py-numpy
module load py-matplotlib
module load py-tqdm


srun python {code_filename} > {output_dir}/{sw}sw_{prob}p_c2.3.log
""")
    
        # Submit the job using sbatch
        subprocess.run(["sbatch", job_filename])
