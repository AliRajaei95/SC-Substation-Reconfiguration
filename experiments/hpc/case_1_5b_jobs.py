import os
import subprocess



sw_list = [0,1,2]

dem=1.0
lm=1
lc=0
prob=0.05



# Directory to save the generated scripts
output_dir = "jobs"
os.makedirs(output_dir, exist_ok=True)

for sw in sw_list:
    
    code_filename = f"{output_dir}/case1.5_PgFix_{sw}sw_14cong.py"
    with open(code_filename, "w") as code_file:
        code_file.write(f"""import pickle

from sc_substation_reconfiguration.original_MIP_model import *

                                                    
data_file = './IEEE_14_bus_Data_CongestionEq.xlsx'

# sys=14
lm=1
dem=1.0
SolverTime=36000



print('sw is:',{sw})

data=read_data_AC(File=data_file,DemFactor=dem,LineLimit=lm,busbar_prob={prob},print_data=False)

coupler_cont=data['coupler_cont']; busbar_cont = data['busbar_cont']
line_cont = [] #data['line_cont_notradial']
lc = len(line_cont)


Pg_market=AC_SC_OPF_lp(data=data,cont_list=line_cont,print_result=True)['Pg']


res = BCC_v60_AC_full(data=data,
                        cont_list=line_cont+coupler_cont+busbar_cont,
                        Probabilistic=True,
                        Zfixdict=None,
                         PgFix=Pg_market,
                           Max_Sw_bus={sw},SolverTime=SolverTime,print_result=True)
    

    
File_name = f'Results/Case1.5/CaseStudy1.5_PgFix_14busCongEq_{sw}sw_{dem}dem_{lm}line_{lc}lc_{prob}prob.pkl'

with open(File_name, 'wb') as f:
    pickle.dump(res, f)

""")
    # Generate the job script (job-i.sh)
    job_filename = f"{output_dir}/{sw}sw_c1.5.sh"
    with open(job_filename, "w") as job_file:
        job_file.write(f"""#!/bin/bash
                       
#SBATCH --job-name="{sw}sw_c1.5"
#SBATCH --output="{output_dir}/{sw}sw_c1.5.out"
#SBATCH --time=11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=2GB
#SBATCH --account=research-eemcs-ese

module load py-numpy
module load py-matplotlib
module load py-tqdm


srun python {code_filename} > {output_dir}/{sw}sw_c1.5.log
""")
    # Submit the job using sbatch
    subprocess.run(["sbatch", job_filename])
