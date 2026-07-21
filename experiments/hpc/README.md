# HPC case-study scripts

These scripts preserve the experiment setup used for the paper. Run them from
the repository root so the package and input paths resolve correctly.

The `*_jobs.py` programs generate Slurm scripts and immediately submit them
with `sbatch`. Review the DelftBlue-specific account, partition, modules, time,
and memory settings before running. Generated jobs, logs, and result pickle
files are ignored by Git.

Expected benchmark workbooks are documented in `data/README.md`.
