# Security-Constrained Substation Reconfiguration

Research code accompanying:

> Ali Rajaei, Olayiwola Arowolo, and Jochen L. Cremer, “Security-Constrained
> Substation Reconfiguration Considering Busbar and Coupler Contingencies,”
> *IEEE Transactions on Power Systems*, 2026.

[Read the open-access preprint](https://arxiv.org/abs/2603.04203)

## Overview

Substation reconfiguration through busbar splitting can relieve congestion and
reduce operating cost. This implementation extends security-constrained
reconfiguration to line, busbar, and coupler contingencies. It includes the
paper's full formulation, heuristic multi-master problem (HMMP), Benders
baseline, and sequential/iterative heuristics.

## Motivation

On **January 8, 2021**, the European power system experienced a major
disturbance that split the continental grid into two areas. The event was
triggered by the **tripping of a highly loaded busbar coupler**, which led to
cascading failures across the network.

The post-event analysis showed that the **substation topology had not been
adjusted after a transmission line outage**, and the **coupler contingency had
not been included in the N-1 security analysis**.

This incident highlights the importance of explicitly considering
**substation elements such as busbars and couplers** when determining secure
grid configurations.

<p align="center">
  <img src="figures/europe_grid_split.jpg" width="650" alt="European grid split">
</p>

*European system split on January 8, 2021 (adapted from the ENTSO-E report).*

<p align="center">
  <img src="figures/substation_topology.jpg" width="500" alt="Substation topology">
</p>

*Illustration of the substation topology involved in the event.*

To address this challenge, our work proposes a **security-constrained
substation reconfiguration framework** that considers **line, coupler, and
busbar contingencies**, while remaining computationally scalable for large
power systems.

## Abstract

Substation reconfiguration via busbar splitting can mitigate transmission grid
congestion and reduce operational costs. However, existing approaches neglect
the security of substation topology, particularly for substations without
busbar splitting (i.e., closed couplers), which can lead to severe
consequences. Additionally, the computational complexity of optimizing
substation topology remains a challenge.

This paper introduces a MILP formulation for security-constrained substation
reconfiguration (SC-SR), considering N-1 line, coupler, and busbar
contingencies to ensure secure substation topology. To efficiently solve this
problem, we propose a heuristic approach with multiple master problems (HMMP).
A central master problem optimizes dispatch, while independent substation
master problems determine individual substation topologies in parallel. Linear
AC power flow equations ensure power-flow accuracy, while feasibility and
optimality subproblems evaluate contingency cases.

The proposed HMMP significantly reduces computational complexity and enables
scalability to large power systems. Case studies on the IEEE 14-bus, IEEE
118-bus, and PEGASE 1354-bus systems show the effectiveness of the approach in
mitigating the impact of coupler and busbar tripping, balancing system security
and cost, and improving computational efficiency.

## Repository status

The research implementation is available and is being prepared for its public
release. The repository currently includes:

- The proposed **SC-SR formulation** implemented in Python
- The **heuristic multi-master problem** solution approach
- Baseline approaches, including Benders decomposition variants and heuristic
  methods
- Case-study scripts for the IEEE benchmark systems

## Repository layout

| Path | Contents |
| --- | --- |
| `sc_substation_reconfiguration/original_MIP_model.py` | Original monolithic MIP baseline, data reader, market dispatch, and AC SC-OPF |
| `sc_substation_reconfiguration/hmmp.py` | Proposed heuristic multi-master method |
| `sc_substation_reconfiguration/benders.py` | Classical Benders baseline |
| `sc_substation_reconfiguration/iterative_heuristic.py` | Iterative heuristic |
| `sc_substation_reconfiguration/sequential_heuristic.py` | Sequential heuristic |
| `sc_substation_reconfiguration/fixed_1354.py` | Fixed-topology PEGASE 1354-bus method |
| `experiments/hpc/` | Slurm job generators and paper experiment scripts |
| `data/` | Input-workbook format and data-source notes |

## Installation

Python 3.10 or newer is recommended. Create an isolated environment and install
the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

The models use the commercial [Gurobi Optimizer](https://www.gurobi.com/). A
working Gurobi license is required; eligible academics can request a free
academic license. The code used Gurobi 10.0.3 for the reported experiments.

## Data and basic use

The original Excel benchmark workbooks are not redistributed here. Their
expected sheets and source information are documented in
[`data/README.md`](data/README.md). Once an input workbook is available:

```python
from sc_substation_reconfiguration.original_MIP_model import (
    solve_ac_sc_opf,
    SC_SR_OrgMIP,
    read_data_AC,
)

data = read_data_AC(
    File="data/IEEE_14_bus_Data_PGLib_ACOPF.xlsx",
    DemFactor=1.0,
    LineLimit=1.0,
    busbar_prob=0.05,
)

market = solve_ac_sc_opf(data=data, cont_list=[], print_result=False)
result = SC_SR_OrgMIP(
    data=data,
    cont_list=data["line_cont_notradial"] + data["busbar_cont"] + data["coupler_cont"],
    Pg_market=market["Pg"],
    FixedCost=True,
    Probabilistic=True,
    Max_Sw_bus=2,
    SolverTime=3600,
    print_result=True,
)
```

The public solver names follow the method labels used in the paper, including
`SC_SR_OrgMIP`, `SC_SR_Benders`, `SC_SR_HMMP`, `SC_SR_1OptH`, and
`SC_SR_SeqH`. Internal model builders and solvers include their owning method
in the function name. HPC scripts are archival research drivers; read
[`experiments/hpc/README.md`](experiments/hpc/README.md) before running them
because the job generators invoke `sbatch`.

## Reproducibility notes

- Solver results can depend on the Gurobi version, license limits, machine, and
  thread configuration.
- Randomized routines should be given an explicit seed when comparing runs.
- Generated jobs, logs, result pickles, and local PDF copies are intentionally
  ignored by Git.
- The open-access paper is linked above instead of bundling a potentially
  publisher-restricted PDF.

## Citation

GitHub can export the citation metadata in [`CITATION.cff`](CITATION.cff). A
BibTeX entry is also provided:

```bibtex
@article{rajaei2026security,
  title   = {Security-Constrained Substation Reconfiguration Considering Busbar and Coupler Contingencies},
  author  = {Rajaei, Ali and Arowolo, Olayiwola and Cremer, Jochen L.},
  journal = {IEEE Transactions on Power Systems},
  year    = {2026}
}
```

## License

The software is released under the [MIT License](LICENSE). Third-party
benchmark data and the paper are governed by their respective licenses.
