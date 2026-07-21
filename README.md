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
paper's full formulation, heuristic multi-master problem (HMMP), optimality-cut
variant, Benders baseline, and sequential/iterative heuristics.

<p align="center">
  <img src="figures/substation_topology.jpg" width="500" alt="Substation topology">
</p>

## Repository layout

| Path | Contents |
| --- | --- |
| `sc_substation_reconfiguration/original_MIP_model.py` | Original monolithic MIP baseline, data reader, market dispatch, and AC SC-OPF |
| `sc_substation_reconfiguration/hmmp.py` | Proposed heuristic multi-master method |
| `sc_substation_reconfiguration/hmmp_optimality.py` | HMMP optimality-cut variant |
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
    AC_SC_OPF_lp,
    BCC_v60_AC_full,
    read_data_AC,
)

data = read_data_AC(
    File="data/IEEE_14_bus_Data_PGLib_ACOPF.xlsx",
    DemFactor=1.0,
    LineLimit=1.0,
    busbar_prob=0.05,
)

market = AC_SC_OPF_lp(data=data, cont_list=[], print_result=False)
result = BCC_v60_AC_full(
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

The source files retain the original function names used in the experiments so
results can be traced back to the paper. HPC scripts are archival research
drivers; read [`experiments/hpc/README.md`](experiments/hpc/README.md) before
running them because the job generators invoke `sbatch`.

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
