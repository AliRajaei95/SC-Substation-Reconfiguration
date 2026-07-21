# Input data

The case studies expect Excel workbooks derived from the IEEE 14-bus and
118-bus systems and the PEGASE 1354-bus system. The workbooks are not included
in this release. Place authorized copies in this directory (or pass their paths
to `read_data_AC`).

The reader in `sc_substation_reconfiguration/original_MIP_model.py` expects the sheets
`Bus`, `Branch`, `DemandSet`, `D2B`, `Gen`, and `G2B`. See that function for the
required columns and indices. Benchmark network data can be obtained from the
[PGLib-OPF project](https://github.com/power-grid-lib/pglib-opf); conversion to
the workbook schema is currently required.

Do not commit third-party datasets unless their licenses permit redistribution.
