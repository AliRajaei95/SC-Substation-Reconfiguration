# Input data

This directory contains the Excel workbooks used for the paper's IEEE 14-bus,
IEEE 118-bus, and PEGASE 1354-bus case studies:

- `IEEE_14_bus_Data_PGLib_ACOPF.xlsx`
- `IEEE_118_bus_Data_PGLib_ACOPF.xlsx`
- `IEEE_1354_bus_Data_PGLib_ACOPF.xlsx`

The reader in `sc_substation_reconfiguration/original_MIP_model.py` expects the
`Bus`, `Branch`, `DemandSet`, `D2B`, `Gen`, and `G2B` sheets. The included
workbooks also provide the `LoadProfile` sheet used by the day-ahead case
study. See `read_data_AC` for the required columns and indices.

The benchmark networks are derived from the
[PGLib-OPF project](https://github.com/power-grid-lib/pglib-opf). Users of the
data should also cite PGLib-OPF and observe its applicable license and citation
requirements.
