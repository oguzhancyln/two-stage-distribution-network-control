# Data dictionary

| Location | Variable / fields | Shape or meaning |
|---|---|---|
| `data/profiles/myloadspvs.mat` | `myloads`, `myPVs` | 33 x 96 load multipliers and 5 x 96 PV traces (inspect source for scaling assumptions) |
| `data/profiles/myloads.mat` | `myloads` | 69 x 96 load multipliers |
| `data/profiles/mypvs.mat` | `PV` | 6 x 96 PV power outputs (MW) |
| `data/networks/case69modified.m` | MATPOWER `mpc` | 69 buses; 73 branches; branch row number is the switch index |
| `data/networks/branch_endpoints_case33bw.csv` | `index,from_bus,to_bus` | Graph-only copy of the 37 standard 33-bus case branch endpoints; not a power-flow case |
| `data/historical_schedules/control_schedule_33bus.csv` | `TimeStep,Mode,OpenLines,CapSteps,...` | Original selected 33-bus 96-step controls |
| `data/historical_schedules/control_schedule_69bus.csv` | same | Original selected 69-bus controls with historical Stage-2 values |
| `results/final_stage2_69bus/control_schedule_69bus_STAGE2_REOPT.csv` | same | New Stage-2 parameters; identical Stage-1 open-line and CAP choices |
| `results/final_stage2_69bus/comparison_69bus.csv` | `Method,JV,...` | Historical corrected-rating replay vs final Stage-2 results |
| `results/final_stage2_69bus/audit_69_full_replay_summary.csv` | per-method metrics | 96-step voltage and integrated P/Q losses, switch counts and Q statistics |
| `results/final_stage2_69bus/audit_69_new_stage2_comparisons.csv` | per-comparator metrics | Fixed VVC, capacitor sensitivity, acceptance gate, alternate capability |

`OpenLines` and `CapSteps` vectors in control CSVs use bracketed numerical sequences. Bus and branch identifiers follow MATPOWER one-based conventions. Profiles, control programs, and raw result files were copied without numerical alteration.
