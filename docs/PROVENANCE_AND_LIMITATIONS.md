# Numerical provenance and important boundaries

The files originated from the authors' research record, including original fixed load/PV matrices, selected historical schedules, the 69-bus revised Stage-2 solution and various published/revision-stage numerical summaries. A separate repository packaging step did not execute or amend optimizer decisions.

## Final 69-bus headline calculation

- Six PV sites: 15, 18, 27, 35, 46, 60. Final inverter nameplate rating assumption: 1.00 MW/site; minimum PF 0.90.
- Five capacitor locations: 7, 23, 30, 39, 55; step = 0.060 MVAr.
- Time series: 96 x 15 minutes; Stage-2 PV-active steps 25–39 and 61–72 (27 intervals).
- Stage-1 discrete open-line indices and capacitor-step settings were retained from the original schedule. Only inverter Volt/VAr parameters were newly tuned using GWO, warm-starting with previous values.
- Final `J_V=1.54885382313137`, 51 residual voltage violations; daily active losses `1797.509191829 kWh` and reactive losses `925.767051073 kVArh` (records in `results/final_stage2_69bus/`).

## Historical 69-bus control and model discrepancy

The historical `hybrid_dnr_cap_droop_69bus_legacy_profiles_safe.m` used an inconsistent 0.20-MW inverter-rating constant with profiles reaching approximately 1 MW/site. The original selected control CSV and the V4 schedule are identical. Replaying this historical schedule under 1.00-MW/site assumptions yields `J_V=1.58073212457191`, 51 violations. The final optimized Stage-2 program is separate and must not be described as the original full two-stage seed-42 optimization.

## Five-seed results

The seed-study algorithm includes infeasible-topology repair/resampling not present in the historical selected-schedule GWO. Therefore the five alternative runs are not fully controlled optimizer reruns with identical candidate mapping. Individual seeds, outcomes and convergence logs are kept for transparency, not reinterpreted as strict same-algorithm reproducibility. See `results/historical/results_test02_multirun_gwo/`.

## Operational validation and scope

The 96-step radiality verification checks only final connectivity and `nb-1` energized lines. It does not demonstrate safe intermediate switching operations. The saved losses come from MATPOWER AC branch-flow replays, not from the graph-only checker. Thermal ratings may be unpopulated in the benchmark and thus meaningful thermal-limit compliance is not established. 51 residual undervoltage events remain in the final 69-bus schedule. Communications, switching transients, measurement errors and closed-loop device dynamics are not reproduced here.

## How to interpret the included code

`src/v4_replay/`: primarily *replay and post-processing* of already selected schedules; does not independently create the final historical controls in one call. `src/historical_comparative/`: historical variants and comparisons with differing algorithms, assumptions, or expected output folder structure. `src/final_stage2_69bus/`: the new, inverter-only Stage-2 code and read-only replays. Results were produced using the authors' computing environment; the repository builder has checked file consistency and graph properties but has not rerun all MATPOWER cases. This is a preservation/verification package, not a universally automated experimental pipeline.
