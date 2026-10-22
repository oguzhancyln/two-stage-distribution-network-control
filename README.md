# Two-stage voltage control in PV-rich radial distribution networks

Research code, fixed input profiles, selected control schedules, and verification records associated with the manuscript “A Two-Stage Voltage Control Framework Combining Network Reconfiguration and Smart Inverter Support in PV-Rich Distribution Networks” (IEEE Access submission, 2026).

Authors: Oğuzhan Ceylan, Alper Savaşçı, Alexis Polycarpou.

> Scope: This repository documents the archived numerical work and supports several targeted checks. It is not a one-command, end-to-end regeneration of every manuscript table and figure. The original selected schedules, later Stage-2 improvement, and alternate five-seed implementation are intentionally distinguished. Before public publication, the authors should approve the data and source-code release and choose an appropriate license.

## Layout

| Location | Contents |
|---|---|
| `data/profiles/` | Original deterministic MATLAB `.mat` load/PV arrays (unchanged) |
| `data/networks/` | Modified 69-bus MATPOWER case; indexed 33-bus branch endpoints |
| `data/historical_schedules/` | Selected original 33-bus and 69-bus 96-step decisions |
| `src/final_stage2_69bus/` | 69-bus Stage-2 1-MW/site optimizer and replay/audit scripts |
| `src/v4_replay/` | Stored-schedule V4 ablation and plotting code |
| `src/historical_comparative/` | Historical optimizer variants, comparisons, multi-seed and timing scripts |
| `src/validation/` | Independent Python radial topology / switching-count audit |
| `results/final_stage2_69bus/` | New inverter schedule, voltages, comparisons and loss audit |
| `results/historical/` | Historical summaries, ablations and five-seed output CSVs |
| `tools/` | Repository integrity check and a separate local execution workspace builder |
| `docs/` | Provenance, validation scope, data inventory and release checklist |

## Quick start: independently check radiality (Python standard library)

From the repository root:

```bash
python3 src/validation/independent_radiality_audit.py \
  --case 33 \
  --branch-file data/networks/branch_endpoints_case33bw.csv \
  --schedule data/historical_schedules/control_schedule_33bus.csv

python3 src/validation/independent_radiality_audit.py \
  --case 69 \
  --branch-file data/networks/case69modified.m \
  --schedule results/final_stage2_69bus/control_schedule_69bus_STAGE2_REOPT.csv
```

Expected: 96/96 radial endpoint checks per system. Recorded initial/subsequent switching counts: 33-bus line 8/8 and capacitor-step 13/8; 69-bus line 10/16 and capacitor-step 14/13. This is graph connectivity only, not an AC load-flow, switching-transient, or thermal-limit validation.

Check release files and their hashes:

```bash
python3 tools/verify_repository.py
```

## Optional MATPOWER/Octave replay: final 69-bus case

The original Octave scripts assume a flat `ALEXIS_PAPER` workspace, rather than the GitHub-friendly folder structure. To avoid silently editing research code, use the separate, generated workspace:

```bash
python3 tools/prepare_workspace.py --output ./local_workspace
export MATPOWER_DIR=/path/to/matpower
cd local_workspace
octave --no-gui --quiet --eval "stage2_reopt_69_octave('check')"
octave --no-gui --quiet --eval "audit_stage2_69_losses_octave()"
octave --no-gui --quiet --eval "audit_stage2_69_final_comparisons()"
```

The MATPOWER package, its `runpf`, `mpoption`, and `loadcase` functions, and a compatible MATLAB/Octave installation are not included. MATPOWER: https://matpower.org/ . The check mode replays the selected *historical* controls with a corrected 1.00-MW inverter rating; its expected `J_V` is 1.580732 with 51 violating bus-time samples. The subsequent read-only audits evaluate the archived updated Stage-2 schedule (`J_V` 1.548854, 51 violations, 1797.509 kWh daily active losses). Not all scripts have been independently rerun in a newly provisioned environment.

For an independent *new* Stage-2 optimization, create a separate empty workspace:

```bash
python3 tools/prepare_workspace.py --output ./local_optimization --mode optimize
export MATPOWER_DIR=/path/to/matpower
cd local_optimization
octave --no-gui --quiet --eval "stage2_reopt_69_octave('optimize')"
```

This computation is substantial (27 steps, 40 wolves x 50 iterations). It preserves Stage-1 decisions and uses the older inverter parameters as a warm start. The script writes new results only into the generated workspace. Do not run it directly from the downloaded source folders; historical scripts may write outputs in-place.

## Key recorded outcomes

| Feeder / condition | Voltage-deviation index `J_V` | Violating bus-time samples |
|---|---:|---:|
| 33-bus selected V4 configuration | 1.5771 | 16 |
| 69-bus historical control *replay* with 1-MW inverter model | 1.580732 | 51 |
| 69-bus new Stage-2-only optimization (retained Stage 1) | 1.548854 | 51 |

For the final 69-bus replay: active losses 1797.509 kWh/day, reactive losses 925.767 kVArh/day, 27 accepted and 0 rejected inverter adjustments, 158 nonzero Q commands and 143 saturation events. See the saved CSVs for full precision and scenario definitions. The 0.95–1.05 pu voltage constraints are not completely met: 51 undervoltage samples remain.

## Essential method distinction

The historical 69-bus stored control parameters were originally obtained with a legacy inverter rating of 0.20 MW/site despite production data reaching 1.00 MW/site. The V4 summary was computed by replaying those decisions with a 1.00-MW/site model. The final 1.548854 result came from an additional 27-step inverter-only GWO optimization at 1.00 MW/site while leaving Stage-1 discrete controls unchanged. It is not a new full two-stage seed-42 optimization.

The five-seed GWO helper repairs infeasible radial candidates by resampling whereas the historical main GWO penalizes them; therefore these five runs are not strict identical-implementation optimizer-repeatability trials. See [docs/PROVENANCE_AND_LIMITATIONS.md](docs/PROVENANCE_AND_LIMITATIONS.md).

## Reuse, citation and release status

Please cite the associated manuscript, once publication information is available. A draft citation record is in `CITATION.cff` (no DOI claimed). Third-party MATPOWER code is not redistributed. No software license has been granted in this package. Authors must choose a license or explicitly define reuse permissions before public release, and confirm rights in the modified benchmark network and supplied profiles. See [docs/PUBLIC_RELEASE_CHECKLIST.md](docs/PUBLIC_RELEASE_CHECKLIST.md).
