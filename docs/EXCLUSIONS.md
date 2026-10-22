# Intentional omissions from the GitHub presentation

The source ZIPs provided with the research contain two ad hoc local utilities with computer-specific Windows home-directory paths:

- `genecheck.m`: local text search tool, not part of the claimed experiments.
- `restore_runtime_benchmark_v2.m`: file-restoration script dependent on an unavailable `.backup_pf095` file; changes another source file in-place, and is not needed for reproducing the reported results.

They are not distributed in this public-ready copy. The original research archive remains available to the authors offline. An identical duplicate of `case69modified.m` in an earlier V4 scripts folder was also omitted (its bytes match `data/networks/case69modified.m`). No included numerical input, experiment algorithm or saved result was edited for this GitHub packaging step.
