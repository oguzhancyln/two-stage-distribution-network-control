# Checklist before making the repository PUBLIC

- [ ] Obtain all co-authors' explicit permission to publish the included source code, raw inputs, control schedules, and supplementary outputs.
- [ ] Confirm redistribution rights and original authorship for the modified 69-bus feeder and any network-derived inputs; MATPOWER itself is not packaged.
- [ ] Choose and add a `LICENSE` file if external reuse is intended. The present ZIP intentionally does not claim a license or grant third parties rights.
- [ ] Review each historical algorithm variant and source file for sensitive information, attribution, and unintended dependencies.
- [ ] Confirm that the associated manuscript and response to reviewers describe historical seed-42 vs new inverter-only Stage-2 calculations accurately.
- [ ] Run `python3 tools/verify_repository.py` on the GitHub clone. Optionally run the MATPOWER replay commands in README on a clean machine.
- [ ] After accepted/publication, update `CITATION.cff` with the correct final journal citation and DOI. Do not add an unpublished DOI placeholder.
- [ ] Decide whether to archive a numbered release on Zenodo for a persistent version-specific DOI (optional).
- [ ] Check for the target journal's peer-review confidentiality and supplementary-code policies; avoid uploading non-public decision letters and reviewer identities.

Repository archive includes no manuscript PDF, review reports, or private correspondence. It contains research code, raw fixed numerical matrices and saved experiment outputs.
