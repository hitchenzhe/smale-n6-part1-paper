# Verification code

The scripts under `independent_review/` replay specific algebraic, modular and interval calculations cited in Appendix A of the manuscript. Read [`appendix_replay_manifest.json`](independent_review/appendix_replay_manifest.json) for the exact scope and limitations of each check. The included `*_certificate.json` files are generated certificate snapshots.

This release omits the local Python virtual environment, machine-specific run logs, exploratory audit notes and the third-party Chang--Chen Mathematica worksheet. The certificate generators and verifiers remain available.

## Windows replay

From the repository root, with Python 3.12 available:

```powershell
.\verification\independent_review\setup_python.ps1
.\verification\independent_review\run_independent_review.ps1
```

`setup_python.ps1` creates `verification/.venv` and installs `requirements.lock.txt`. The replay wrapper runs the listed jobs sequentially with a 1,024 MB memory limit and a 900 second timeout per job. Generated logs and summaries stay local and are ignored by Git.

`verify_n6_structure.py` checks the 117-label accounting and the seven-product exponent matrix. `generate_117_diagrams.py` requires the separately obtained Chang--Chen worksheet described in the root README and is not part of the standard replay.

The scripts check their stated computations. They do not independently certify the cited 117-diagram classification or all analytic and geometric arguments in the manuscript.
