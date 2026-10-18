# Planar six-body central configurations: Part I

This repository contains the manuscript **Generic finiteness of planar six-body central configurations: Part I, exceptional-mass first limits** by Zhe Chen and Zhicheng Tong, plus the authors' supporting verification code and selected generated certificates.

- [Manuscript (PDF)](Smale6-CT-v3.pdf)
- [LaTeX source](Smale6-CT-v3.tex)
- [Verification code and replay instructions](verification/README.md)
- [References and provenance](CITATIONS.md)

The manuscript proves a first-limit confinement statement and generic finiteness outside an exceptional algebraic mass locus, subject to the cited Chang--Chen classification and the stated computational certificates. It does **not** claim finiteness on every exceptional positive mass vector; that is left to Part II. This manuscript is a preprint and has not been peer reviewed. The presence of replay code does not mean that every analytic step has been formally verified.

The source here is byte-for-byte identical to the author's `Smale6-CT-v3.tex` dated 18 October 2026 (SHA-256 `4cfa49c9999d29f7a84765bfa73070396c0c2ebde673616c06e4f120bd7dc055`). The PDF was built from that source with MiKTeX pdfLaTeX in three passes; it has 60 A4 pages. The paper itself was not edited for this release.

## Build the PDF

With a TeX distribution providing `pdflatex`, run from the repository root:

```powershell
pdflatex -halt-on-error -interaction=nonstopmode Smale6-CT-v3.tex
pdflatex -halt-on-error -interaction=nonstopmode Smale6-CT-v3.tex
pdflatex -halt-on-error -interaction=nonstopmode Smale6-CT-v3.tex
```

## Scope and provenance

The `verification/independent_review/appendix_replay_manifest.json` file maps the manuscript appendix to the implemented checks and records each check's limits. The published Chang--Chen 117-diagram classification and its symbolic computations are cited inputs; the included scripts do not reproduce their entire classifier.

The `generate_117_diagrams.py` script can redraw the 117 cached diagrams from the Chang--Chen author-distributed Mathematica worksheet. That third-party worksheet is **not** redistributed here. Obtain it from the [authors' supplement page](https://www.math.nthu.edu.tw/~kchen/SSCC.html) and place `SingSeqCCver202303.nb` in `verification/` before running the redrawer. The redraw checks cached graph data; it does not rerun Algorithm I.

No licence is granted for the manuscript or code by this repository's publication alone. For reuse beyond applicable law, contact the authors.
