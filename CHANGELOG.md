# Changelog

## 1.2.0 — November 2026 (unreleased)

Accompanies version 4 of the WHALEFALL preprint (`paper/zenodo_whale/main-v4.tex`), which adds a second realisation of
the birth with Einstein gravity and a cuscuton and narrows the wording of several statements, without changing any
result of version 3.

- U57–U65: the cuscuton birth (background, inflation, linear scalars, anisotropy, curvature and imprint, Bogoliubov
  coefficients, the CMB against Planck 2018, cubic order, tensors and vectors), with the specification, code, result
  file and report of each test; `chain/explora_cuscuton.py`, the background integrator they import; the reference
  result of the exploratory run C1 (`data/explorador/cuscuton_c1/`) and the Planck 2018 TT spectrum
  (`data/planck2018/`). The intermediate rows of the long runs are not exported.
- `paper_numbers_v4.py` and `figures/make_figures_v4.py` (fig_cuscuton); `make numbers` checks v4 and v3, and
  `make paper` builds v4. Version 3 (`main-v3.tex`) stays.

## 1.1.0 — 26 October 2026

Accompanies version 3 of the WHALEFALL preprint (`paper/zenodo_whale/main-v3.tex`), which adds the coherence of the birth.

- U27: coherence of the direct (non-eternal) children with the pole plateau, in the dust model of the parent, with the
  original code of that model (`components/c45j_ltb/`).
- U28: the parent evolved through its bounce with dust and an inhomogeneous inflaton (`chain/u28_comoving.py`,
  `chain/u28_c45k.py`), with its validation gates and the 780 runs (`results/u28_runs/`); the double-null code its
  equations come from is in `components/c45b_dn/`.
- The contract's link 2 (D11) carries the verdicts of U27 and U28.
- Version 2 of the preprint (`main-v2.tex`, `paper_numbers_v2.py`) is no longer exported; it stays in 1.0.0 and on Zenodo.

## 1.0.0 — 25 October 2026

Initial public release, accompanying version 2 of the WHALEFALL preprint (`paper/zenodo_whale/main-v2.tex`).

- Chain S1–S6 and its contract, with the pole plateau (default) and the ten-level hierarchy of version 1.
- String–wall collapse on the lattice (U19–U19d) with every network; summaries of the emission tests (U2–U25).
- CLASS and Lyman-α pipelines, and the full-likelihood harness (U15b–U15d) with its summaries.
- Preregistered specification and report of each test.
- Scripts that recompute every result number of the preprint and its figures from the stored results.
