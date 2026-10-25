# WHALEFALL: the chain behind the preprint

[![Preprint DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23478896.svg)](https://doi.org/10.5281/zenodo.23478896)
[![reproducibility](https://github.com/indigenica/whalefall-chain/actions/workflows/ci.yml/badge.svg)](https://github.com/indigenica/whalefall-chain/actions/workflows/ci.yml)

Code, frozen inputs, results, preregistrations and paper source for

> J. Bautista, *WHALEFALL: a preregistered, falsifiable toy cosmology from a quasi-topological bounce to cold dark matter
> from the collapse of an RP³ string–wall network*, version 2, preprint (2026), Zenodo, [doi:10.5281/zenodo.23478896](https://doi.org/10.5281/zenodo.23478896).
> Source in `paper/zenodo_whale/main-v2.tex`, PDF in `paper/zenodo_whale/main-v2.pdf`.

WHALEFALL is a toy cosmology built as one chain of computations, from a non-singular beginning to a dark-matter
prediction. A small closed contracting patch bounces through a quasi-topological gravity core and kicks a scalar onto an
order-four pole plateau (n_s = 0.9728, r = 0.0050). A symmetric-tensor order parameter orders into an RP³ vacuum and leaves
a network of global ℤ₂ strings, which collapses as string–wall systems when H ≃ m_G/3. The collapse, measured on the
lattice and passed through a patched CLASS, Lyman-α gates and a full Lyman-α likelihood, leaves a band of the decay
constant F ≈ 0.06–0.5 F₀ (F₀ = 1.07×10¹⁵ GeV) with cold pseudo-Goldstone dark matter, m_G ≈ 10⁻¹⁷–10⁻¹³ eV.

Every computation that yields a verdict was preregistered, with its thresholds and verdict logic, before it was run.
Amendments are dated and failures are reported as such. Appendix A of the preprint lists every test by its code (U1,
U2, …); the code indexes this repository.

## What is here, and what is not

This repository is a **one-way export** of the working repository (`indigenica/whalefall`, private) with only what is
needed to reproduce the preprint: the code, the frozen inputs, the results, and the preregistered specification and report
of each test. It carries no git history, superseded material, internal notes, audits or reviews; the audits, and the
external reviews with their answers, are available from the author on request. `EXPORT_PROVENANCE.json` gives the source
commit and the SHA-256 of every file.

The RP³ and U(1) **string-network emission** measurements (U2–U25: density ratio R_ξ, emission ratio R_Γ, spectral
slope q, resolution, drain and time step) are the subject of a companion paper with its own repository,
[indigenica/whalefall-rp3-string-networks](https://github.com/indigenica/whalefall-rp3-string-networks)
([doi:10.5281/zenodo.23451009](https://doi.org/10.5281/zenodo.23451009)). Here they enter through the summary of each
test (`network/results/*/U*_RESULT.json`), which is what the chain and its contract read. The networks one by one are
there, except those of U2 and U3 (used by `make test`) and the 576³ networks drawn in `fig_bracket.pdf`. Re-running
those analyses needs the per-network files of the companion repository, which also holds the independent audit of the
emission estimators (`audit/u4_u11/`; its summary `a1_network.json` is kept here because the contract cites it).

The **string–wall collapse** (U19–U19d) is complete here: code, every network, validations and run logs.

## Reproduce

Requirements: Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib` and `sympy` (`pip install -r requirements-action.txt`),
and `gcc` with OpenMP (the lattice kernel and the homogeneous birth integrator are compiled on first use). TeX Live with
`latexmk` only for the PDF.

```bash
make test       # fast tests: chain pieces, small lattice runs, Boltzmann solver, likelihood-chain tools (~3 min)
make reproduce  # chain S1-S6 + contract with the pole plateau and with the ten-level hierarchy of version 1,
                # then check that both reproduce the stored results (rtol 1e-6; ~10 min)
make action     # common action: symbolic residuals, tree-level decays, late fluid (sympy)
make numbers    # every result number quoted in main-v2.tex, recomputed from the stored results
make figures    # preprint figures from the stored results
make paper      # paper/zenodo_whale/main-v2.pdf
```

`make all` runs `test reproduce action numbers figures` (what CI runs, ~15 min on a laptop). Without `gcc` some tests
are reported as `SKIP`, never as passed; CI sets `REQUIRE_ALL=1`, so there a skipped test fails.

**Heavy computations** (not run by `make`; each script states its command, cost and outputs in its docstring):

- **Lattice.** `network/kern.c` (C + OpenMP) with `network/lattice.py`; an optional JAX port (`network/jax_kernel.py`,
  `WHALEFALL_BACKEND=jax`) was used on an A100 GPU for the largest collapse networks. The collapse networks
  (`network/u19*_produce.py`, resumable, skip finished networks) took 0.4–2.6 h each on a 4-core machine or the GPU at
  256–576³ and 7–11 h each at 512³ on a 44-thread workstation (~30 GB of RAM per RP³ network at 512³).
  `make network-u2` and `make network-u3` reproduce the 256³ networks of U2 and U3b in a few CPU hours.
- **CLASS.** CLASS v3.4.1 (commit `814ce59`, `env/CLASS_VERSION.txt`) with the quadrature patch
  `env/class_qm_trapz_log.patch`; `env/cloud_setup.sh` builds it and `make env-check` verifies it. The Lyman-α tests are
  `lya/u*_class.py`; each CLASS row takes 20–40 min on one core.
- **Full Lyman-α likelihood.** MontePython with Planck lite, BAO and the Lyman-α likelihood of `lya/u15b/` (tables in
  `lya/u15b/tables/`, harness and chain tools in `lya/u15b/` and `lya/u15c/`; setup in `env/INSTALAR_WSL.md`). The
  chains took several days per campaign on the workstation; they are kept by the author, and their summaries are in
  `lya/results/u15*/`.

## Layout

```
chain/              S1 birth, S2 pocket, S3 sectors, S4 dark ledger, S5 closure, S6 measured closure, the contract
                    (u_contract.py) and the dark-sector tests run inside the chain (U17, U20, U20b, U22, U22b, U23, U26)
components/         frozen components from earlier campaigns (birth integrator, hierarchy, relic, pseudo-Goldstone)
data/               frozen inputs, with SHA-256 manifest (data/MANIFEST.json)
results/            chain outputs and the contract (CONTRATO_CIERRE.json); results/plateau_sum/ with the hierarchy
action/             common action (Q + χ + radiation + quasi-topological gravity): symbolic checks, decays, late fluid
lya/                Boltzmann solver (boltz.py), CLASS pipelines, the full-likelihood harness (u15b/, u15c/) and results
network/            lattice kernel and the string–wall collapse (U19-U19d) with every network; emission summaries
docs/               preregistered specifications (preregistros/ESPEC_U*.md, index in its README) and the report with
                    the verdict of each test (informes/INFORME_U*.md)
env/                environment check, CLASS patch and setup
paper/zenodo_whale/ the preprint (main-v2.tex, refs.bib), figure scripts and paper_numbers_v2.py
tests/              fast tests (make test)
```

The preregistrations and reports are in Spanish: they are the original dated records and are never rewritten
(corrections are dated notes appended to them). The preprint, this README and the code comments are in English.

**Plateau.** Since U26 the chain uses the order-four pole plateau by default; `WHALEFALL_PLATEAU=sum` selects the ten-level
hierarchy of version 1 and writes to `results/plateau_sum/` (chain) or `action/results/` (common action).

## Licence

Code under GPL-3.0-or-later; text, figures and data under CC-BY-NC-4.0. See `LICENSE.md`.
