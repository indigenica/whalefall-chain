# whalefall-chain: code, frozen inputs, results and preregistrations behind the WHALEFALL preprint, version 2.
# Generated from indigenica/whalefall (see EXPORT_PROVENANCE.json). Do not edit here: edit whalefall and re-export.
PY ?= python3
export PYTHONDONTWRITEBYTECODE = 1
export OMP_NUM_THREADS ?= 2

.PHONY: all test chain chain-sum reproduce action numbers figures paper env-check network-validate network-u2 network-u3 clean

all: test reproduce action numbers figures ## everything that runs in ~15 min on a laptop (what CI runs)

test:            ## fast tests: chain pieces, small lattice runs, Boltzmann solver, likelihood-chain tools (~3 min; gcc + OpenMP)
	$(PY) -B tests/run_tests.py

chain:           ## S1-S6 + contract + figures from data/ and network/results (~4-6 min); pole plateau (default since U26)
	$(PY) -B chain/run_all.py

chain-sum:       ## the same chain with the ten-level hierarchy of version 1 (writes results/plateau_sum/)
	WHALEFALL_PLATEAU=sum $(PY) -B chain/run_all.py

reproduce: chain chain-sum ## both chains, then check that they reproduce the stored results
	$(PY) -B chain/check_reproduction.py results
	$(PY) -B chain/check_reproduction.py results/plateau_sum

action:          ## common action: symbolic residuals, decays, late fluid (needs sympy)
	mkdir -p action/logs action/results action/logs_pole action/results_pole action/runtime && cd action && $(PY) -B src/run_all.py

numbers:         ## every result number quoted in main-v2.tex, recomputed from the stored results (exit 1 on any mismatch)
	$(PY) -B paper/zenodo_whale/paper_numbers_v2.py

figures:         ## preprint figures from the stored results (fig_collapse, fig_band, fig_lyafull, fig_bracket)
	$(PY) -B paper/zenodo_whale/figures/make_figures_v2.py
	$(PY) -B -c "import sys; sys.path.insert(0, 'paper/zenodo_whale/figures'); import make_figures as M; M.fig_bracket()"

paper: figures   ## paper/zenodo_whale/main-v2.pdf (needs latexmk and a TeX distribution)
	cd paper/zenodo_whale && latexmk -pdf -interaction=nonstopmode main-v2.tex

env-check:       ## checks numpy, sympy, gcc+OpenMP and, for the Lyman-alpha tests, classy with the quadrature patch
	$(PY) env/check_env.py

# ---- production (hours of CPU; see README) ----
network-validate: ## kernel validations V1-V4 (~1 min)
	cd network && $(PY) -B validate.py

network-u2:      ## the 20 networks of U2 (256^3; ~3 h with 2 cores, ~1.5 GB) + analysis
	cd network && for S in D A_U1 A_RP3a E B C; do $(PY) -B produce.py $$S; done && $(PY) -B analyze.py

network-u3:      ## U3b core-preparation networks (~2 h) + gates G1/G2
	cd network && for S in F1_U1 F2 F1_RP3; do $(PY) -B produce_u3.py $$S; done && $(PY) -B analyze_u3.py

clean:
	rm -rf components/c45b_hom/_build network/_build results/cadena_rerun logs action/logs action/results action/logs_pole action/results_pole action/runtime
	cd paper/zenodo_whale && rm -f *.aux *.log *.fls *.fdb_latexmk *.out *.bbl *.blg
