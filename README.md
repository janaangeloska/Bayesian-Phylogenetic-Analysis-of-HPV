# Bayesian Phylogenetic Analysis of HPV-16

A reproducible, end-to-end pipeline for Bayesian phylodynamic inference of Human papillomavirus type 16 (HPV-16). The pipeline covers sequence retrieval, quality filtering, multiple sequence alignment, alignment validation, substitution model selection, BEAST configuration, MCMC convergence assessment, and phylogenetic tree analysis.

**Authors:** Oliver Buteski, Jana Angeloska, Rebeka Maneva, Biljana Tojtovska Ribarski PhD  
**Institution:** Faculty of Computer Science and Engineering, Ss Cyril and Methodius University, Skopje, N. Macedonia

---

## Overview

Starting from 13,161 raw HPV-16 nucleotide sequences retrieved from NCBI Virus, the pipeline produces a curated dataset of 528 complete or near-complete genomes (spanning 2001-2024) and submits them to Bayesian inference in BEAST v2.7.7 under the HKY+G+I substitution model with a Coalescent Exponential Population prior. Key findings are summarized in the Results Summary section below.

---

## Dependencies

### Python

Python 3.13 is required.

| Package    | Version |
|------------|---------|
| biopython  | 1.87    |
| matplotlib | 3.10.9  |
| numpy      | 2.4.6   |
| pandas     | 3.0.3   |
| scipy      | 1.17.1  |

### External Tools

These tools must be installed separately and available on your system PATH (or invoked directly by path).

| Tool            | Version  | Purpose                                      |
|-----------------|----------|----------------------------------------------|
| MAFFT           | 7+       | Multiple sequence alignment                  |
| MEGA            | 12       | Substitution model selection (AIC/BIC)       |
| BEAST           | 2.7.7    | Bayesian phylogenetic inference              |
| BEAUti          | 2.7.7    | BEAST XML configuration (bundled with BEAST) |
| Tracer          | 1.7      | MCMC convergence diagnostics                 |
| TreeAnnotator   | 2.7.7    | MCC tree summarization (bundled with BEAST)  |
| FigTree         | 1.4.4    | Phylogenetic tree visualization              |
| TempEst         | 1.5.3    | Temporal signal assessment (recommended)     |

---

## Repository Structure

```
.
├── data/          # Sequence data at each pipeline stage (raw, filtered, aligned, MCC tree)
├── notebooks/     # Jupyter notebooks 01–04 (preprocessing → alignment → validation → visualization)
├── scripts/       # Post-BEAST Python scripts (temporal signal, clade assignment, population dynamics, etc.)
├── results/       # Output figures and summary files
└── literature/    # Reference PDFs
```

---

## Pipeline

The analysis follows these sequential steps:

1. **Data retrieval** - Download HPV-16 sequences and metadata from NCBI Virus
2. **Preprocessing** (`notebooks/01`) - Parse FASTA, join metadata, filter sequences lacking collection dates, standardize date formats
3. **Quality filtering + alignment** (`notebooks/02`) - Filter by length (>=7000 bp), ambiguity fraction (<5% N), and duplicates; align with MAFFT `--auto`
4. **Alignment validation** (`notebooks/03`) - Quantitative checks: gap fraction, GC content consistency, column conservation, pairwise divergence, temporal signal proxy
5. **Model selection** (MEGA 12) - Evaluate five candidate substitution models by AIC/BIC on the aligned dataset
6. **BEAST configuration** (BEAUti) - Configure tip dates, substitution model (HKY+G+I), strict clock, Coalescent Exponential Population prior, and MCMC settings; export XML
7. **Bayesian inference** (BEAST v2.7.7) - Run five parallel BEAST analyses (one per candidate model) for cross-model convergence comparison
8. **Convergence assessment** (Tracer v1.7) - Inspect ESS values and trace plots; confirm stationarity
9. **Tree summarization** (TreeAnnotator) - Generate MCC tree with 10% burn-in and median node heights
10. **Post-processing** (`scripts/`) - Temporal signal, TMRCA by clade, lineage assignment, phylogenetic diversity, branch outlier detection, population dynamics

---

## Reproducing the Analysis

```bash
# 1. Clone the repository
git clone https://github.com/janaangeloska/Bayesian-Phylogenetic-Analysis-of-HPV.git
cd Bayesian-Phylogenetic-Analysis-of-HPV

# 2. Set up the Python environment
python -m venv .venv
.venv\Scripts\activate       # Windows
source .venv/bin/activate    # macOS/Linux
pip install biopython==1.87 matplotlib==3.10.9 numpy==2.4.6 pandas==3.0.3 scipy==1.17.1

# 3. Run the notebooks in order (01 through 04)
jupyter notebook

# 4. Run MEGA model selection on data/hpv_aligned_mafft.fasta (GUI)

# 5. Configure and run BEAST using BEAUti + the exported XML

# 6. Run post-processing scripts after BEAST completes
python scripts/00_make_lineage_meta.py
python scripts/01_temporal_signal.py
python scripts/02_clade_assignment.py
python scripts/03_tmrca_by_clade.py
python scripts/04_phylogenetic_diversity.py
python scripts/05_branch_outliers.py
python scripts/06_population_dynamics.py
```

Results are written to `results/`.

---

## Results Summary

| Parameter | Estimate | 95% HPD / CI |
|-----------|----------|--------------|
| TMRCA | 1709 AD | 1609–1797 AD |
| Tree height | 315.5 yrs | 235.6–428.1 yrs |
| Clock rate (μ) | 4.04 × 10⁻⁵ subs/site/yr | 2.87×10⁻⁵ – 5.34×10⁻⁵ |
| Growth rate (r) | 0.023 yr⁻¹ | 0.016–0.031 |
| Ne at present | 6,812 (coalescent units) | — |
| κ (kappa) | 3.32 | 3.09–3.56 |
| p_inv | 0.713 | 0.692–0.733 |
| α (gamma) | 0.539 | 0.453–0.637 |

### MCMC Convergence (HKY+G+I, best model)

| Parameter | ESS |
|-----------|-----|
| treeLikelihood | 20,039 |
| gammaShape | 4,028 |
| proportionInvariant | 3,555 |
| posterior | 1,776 |
| growthRate | 312 |
| ePopSize | 304 |
| Tree.height | 278 |
| clockRate | 256 |

All parameters above the ESS ≥ 200 adequacy threshold.

---

## Important Caveats

- **Temporal signal:** 94.3% of tips carry year-only imputed dates (July 1st), creating a binning artefact that inflates root-to-tip r to 0.998. This does not reflect genuine per-sequence clock signal. TempEst analysis on an ML tree is recommended before treating clock estimates as fully data-driven.
- **Strict clock:** All results assume a constant substitution rate across lineages. A relaxed lognormal clock is the recommended next step.
- **Preliminary convergence:** Tree.height and clockRate ESS values (~256–278) are marginal. Results should be treated as preliminary pending longer MCMC chains.
- **Sampling bias:** The 528 sequences are biased toward complete, well-annotated genomes from high-income countries. Sub-Saharan Africa and Southeast Asia are underrepresented.

---

## Citation

If you use this pipeline, please cite the key tools:

- Katoh & Standley (2013) - MAFFT v7, *Mol. Biol. Evol.* 30(4):772-780
- Kumar et al. (2018) - MEGA X, *Mol. Biol. Evol.* 35(6):1547-1549
- Bouckaert et al. (2019) - BEAST 2.5, *PLoS Comput. Biol.* 15(4):e1006650
- Rambaut et al. (2016) - TempEst, *Virus Evol.* 2(1):vew007