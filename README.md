# Guide selection across noncanonical ORF screens

Data and analysis pipeline for **Guide selection transfers depletion effects across noncanonical ORF screens**, by Robin Gounder (University of Cambridge) and Samuel Schmidt (Vaionex Corporation).

The study matches complete reference ORFs and exact CRISPR guide identities across published human fitness screens. It quantifies persistent within-product guide responses and a retrospective guide-selection consequence. It does not establish a new peptide mechanism or medical target. Shared on-target efficacy acting on true dependencies can generate both patterns; the synthetic positive control demonstrates compatibility and does not fit an explanation of the observed data.

## Reproduce the released analysis

Use Python 3.11 in a fresh virtual environment:

```bash
python -m pip install -r requirements.txt
python reproduce.py --output results/reproduction
```

This verifies input hashes, reruns the shared-guide and guide-selection analyses, independently reconstructs the additional-assay and third-publication statistics, reruns the context-exclusion and synthetic controls, and rebuilds all eight figures. Outputs and logs are written to the selected output folder; released data remain unchanged. The first run should take a few minutes. Arial produces the supplied artwork typography; another installed font may change rendering, while numerical reproduction is unchanged.

To also regenerate the full 2,000-cluster by 1,000-guide-plan balanced resampling calculation:

```bash
python reproduce.py --output results/full --full-resampling
```

The default run recalculates its point estimates and intervals from the released fixed guide-plan statistics and bootstrap draws. The full run regenerates those draws from the frozen plans and processed guide outcomes. Guide-count balancing does not resolve the uncertain product-level comparison.

## Contents

| Path | Contents |
| --- | --- |
| `supplementary/` | All 34 CSV supplements: products, guide eligibility, outcomes, uncertainty and synthetic controls. |
| `evidence/analysis_extensions/` | Product/reference models, analysis-ready inputs, original analysis/verification drivers and frozen protocols. |
| `evidence/selection_context_sensitivity/` | Context-exclusion outputs and protocols. |
| `evidence/guide_efficacy_control/` | Fixed synthetic positive-control protocol and complete outputs. |
| `figure_source_data/`, `figures/` | Numerical figure data and vector PDF, editable SVG and 600-dpi PNG masters. |
| `code/`, `reproduce.py` | Portable frozen-input statistical replay and figure rebuild. |
| `provenance/EXTERNAL_SOURCE_FILES.json` | Original provider URLs, exact sizes, hashes and the pinned third-study commit. |

## Scope of reproduction

The tested entrypoint starts from the released analysis-ready data and frozen product/guide geometry. It does not redo genome alignment or independently establish that each reference product is translated in the screened cells. Historical source-preparation and verification drivers are supplied for transparency; running source-level reconstruction requires retrieving their pinned external sources and restoring the documented source layout. These stages are not presented as a tested, one-command raw-genome pipeline. No new model scores, private working notes, credentials or unpublished patient data are included.

Published source workbooks and the full third-study source table are retrieved from the original providers, rather than redistributed. See [DATA_PROVENANCE.md](DATA_PROVENANCE.md). This is a data/code repository; no SSRN submission is performed by these scripts. A versioned archival DOI is distinct from the GitHub URL and is not claimed here.

## Generative AI disclosure

OpenAI Codex agents contributed substantially to source discovery, literature synthesis, and mathematical exposition

## Citation and licence

See `CITATION.cff`. Original analysis code uses the MIT licence; source data retain their provider attribution and terms. Interpretive limits, fixed cohorts, nominal intervals and retrospective chronology are recorded in the frozen protocols and manuscript.
