# Data provenance and reuse

The CSV supplements and frozen analysis tables are study-specific reconstructions, joins and summaries of published observations. They include source identifiers and sequences needed to inspect product and guide correspondence. They are not new experimental measurements. Source observations and reference sequences retain their upstream attribution and terms; the software licence does not replace those terms.

Sources: Prensner et al. (2021), https://doi.org/10.1038/s41587-020-00806-2 ; Deutsch et al. (2026), https://doi.org/10.1038/s41586-026-10459-x ; Schlesinger et al. (2025), https://doi.org/10.1016/j.isci.2025.111884 ; GENCODE Phase-1 catalogue and GENCODE v50 annotations, https://www.gencodegenes.org/pages/riboseq_orfs/ and https://www.gencodegenes.org/human/ ; hg38 reference genome, UCSC.

`provenance/EXTERNAL_SOURCE_FILES.json` pins the original workbooks and author table by URL, size and SHA-256. They are not redistributed. The source-level scripts identify additional genome/annotation inputs and their recorded hashes. Users should retrieve provider materials under the provider's terms. Frozen-input statistical replay does not require these larger downloads.

Source authors performed molecular follow-up beyond the pooled outcomes reanalysed here. This repository neither re-evaluates every source experiment nor identifies the mechanism responsible for the observed guide responses. Synthetic outputs are explicitly identified in Tables S33-S34 and the efficacy-control directory.
