# Ecological evidence and transition curation

Ecological tests require species-level evidence and independently supported transitions. A database match to a genus is useful for prioritizing curation, but it is not a validated species state, a negative result for unrecorded traits, or evidence of an independent evolutionary transition.

## Genus-level candidate catalogue

`import_ecology_candidates.py` pins the publisher-hosted FungalTraits supplementary workbook by SHA256 and matches the first word of the selected fungal name to an exact `GENUS` label. It preserves primary and secondary lifestyles, other source trait strings, source rank, spreadsheet row, source taxonomy and links to taxon-identity flags. It neither resolves synonyms nor fills missing trait values. Duplicate source-genus rows are treated as ambiguous rather than choosing one. Non-fungal outgroups are not assigned fungal traits.

The retrieved workbook is 990,658 bytes and contains 10,771 genus rows with 10,765 distinct labels and six duplicate labels. These are observed snapshot contents, not a claim that this download reproduces the historical row counts in the article. The source is [Põlme et al., FungalTraits](https://doi.org/10.1007/s13225-020-00466-2), publisher Supplementary material 4, `data` worksheet. The paper's volume year is 2020 and online publication date is 19 January 2021. The source snapshot and interpreter version are recorded in the receipt.

The working catalogue has exact, unambiguous genus candidates for 482 fungal entries, no exact match for 19 fungal entries, and 25 unassigned outgroups. All candidate rows have `species_verified=False` and are ineligible for confirmatory tests. Source comments are flagged for consultation in the workbook; they are not replaced by guessed numerical confidence. The main analysis manifest's ecology fields are not automatically overwritten.

## Species-level literature evidence

`config/ecology_species_evidence.json` contains six manually reviewed statements. Five derive from the focal-species classification on page 2787 of [Hess et al. 2018](https://doi.org/10.1093/molbev/msy179), and one from the primary [Laccaria bicolor genome study](https://doi.org/10.1038/nature06556). `build_species_ecology_evidence.py` requires each statement to match exactly one selected fungal entry, attaches assembly identity, and records disagreement with genus-level candidates.

Two projections require correction: the exact-name genus match gives an ectomycorrhizal candidate for Amanita inopinata and A. thiersii, whereas the species study classifies them as asymbiotic. This is a conflict in genus-to-species transfer, potentially involving genus circumscription and name changes; it is not evidence that the source database itself is wrong. “Asymbiotic” is retained as the source's classification, without inferring a particular decay substrate or absence of every kind of interaction.

The Amanita muscaria statement specifically concerns var. guessowii; equivalence to the selected assembly's taxonomic concept remains flagged. Other species classifications also do not establish direct experimental validation of the selected isolate. The sampled symbiotic Amanita represent one origin in the source study and must not be counted as independent replicate origins. Transition group labels are provisional curation aids, not final mapped branches.

## Reproduction and next evidence gates

```bash
python scripts/import_ecology_candidates.py
python scripts/build_species_ecology_evidence.py
```

The workbook stays outside Git under `data/traits/`; curated configurations, candidate/evidence tables and hash receipts are versioned. The candidate importer uses openpyxl 3.1.5. Tests check ambiguous names, missing matches, outgroup handling and genus-only status. The species evidence table distinguishes published classifications from exact-isolate validation and leaves confirmatory-test status pending taxonomy, phylogeny and replicated-transition review.

Continue species and strain curation using primary descriptions and experiments, preserve mixed or context-dependent states, and match accepted names explicitly. Include source conflicts and uncertainty in trait-coding sensitivities. Count evolutionary replication on supported trees, not by the number of related genomes sharing a state. The six statements do not complete the project's ecological coverage or its association tests.
