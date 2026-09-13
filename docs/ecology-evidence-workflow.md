# Ecological evidence and transition curation

Ecological tests require species-level evidence and independently supported transitions. A database match to a genus is useful for prioritizing curation, but it is not a validated species state, a negative result for unrecorded traits, or evidence of an independent evolutionary transition.

## Genus-level candidate catalogue

`import_ecology_candidates.py` pins the publisher-hosted FungalTraits supplementary workbook by SHA256 and matches the first word of the selected fungal name to an exact `GENUS` label. It preserves primary and secondary lifestyles, other source trait strings, source rank, spreadsheet row, source taxonomy and links to taxon-identity flags. It neither resolves synonyms nor fills missing trait values. Duplicate source-genus rows are treated as ambiguous rather than choosing one. Non-fungal outgroups are not assigned fungal traits.

The retrieved workbook is 990,658 bytes and contains 10,771 genus rows with 10,765 distinct labels and six duplicate labels. These are observed snapshot contents, not a claim that this download reproduces the historical row counts in the article. The source is [Põlme et al., FungalTraits](https://doi.org/10.1007/s13225-020-00466-2), publisher Supplementary material 4, `data` worksheet. The paper's volume year is 2020 and online publication date is 19 January 2021. The source snapshot and interpreter version are recorded in the receipt.

The working catalogue has exact, unambiguous genus candidates for 482 fungal entries, no exact match for 19 fungal entries, and 25 unassigned outgroups. All candidate rows have `species_verified=False` and are ineligible for confirmatory tests. Source comments are flagged for consultation in the workbook; they are not replaced by guessed numerical confidence. The main analysis manifest's ecology fields are not automatically overwritten.

## Species-level literature evidence

`config/ecology_species_evidence.json` contains 21 manually reviewed species statements (17 ectomycorrhizal, two asymbiotic and two saprotrophic). The original five Amanita statements derive from the focal-species classification on page 2787 of [Hess et al. 2018](https://doi.org/10.1093/molbev/msy179), and the original Laccaria statement from the primary [Laccaria bicolor genome study](https://doi.org/10.1038/nature06556). `build_species_ecology_evidence.py` requires each statement to match exactly one selected fungal entry, attaches assembly identity, and records disagreement with genus-level candidates.

Two projections require correction: the exact-name genus match gives an ectomycorrhizal candidate for Amanita inopinata and A. thiersii, whereas the species study classifies them as asymbiotic. This is a conflict in genus-to-species transfer, potentially involving genus circumscription and name changes; it is not evidence that the source database itself is wrong. “Asymbiotic” is retained as the source's classification, without inferring a particular decay substrate or absence of every kind of interaction.

The Amanita muscaria statement specifically concerns var. guessowii; equivalence to the selected assembly's taxonomic concept remains flagged. Other species classifications also do not establish direct experimental validation of the selected isolate. The sampled symbiotic Amanita represent one origin in the source study and must not be counted as independent replicate origins. Transition group labels are provisional curation aids, not final mapped branches.

## Reproduction and next evidence gates

```bash
python scripts/import_ecology_candidates.py
python scripts/build_species_ecology_evidence.py
python scripts/build_suillus_host_evidence.py
python scripts/build_ecology_structure_coverage.py
```

The workbook stays outside Git under `data/traits/`; curated configurations, candidate/evidence tables and hash receipts are versioned. The candidate importer uses openpyxl 3.1.5. Tests check ambiguous names, missing matches, outgroup handling and genus-only status. The species evidence table distinguishes published classifications from exact-isolate validation and leaves confirmatory-test status pending taxonomy, phylogeny and replicated-transition review.

Continue species and strain curation using primary descriptions and experiments, preserve mixed or context-dependent states, and match accepted names explicitly. Include source conflicts and uncertainty in trait-coding sensitivities. Count evolutionary replication on supported trees, not by the number of related genomes sharing a state. The 21 statements do not complete the project's ecological coverage or its association tests.

## Expanded species and host evidence

[Lofgren et al. 2021](https://doi.org/10.1111/nph.17160), Table 1 (journal p. 776), supplies nine selected Suillus classifications and Laccaria amethystina. [Kohler et al. 2015](https://doi.org/10.1038/ng.3223), p. 413, supports Piloderma croceum. [Peter et al. 2016](https://doi.org/10.1038/ncomms12662), abstract and genome-comparison sections, supports ectomycorrhizal Cenococcum geophilum and its saprotrophic comparators Glonium stellatum and Lepidopterella palustris. [Martin et al. 2010](https://doi.org/10.1038/nature08867), abstract, supports Tuber melanosporum. These are reviewed published classifications, with DOI and passage locators in the configuration; no local PDF checksum is claimed.

Piloderma requires name review: the 2015 study uses P. croceum, while Lofgren's table uses P. olivaceum for project Pilcr1. The configuration preserves this discrepancy without assuming synonymy or selected-isolate equivalence. Cenococcum strain 1.58 is explicit in the source methods, but its relationship to the selected assembly still requires verification.

The separate Suillus host table covers all ten selected entries: eight have reported host classifications, S. weaverae has explicitly uncertain host assignment, and S. discolor lacks a verified matching Table 1 assignment. The source's red-pine category denotes Pinus subgenus Pinus, not the single species P. resinosa. Multiple hosts remain categorical sets. The genomics study compiles host literature; these rows do not establish colonization experiments for our isolates. No number of independent host switches is inferred from tip counts.

## Predictor-specific coverage gate

The frozen AlphaFold and ESMFold paired-input receipts are joined separately to curated taxa. All 21 have an eligible marker in at least one source. However, the five-taxon Amanita group shares zero eligible markers across every member within either source. The three-taxon Cenococcum/comparator group shares one AlphaFold marker (4974767at2759) and zero ESMFold markers. Cenococcum has 1 AlphaFold/28 ESMFold eligible markers, versus 118/0 for Glonium and 122/0 for Lepidopterella.

These stringent group intersections identify missing matched-source coverage; they do not rule out every pairwise comparison or quantify power. Prediction method is strongly confounded with ecological state in the current Cenococcum contrast. Targeted same-method prediction and renewed confidence filtering are needed before ecological inference, alongside orthology, taxonomy and phylogenetic replication checks. Coverage receipts pin both source inputs, both curated evidence tables, the builder and output tables.

## Same-method prediction inputs

`python scripts/prepare_ecology_prediction_inputs.py --output data/prediction_inputs/ecology-markers-v1`
freezes all eligible missing ESMFold marker sequences across the 21 curated
species, including proteins already represented by AlphaFold. The immutable queue
contains 675 full canonical proteins (212,983 residues) across 18 species. Another
739 sequences are reserved in existing prediction queues; this is not a claim of
completed predictions. There are 1,040 length-deferred and 27 noncanonical
sequences, retained explicitly without truncation or residue substitution.

Source receipts and full FASTA readback pass, and the new queue is disjoint from
all three existing input queues. This queue is prepared but not launched: both
GPUs currently run earlier batches. A resource estimate and device-availability
check precede execution. More predictions do not guarantee confidence-qualified
coverage, verified orthology or ecological replication.
