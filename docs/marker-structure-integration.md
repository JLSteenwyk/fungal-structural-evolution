# Connecting structures to phylogenetic sites

Run `python scripts/map_marker_structures.py --output results/structural_markers/snapshot-v1` to make an immutable snapshot from currently verified AFDB models and the complete 125-marker panel. Choose a new output directory as acquisition advances.

Each linked model must have the exact complete amino-acid sequence of a marker protein. The script rechecks its CIF checksum, polymer sequence and alpha-carbon coverage. Identical sequences can reuse a model across taxa, with species and original protein identifiers retained; shared coordinates are not independent experimental observations. Where multiple candidates have identical sequences, the snapshot chooses the highest mean full-protein pLDDT, then version and model ID as deterministic tie-breaks. Candidate count and selected-model provenance are recorded. This rule is a baseline, not proof of superior biological accuracy.

The full Stockholm alignment connects each original protein residue to its BUSCO profile position. Retained profile positions then map into the concatenated species matrix. Inserted residues advance protein coordinates; gaps do not. The compressed residue table records marker, taxon, protein, model/version, matrix column, protein residue number and CA pLDDT. This makes later site-level sequence/structure comparisons traceable without treating raw structural distances as additive evolutionary branch lengths.

The first executed snapshot screened all 59,840 marker proteins and linked 151 proteins in two taxa, with 57,413 matrix-to-structure residue links. No marker had four linked taxa, so this snapshot does not support branch-specific structural inference. Confidence summaries describe the retained marker residues. PAE, domain boundaries, prediction-source effects, orthology and phylogenetic uncertainty still require assessment.

## Retrieval ordering

Current acquisition first visits marker candidates in round-robin taxon order, then processes every other exact-sequence candidate. No atlas candidates are discarded. The queue snapshot records its size and taxon coverage. Two workers and four outstanding futures keep concurrent requests bounded. New stop handling finishes the few active requests and records their receipts before exiting. Completed CIF checksums are verified on resume.

The earlier queue processed complete proteomes consecutively, delaying comparative coverage. It was intentionally terminated after preserving 6,924 valid model receipt records, then replaced only after its process was confirmed terminal. The new queue started with 12,007 priority marker accessions across 113 taxa, followed by the remaining candidates in a 1,043,924-accession snapshot. Matching continues across the full cohort, so later snapshots will add newly available candidates. The retrieval rate and full-data objective are unchanged.

The second executed snapshot, `results/structural_markers/snapshot-v2`, confirms improved comparative coverage: 452 marker proteins across 118 taxa use 422 distinct models, with 159,837 matrix-to-structure residue links. Fifty-four markers have at least four represented taxa. Reuse across identical sequences explains why linked taxon count can exceed the number of taxa in the priority download set. This is an acquisition/mapping checkpoint, not evidence of accelerated structural evolution. The latest receipt and marker counts are copied into `metadata/marker_structure_mapping_snapshot.json` and `metadata/marker_structure_coverage.tsv`.

## Direct matched-position comparisons

```bash
python scripts/compare_marker_structures.py --snapshot results/structural_markers/snapshot-v2 --output results/structural_comparisons/snapshot-v1
python scripts/plot_structural_comparisons.py --comparisons results/structural_comparisons/snapshot-v1
```

These descriptive comparisons use the same homologous matrix positions for sequence differences and structural geometry. Both models must meet each pLDDT threshold (50, 70 or 90) at a position, and both amino acids must be unambiguous. A pair must retain at least 50 qualified residues and at least half of its shared profile positions. Excluded comparisons remain in an audit table.

Global geometry is Cα RMSD after a least-squares proper rotation and translation; mirror reflections are prohibited. A separate local metric measures changes in Cα pair distances for unordered residue pairs within 15 Å in either model and separated by at least three sequence positions in both. Mean absolute and RMS distance changes are reported. This custom local-distance metric is not standard lDDT. Sequence difference is the uncorrected fraction of mismatches at exactly the compared positions. Shared coordinate files are flagged rather than counted as independent experimental confirmations.

The first executed batch contains 2,209 qualifying rows across three thresholds, representing 865 distinct within-marker taxon pairs across 111 markers; 386 threshold-specific comparisons were excluded. At pLDDT ≥70, 858 comparisons qualified. The figure contrasts global RMSD and local geometry but fits no correlation or significance test because observations share proteins, taxa and ancestry. Some large global RMSDs coexist with much smaller local changes; domain motion, uncertain interdomain orientation, alignment error and annotation artifacts remain possible explanations. PAE/domain checks and supported phylogenies are required before interpreting outliers as evolutionary changes.

Results are under `results/structural_comparisons/snapshot-v1`; the displayed SVG is copied to `docs/figures/direct_comparisons.svg`. Numpy, Biopython and plotting dependencies are pinned in `environments/structural-comparisons.yml`. Matrix/mapping/model checksums and rotation/reflection tests validate computation, not the biological assumptions. These distances must not be inserted directly as additive structural branch lengths.

## PAE confidence sensitivity

`retrieve_marker_pae.py` retrieves the version-specific PAE URL recorded for every distinct model in a mapping snapshot. It checks model/version identity, protein length, square matrix dimensions, finite nonnegative entries and the declared maximum (allowing export rounding). Raw JSON is gzip-compressed reproducibly outside Git; receipts preserve compressed and original checksums, sequence identity, retrieval time and source URL. Failed downloads are explicit, never treated as zero uncertainty. A shared lock and hash-checked cache make interrupted acquisition restartable into a new output snapshot.

`assess_pae_sensitivity.py` reuses the exact pLDDT ≥70 residue sets from direct comparisons and independently reproduces their RMSDs and local residue-pair counts. For each unordered residue pair, confidence requires PAE ≤5, ≤10 or ≤15 Å in **both directions in both models**. The analysis reports counts and mean absolute Cα distance changes separately for confident and uncertain pairs, for all nonadjacent residue pairs and for the existing ≤15 Å neighborhood definition. Sequence separation must be at least three residues in each protein. Empty sets yield missing values, not zeros. PAE filters change the residue-pair composition, so changes in averages do not by themselves establish prediction artifacts or biological effects.

PAE describes expected aligned positional error; it is not an experimental error measurement or a calibrated confidence interval on a Cα distance. These descriptive sensitivities do not assign domains, establish domain movement, yield branch rates or resolve prediction circularity. Domain annotation and independent prediction/experimental controls remain required. See the [AlphaFold DB FAQ](https://alphafold.ebi.ac.uk/faq) for metric interpretation and the current two-dimensional JSON format.

```bash
python scripts/retrieve_marker_pae.py --snapshot results/structural_markers/snapshot-v2 --output results/structural_pae/snapshot-v1
python scripts/assess_pae_sensitivity.py --snapshot results/structural_markers/snapshot-v2 --comparisons results/structural_comparisons/snapshot-v1 --pae results/structural_pae/snapshot-v1 --output results/structural_pae/comparisons-v1
```

Use a new output directory for each snapshot. Missing PAE remains an explicit exclusion. The acquisition and analysis scripts use the existing structural-comparisons environment.

`python scripts/plot_pae_sensitivity.py --comparisons results/structural_pae/comparisons-v1` generates SVG/PNG/PDF confidence-sensitivity plots and orders all assessed taxon pairs by global RMSD in `global_rmsd_review_order.tsv`. The ordering is for manual inspection, not a list of statistically supported changes. Species names are joined from the checksum-recorded working manifest. Copies of summary receipts and the SVG are retained in Git; raw PAE and detailed comparisons remain outside Git.

## Expanded mapping and confidence integration

The GDM-only expanded mapping contains 13,153 models, 13,654 marker links and 322 taxa. Its prefetch catalog now has all 13,153 PAE downloads verified with zero failures. The completed acquisition receipt is versioned as `metadata/expanded_marker_pae_prefetch_receipt.json`. Native coordinate encoding is already complete, but coordinate confidence alone does not meet the joint structural-feature filter.

`scripts/advance_expanded_structural_confidence.py` records the original prefetch process identity (PID, start time and command) and waits for actual termination. It requires complete verified acquisition, checks catalog/final-mapping agreement, revalidates cached PAE into a receipt bound to the final mapping, then runs `qualify_native_pae.py`. It checks qualified output hashes, the full model universe, valid/invalid PAE entries and every per-model joint pLDDT ≥70 / PAE ≤10 count. Producer hashes and immutable output paths are pinned in `metadata/expanded_confidence_pipeline_config.json`.

The current controller is running under `results/structural_alphabet/expanded-confidence-control-v1`; its stage logs are there and its main log is `logs/expanded_confidence_pipeline_v1.log`. Mapping-bound PAE goes to `results/structural_pae/gdm-expanded-v1`, and qualified encodings to `results/structural_alphabet/audited-gdm-expanded-v1`. These downstream stages are not yet recorded as complete. A failed or interrupted partial child stage requires inspection and a new output location; do not overwrite partial outputs or substitute the catalog receipt for final residue-mapping provenance. Existing completed PAE cache files remain reusable.
