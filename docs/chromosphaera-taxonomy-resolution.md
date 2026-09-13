# Chromosphaera outgroup identifier resolution

The project taxon **OFS5426494**, Chromosphaera perkinsii, now has a verified database query identifier: **NCBI taxon 1932427**, whose database scientific name is **Ichthyosporea sp. XGB-2017a**. Both names are retained rather than silently replacing the publication name.

The [NCBI taxonomy record](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1932427) explicitly includes Chromosphaera perkinsii among its alternate names. The returned alias classification is “unpublished name”; that database field is preserved and is not treated as a conclusion about nomenclatural validity. A direct genus-name search did not discover the record, which is why a name-only lookup had left the identifier unresolved.

Independent deposition evidence comes from [BioSample SAMN06200032](https://www.ebi.ac.uk/ena/browser/view/SAMN06200032) in BioProject PRJNA360047: the sample alias and title identify Chromosphaera perkinsii genome sequencing, and runs SRR5195170/SRR5195174 use taxon 1932427. The [original genome study](https://doi.org/10.7554/eLife.26036) links its genome data to that BioProject. The NCBI taxonomy XML, ENA study-run table and public article XML are archived with URLs and SHA256 hashes in `metadata/chromosphaera_taxonomy_resolution.json`.

`config/taxonomy_overrides.json` supplies the identifier for current inventory queries. `match_uniprot_structures.py` verifies the project taxon/name, evidence receipt, source-file checksums and any existing identifier before applying it. Historical sampling manifests and inputs of active analyses remain unchanged. A conflicting name or identifier fails validation; a cached UniProt query for another identifier cannot be silently reused. A focused test covers accepted linkage, conflicting identity and altered evidence.

The completed UniProt query for taxon 1932427 returned **zero UniProt records**, including zero AlphaFoldDB cross-references, against the project's 12,463-protein input. This closes the missing-query state; it is not evidence that no structure exists in any database. All 526 working taxa now have completed query/exact-matching records. The completed inventory is frozen at `results/structure_inventory/completed-query-snapshot-v2` and summarized in `metadata/uniprot_inventory_completion.json`.

Reproduce the evidence retrieval with `python scripts/resolve_chromosphaera_taxonomy.py` in a fresh evidence directory/version, then run `python scripts/match_uniprot_structures.py`. Existing pinned evidence is deliberately not overwritten. The latter command resumes only unresolved or identifier-changed query states. Tests: `python -m unittest discover -s tests -p test_taxonomy_overrides.py`.

## Complete-inventory prediction candidates

A new frozen candidate inventory retains all 58,883 unique marker sequences and 59,840 taxon/marker links. It classifies 14,516 sequences as having a verified reuse receipt, 16,106 as nominated retrievals still pending, and 28,261 as having no candidate in the completed inventory. No sequence remains in the earlier “inventory pending” category. These reuse-receipt labels do not replace model-level coordinate/source validation.

Compared with the initial 20,480-candidate queue, 7,799 candidates are additional, 20,462 remain candidates, and 18 now have reuse candidates or receipts. Of the additional sequences, 4,252 are canonical and within the existing ≤512-residue execution configuration. The delta excludes every sequence in the previous queue, including its length-deferred sequences; it neither resubmits them nor changes the active GPU job. Longer/noncanonical proteins and validation/reuse of existing local predictions remain separate work before future execution. This preparation does not launch the new candidate set.

```bash
python scripts/prepare_marker_prediction_inputs.py \
  --output data/prediction_inputs/markers-full-inventory-v2
python scripts/compare_prediction_queues.py \
  --previous data/prediction_inputs/markers-v1 \
  --current data/prediction_inputs/markers-full-inventory-v2 \
  --output data/prediction_inputs/full-inventory-delta-v2
```

Use new output versions if these immutable snapshots already exist. Receipts are tracked as `metadata/full_inventory_prediction_input_receipt.json` and `metadata/full_inventory_prediction_delta_receipt.json`; candidate FASTAs and detailed per-sequence records remain outside Git. The full fungal/outgroup sampling and all evolutionary objectives remain active and incomplete.
