# Species-level ecological source review

Reviewed all 15 species imported from Miyauchi Dataset 2 that were absent from
the 32-statement curated table. A new version preserves those 32 statements
unchanged and adds 13 published ecological classifications, for **45 species**.
This records what the source reports; it does not validate exclusive lifestyles,
the phenotype of selected isolates, or independent transitions.

| Published class retained | Added species |
|---|---|
| Wood decayer | Neolentinus lepideus; Jaapia argillacea; Auricularia subglabra; Dacryopinax primogenitus |
| Saprotroph | Amorphotheca resinae; Aureobasidium pullulans; Morchella importuna |
| Pathogen | Rhizoctonia solani; Botrytis cinerea; Aspergillus flavus; Taphrina deformans |
| Endophyte | Xylona heveae |
| Arbuscular mycorrhizae | Rhizophagus irregularis |

These are species classifications in [Miyauchi et al., Supplementary Dataset 2](https://www.nature.com/articles/s41467-020-18795-w).
Exact species names and selected assembly accessions come from the previously
independently verified source import. Dataset 1 provides no unique sample
identifier for these additions; selected-isolate equivalence remains unresolved.
Wood decay, pathogenicity and endophytism are not silently mapped to absence of
ectomycorrhizal capability. New rows have no inferred transition group; their
`unassigned_<taxon>` labels are bookkeeping identifiers, not replicate origins.

Two source records remain separate:

- **Saccharomyces cerevisiae:** the imported yeast label describes growth form.
  It is retained in the review disposition table, not added as an ecological
  guild. The source's M3837 strain is not equated to the selected assembly.
- **Punctularia strigosozonata:** the source's pathogen label is preserved but
  withheld from ecological coding. A search result for the primary
  [JGI genome project](https://myco-lb.jgi.doe.gov/Punst1/Punst1.home.html)
  reports white rot; direct page retrieval failed during this review. This is
  a follow-up lead, not a fully adjudicated correction, and neither label is
  substituted into the curated table.

The distinction between reported species labels and selected-isolate evidence
is retained for every new record. Rhizoctonia source AG-1 IB, Aureobasidium source
var. pullulans EXF-150 and other source strain descriptions remain available in
the immutable import; no strain or taxonomic-concept equivalence is inferred.

## Reproduction and validation

```bash
python scripts/extend_reviewed_species_ecology.py
python scripts/check_reviewed_species_ecology.py
```

The explicit decisions and pinned inputs are in
`config/ecology_species_source_review_20260927.json`. The new table is
`metadata/species_ecology_evidence_reviewed_20260927.tsv`; all 15 dispositions,
the producer receipt and source-join readback are under
`metadata/ecology_species_source_review_*_20260927.*`.

The checker verifies all 32 preserved rows, 13 added identities and source
fields, retained states, unresolved-status fields and all 15 review decisions.
It relies on the existing full raw-XLSX XML verification of the imported source.
This verifies transcription and linkage, not the truth of ecological labels.
Previously computed parsimony and coverage results retain their frozen 32-taxon
inputs. Extending those diagnostics and reviewing independent transitions are
next steps; no new ecological association result is claimed.

## Primary wood-decay evidence added

[Figure 1 of Riley et al. (2014)](https://doi.org/10.1073/pnas.1400592111)
was visually reviewed and transcribed as a separate trait table: 12 white-rot,
seven brown-rot and three uncertain classifications. Seven source names match
selected ingroup names exactly: four white-rot, one brown-rot and two uncertain.
All 15 unmatched names remain visible; no synonym transfer was inferred.
In particular, Jaapia argillacea and Botryobasidium botryosum retain the figure's
uncertain decay-mode category. The study cautions against a simple decay-mode
dichotomy. These states are not binary ECM assignments or validated origins.

The reviewed source PDF is stored outside Git at
`data/traits/riley2014-primary.pdf`, with its URL/hash and 22 source statements in
`config/wood_decay_primary_species_20260927.json`. Reproduce the exact-name join
with `scripts/import_primary_wood_decay_species.py` using fresh output paths.
Every exported identity and accession was checked against a separate scalar
manifest lookup. Outputs are `metadata/wood_decay_primary_species_20260927.tsv`
and its adjacent receipt. An initial import used the wrong manifest field name;
it was corrected to `study_role` before any table was written.

### Punctularia source discrepancy

An independent [primary culture study](https://doi.org/10.1371/journal.pone.0130381)
also describes Punctularia strigosozonata as white rot and identifies strain
HHB-11173 SS-5. Deposited assembly metadata use HHB-11173 SS5; these strain
strings agree after punctuation normalization, but complete experimental-isolate
provenance is not established. The primary classification is now recorded
alongside the earlier supplementary pathogen label. These roles are not assumed
mutually exclusive, so the new evidence does not justify declaring the earlier
label erroneous. Neither source changes the frozen ecological reconstructions.

The publisher HTML snapshot is stored outside Git, with source location, hash,
strain fields and review decision in
`metadata/punctularia_primary_ecology_review_20260927.json`. This replaces the
previous reliance on a search-result-only lead with directly inspected primary
literature. It does not establish absence of pathogenicity or validate a binary
trait coding for ecological effect tests.

## Qualified structural coverage for the seven exact-name matches

Both source-specific coverage analyses and full independent matrix readbacks
have completed. All seven taxa have eligible AlphaFold markers, spanning 19 of
21 taxon pairs with at least one marker sharing 50 qualified columns. ESMFold
covers three taxa and three pairs. Its three covered taxa are all classified
as white rot in the reviewed figure. The sole matched brown-rot taxon,
Gloeophyllum trabeum, has AlphaFold coverage but no eligible ESMFold markers.
Thus ESMFold alone cannot currently supply a between-class comparison here.

AlphaFold marker counts range from two (Fomitiporia mediterranea) to 122
(Gloeophyllum trabeum); all four white-rot taxa have no marker in common as a
complete group. Pairwise availability does not establish balanced sampling or
independent ecological replication. A single brown-rot tip is insufficient to
claim independently replicated brown-rot origins. The two uncertain taxa remain
uncertain and are not forced into either class.

Plans and proofs are `metadata/qualified_{afdb,esmfold}_wood_decay_*_20260927.*`;
completion is recorded in `metadata/wood_decay_coverage_completed_20260927.json`.
Full tables remain under `results/ecology/qualified-{afdb,esmfold}-wood-decay-overlap-20260927-v1`.
Reproduction uses `scripts/prepare_wood_decay_coverage.py`, then the existing
source-specific `assess_qualified_*ecology_overlap.py` and
`readback_qualified_ecology_overlap.py` with these plans and fresh output paths.
All 21 pairs, per-marker counts and three descriptive source-category groups
were checked independently. Category grouping is bookkeeping for coverage,
not a grouping of inferred evolutionary origins. Each stage had a one-CPU,
4-GiB, 1–10-minute planning allowance and completed in seconds; no GPU was used.

## Targeted name review

A three-taxon review is recorded in
`metadata/wood_decay_name_review_20260927.json`:

- The [Auricularia BioProject](https://www.ncbi.nlm.nih.gov/bioproject/60553)
  binds selected assembly GCA_000265015.1 and strain TFB-10046 SS5 to the old
  delicata label and current subglabra label. The precise Riley study-sample
  link remains pending; no species-wide synonym transfer was made.
- The [Dacryopinax genome project](https://mycocosm.jgi.doe.gov/Dacsp1/Dacsp1.home.html)
  directly describes named D. primogenitus as brown rot. This supplies independent
  species evidence without equating every unnamed Dacryopinax to that species.
- The [Heterobasidion strain record](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=747525)
  places TC 32-1 in H. irregulare within the H. annosum species complex. Complex
  membership alone is insufficient to transfer the older source classification.

All three selected names and assemblies were checked against the manifest.
The two Auricularia pages have hashed local snapshots. JGI's full description
was inspected in the browser, but its separate local download returned HTTP
403; no local snapshot is claimed. Existing classification and coverage tables
remain frozen. Incorporating the new independent Dacryopinax statement into a
versioned expanded analysis is the next step.

## Eight-taxon coverage update

A versioned expansion now adds the independently sourced Dacryopinax primogenitus
brown-rot statement. The original seven evidence rows are preserved exactly.
The new species has 116 eligible AlphaFold markers and no eligible ESMFold
markers. Across all eight taxa, AlphaFold supplies at least one 50-column marker
for 26 of 28 pairs; ESMFold still supplies three pairs, all among white-rot taxa.
The table now has four white-rot, two brown-rot and two uncertain classifications.
Two brown-rot labels are not themselves evidence of two independent origins.

Both full matrix readbacks passed. A separate comparison verified exact
preservation of all original taxon, pair and pair-marker rows. Category-level
summaries were recalculated because the brown-rot category gained a member.
The frozen seven-taxon outputs remain intact. Preparation is reproducible with
`scripts/extend_wood_decay_coverage.py`; use the
`metadata/qualified_{afdb,esmfold}_wood_decay_v2_plan_20260927.json` plans for
coverage and readback, with fresh paths when rerunning. Completion evidence is
`metadata/wood_decay_coverage_v2_completed_20260927.json`.
Full outputs are under `results/ecology/qualified-{afdb,esmfold}-wood-decay-overlap-20260927-v2`.
As before, the one-CPU, 4-GiB, 1–10-minute per-stage allowance used no GPU.

The [complete edge-placement diagnostic](wood-decay-phylogenetic-diagnostic-20260927.md)
now maps these labels onto both completed tree alternatives, including all four
binary sensitivity assignments of the two uncertain taxa. Full independent
network-flow verification passed. The primary coding requires two changes but
no particular edge; a required edge appears only under a forced Jaapia coding.
This does not establish independently replicated ecological origins.
