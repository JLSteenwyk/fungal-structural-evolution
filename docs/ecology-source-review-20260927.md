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
