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
