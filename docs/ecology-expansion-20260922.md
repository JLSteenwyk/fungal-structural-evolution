# Additional species evidence and structural coverage

The expanded evidence table has 26 species statements, retaining the previous
21 unchanged and adding five classifications reviewed in [Miyauchi et al.
2020](https://doi.org/10.1038/s41467-020-18795-w), Results, “Losses of PCWDEs”.
The Cantharellales paragraph identifies Cantharellus anzutake and Hydnum
rufescens as ectomycorrhizal, Botryobasidium botryosum as saprotrophic, and
Tulasnella calospora as an orchid symbiont. The next paragraph identifies
Sphaerobolus stellatus as saprotrophic. These statements broaden the curated
comparison set beyond the existing Amanita and Suillus sampling.

| Species | Published state | ESMFold only | AlphaFold only | Recovered sequences without models |
|---|---|---:|---:|---:|
| Cantharellus anzutake | Ectomycorrhizal | 97 | 0 | 9 |
| Hydnum rufescens | Ectomycorrhizal | 0 | 116 | 2 |
| Botryobasidium botryosum | Saprotrophic | 0 | 121 | 2 |
| Tulasnella calospora | Orchid mycorrhizal | 0 | 104 | 1 |
| Sphaerobolus stellatus | Saprotrophic | 42 | 57 | 2 |

Counts are marker/protein links from the current frozen availability table,
before confidence qualification. None of these five has a link in both source
catalogs. Unrecovered sequences are separately recorded in the output table;
the displayed counts do not imply 125 recovered markers per taxon.

Cantharellus comparisons currently confound taxon identity with prediction
source. Hydnum and the two Cantharellales comparators have AlphaFold coverage,
but taxon totals do not establish common-marker coverage or qualified paired
sites. Orchid symbiosis remains a distinct state; it must not silently become
an absence-of-symbiosis control. The species classifications do not validate
the selected isolates, assign transitions to branches, establish their
independence, or complete the ecological association tests. Rooting, trait
mapping, matched-source coverage and phylogenetic uncertainty remain gates.

Reproduce with `python scripts/extend_ecology_species_evidence.py`.
Inputs and outputs are hashed in
`metadata/ecology_evidence_expansion_20260922_receipt.json`. The added source
statements are in `config/ecology_species_evidence_additions_20260922.json`;
the combined table is `metadata/species_ecology_evidence_expanded_20260922.tsv`.
The script checks exact species/assembly joins and preserves the earlier
evidence table used by frozen analyses. The new table has not yet been used
to fit ecological effects. Source passages were reviewed online; no local
article checksum or experimental verification is claimed.
