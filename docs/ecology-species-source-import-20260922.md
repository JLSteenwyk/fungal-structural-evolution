# Species-level supplementary ecology source

Supplementary Dataset 2 of [Miyauchi et al. 2020](https://doi.org/10.1038/s41467-020-18795-w)
contains 135 genome entries. Matching its species names exactly to the fungal
sampling manifest, with only surrounding whitespace removed, links 36 entries.
Fifteen already have records in the 26-species curated table; 21 provide new
species-level source links. These imports remain separate pending identity
and source-conflict review; the curated table still contains 26 species.

The matched source labels comprise 14 ectomycorrhizal, eight saprotrophic,
five wood-decayer, five pathogen, and one each orchid-mycorrhizal, endophyte,
yeast and arbuscular-mycorrhizal classifications. Labels are preserved
literally in the machine-readable table. They are not a mutually exclusive
trait ontology: for example, yeast describes growth form. No automatic
binary ecological coding is applied.

The source classifies Ramaria rubella and Amanita thiersii as saprotrophic,
disagreeing with their genus-derived ectomycorrhizal candidates. The existing
curation calls A. thiersii asymbiotic based on another source; those labels
need not be contradictory. Ramaria's source strain field additionally names
R. acris, requiring taxonomic review. The table's Punctularia strigosozonata
pathogen label is retained as reported and requires review before biological
interpretation. Exact species-name matching does not establish selected-strain
or assembly equivalence, and the imported rows explicitly leave that pending.

`metadata/miyauchi_species_ecology_matches.tsv` retains the source strain
description, JGI identifier, worksheet row, DOI, download URL, source checksum,
selected assembly, earlier classification and identity flags. A separate
135-row disposition table preserves unmatched and ambiguous records; synonyms
are not silently resolved. Both are produced by
`python scripts/import_miyauchi_species_ecology.py`.

Source workbooks are outside Git under `data/traits/miyauchi-2020/`.
`metadata/miyauchi_ecology_source_downloads.json` records their retrieval URLs,
sizes and hashes; download those URLs to the recorded paths before running
the importer on another machine. Dataset 1 was also retrieved for subsequent
strain/accession review but does not supply classifications to this import.
The importer verifies the downloaded bytes before reading Dataset 2 with
openpyxl 3.1.5. Its read-only parser reports an unsupported Excel extension;
the file is never rewritten. An independent standard-library ZIP/XML
readback verifies all 135 source name/ecology cells, all 36 matched strain
and identifier cells, the complete name intersection and assembly joins.

These additional source links support expanded trait curation, not completed
ecological tests or a count of independent transitions. In particular, the
earlier six-change diagnostic still uses its frozen 26-species table.

## Published-sample linkage

`python scripts/link_ecology_source_biosamples.py` joins Dataset 1's explicit
identifiers to the exact selected assembly accessions in the cached NCBI
assembly catalogue. It also checks catalogue BioSamples against the
independently parsed assembly-statistics table (treating catalogue `na` and
an empty report field as missing). Eight species match both the published
BioSample and WGS accession root:

- Amanita rubescens
- Cantharellus anzutake
- Gautieria morchelliformis
- Hydnum rufescens
- Hysterangium stoloniferum
- Ramaria rubella
- Thelephora ganbajun
- Thelephora terrestris

Cantharellus has a different selected BioProject identifier, PRJNA691513,
from the paper's PRJNA245611, while its BioSample SAMN02745813 and WGS
project WJDU00000000 agree. This discrepancy is retained. WGS version
suffixes are explicitly ignored for project linkage; Thelephora terrestris
has selected WGS version `.2`. No claim of identical assembly sequences or
annotation versions follows. Amanita's paper strain spelling includes an
accent; identifiers, rather than a guessed transliteration, establish the
sample linkage.

The 36-row `metadata/miyauchi_ecology_sample_identity.tsv` preserves all
source and selected identifiers. Twenty-eight rows have no unique published
sample identifier in Dataset 1 and remain unresolved, rather than being
declared mismatches. Raw ZIP/XML parsing independently checked the eight
matched source rows and identifier comparisons. Receipts pin the workbook,
catalogue, statistics table, imported species table and implementation.
These sample links do not resolve source ecological-label disagreements,
taxonomic concepts or selected-isolate experimental phenotype validation.
