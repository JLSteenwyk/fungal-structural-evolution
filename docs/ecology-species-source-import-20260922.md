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
