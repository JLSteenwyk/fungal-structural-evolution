# Fungal structural evolution

Comparative structural genomics of approximately **500 fungal species plus 25 non-fungal outgroups**. The central question is where structural evolution accelerates or decouples from sequence evolution, and how these changes relate to duplication and ecology.

Status: project initialized; public assembly catalogs and taxonomic coverage inventoried; species selection underway. A provisional draft of 502 fungi and 25 genome-backed outgroups is under QC; 499 NCBI proteomes plus two published gap-fill proteomes have been acquired. Final taxa, structural predictions, and evolutionary results are not yet established. Full-scale design; no separate pilot.

## Project records
- [Original objective](docs/objective.txt)
- [Research plan](docs/research-plan.md)
- [Progress and completion evidence](docs/progress.md)
- [Decisions and questions](docs/decisions.md)
- [Literature](docs/bibliography.md)
- [Resource assessment](docs/resources.md)
- [Metadata definitions](metadata/README.md)

## Reproduction
Python 3.10+ standard library suffices for NCBI metadata inventory; openpyxl is used to import published supplementary tables. Run `make inventory` to retrieve public NCBI fungal catalogs, record checksums, and generate candidate tables. These are discovery catalogs, not a final sample. Raw downloads remain in ignored `data/`; versioned metadata records source URLs and hashes.

Downstream workflow implementation is pending sample selection and data inventory. Large files must not enter Git history. GitHub: public repository https://github.com/JLSteenwyk/fungal-structural-evolution (user-confirmed).
