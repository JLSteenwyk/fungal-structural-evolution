# Decisions and unresolved questions

- User-directed scope: approximately 500 unique fungal species and 25 additional non-fungal outgroups, no separate pilot.
- Repository name: fungal-structural-evolution. User confirmed public repository JLSteenwyk/fungal-structural-evolution; origin configured and remote verified empty before initial push.
- NCBI RefSeq and GenBank fungal assembly catalogs are initial discovery sources, not sufficient coverage by themselves. Augment sparse lineages with published datasets and repositories.
- Fungal circumscription is unresolved. Szánthó et al. (2025) use a narrower osmotrophic Fungi definition, treating Aphelida, Rozellida and Microsporidia as close relatives. Record NCBI classification separately from study role; do not silently count boundary taxa as non-fungal outgroups. Assess alternative definitions before freezing roles.
- Shared host: do not infer exclusive resource allocation from physical capacity. Use conservative metadata jobs initially; no GPU jobs launched.
- Existing `gh` executable is not the GitHub CLI (it rejects `auth status`); a GitHub configuration file exists. Do not expose credentials. Establish an appropriate client before remote creation.

Open: outgroup availability and boundary definition; trait sources/replicate transitions; predictor/database availability; actual reusable structure coverage; sustained shared-machine allocation.
