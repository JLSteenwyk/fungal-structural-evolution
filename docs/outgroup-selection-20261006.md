# Outgroup selection record — 2026-10-06

This record makes the frozen 25-taxon non-fungal panel auditable. It is a
sampling rationale, not a new phylogenetic result. The selected taxa remain in
the 526-entry analysis manifest; protein-family analyses may use a suitable
subset when missingness, homology, or alignment coverage requires it.

## Evidence and design

The close anchor is *Fonticula alba* (Nucleariida). Brown et al. recovered
*Fonticula* plus *Nuclearia* as sister to Fungi, and explicitly named the
resulting clade Nucletmycea ([Mol. Biol. Evol. 2009](https://doi.org/10.1093/molbev/msp185)).
This makes it a useful close non-fungal comparison, while a single genome is
not treated as proof of a root.

The panel also retains broad sampling across unicellular holozoan and animal
lineages. Grau-Bové et al. analyzed genomes of animals and their closest
unicellular holozoan relatives, including three ichthyosporeans and
*Corallochytrium* ([eLife 2017](https://doi.org/10.7554/eLife.26036)). Their
phylogenomic analysis provides a direct precedent for retaining those
comparators. A broader fungal timetree used 110 fungal and 43 non-fungal taxa
and retained markers with coverage in both fungal groups and outgroups
([Szánthó et al. 2025](https://doi.org/10.1038/s41559-025-02851-z)); it informs
the broad-coverage and marker-occupancy safeguards here rather than supplying
this study's topology or dates.

## Frozen membership and roles

| Selection role | Lineages and count | Taxa |
| --- | --- | --- |
| Close holomycotan anchor | Nucleariida: 1 | *Fonticula alba* |
| Unicellular holozoan comparisons | Choanoflagellatea: 2; Filasterea: 1; Ichthyosporea: 5; Corallochytrea: 1 (9 total) | *Monosiga brevicollis*, *Salpingoeca rosetta*, *Capsaspora owczarzaki*, *Sphaeroforma arctica*, *Chromosphaera perkinsii*, *Pirum gemmata*, *Abeoforma whisleri*, *Creolimax fragrantissima*, *Corallochytrium limacisporum* |
| Animal lineage coverage | Metazoa: 9 | *Amphimedon queenslandica*, *Nematostella vectensis*, *Trichoplax adhaerens*, *Mnemiopsis leidyi*, *Drosophila melanogaster*, *Caenorhabditis elegans*, *Capitella teleta*, *Lottia gigantea*, *Branchiostoma floridae* |
| Broader rooting sensitivity | Apusomonadida: 1; Amoebozoa: 5 (6 total) | *Thecamonas trahens*, *Acanthamoeba castellanii*, *Dictyostelium discoideum*, *D. purpureum*, *Polysphondylium violaceum*, *Planoprotostelium fungivorum* |

The exact members, assemblies, source versions, inclusion reasons and
provisional statuses are in the [analysis manifest](../metadata/analysis_manifest.tsv)
and [outgroup draft](../metadata/outgroup_sampling_draft.tsv). The
[selection audit](../metadata/outgroup_selection_audit_20261006_v1.json)
checks that this table remains exactly concordant with
[`config/outgroups.json`](../config/outgroups.json), including the cited DOI
links. It reports one close holomycotan anchor, nine unicellular holozoan
comparators, nine metazoans, and six broader rooting-sensitivity taxa.

## Interpretation limits and next gate

The configured labels are practical sampling roles. They do not demonstrate a
species-tree topology, locate a root along an edge, establish branch times, or
show that an assembly is fit for all gene families. The direct root-split
diagnostic only checks the manifest role partition and is explicitly not a
rooting result. The species-tree workflow will assess marker, taxon and model
sensitivities before accepting a topology; family analyses will record their
actual retained outgroups rather than force all 25 into every analysis.
