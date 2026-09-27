# Taxon identity and hybrid sensitivity inputs

The working dataset contains 501 fungal entries and 25 outgroups. Distinct labels
or taxonomy IDs do not establish distinct biological species. The original
name-only screen flagged 21 incompletely identified labels and one explicit
hybrid; it could not detect hybrid ancestry hidden behind an ordinary binomial.

Two selected assemblies are now linked to primary hybrid genome studies through
exact accession, strain and publication ID in the frozen NCBI assembly catalogue:

- VIN7, GCA_000326105.1, PMID 22136070: [Borneman et al.](https://doi.org/10.1111/j.1567-1364.2011.00773.x) report an allotriploid S. cerevisiae–S. kudriavzevii genome.
- CBS 1483, GCA_011022315.1, PMID 31791228: [Salazar et al.](https://doi.org/10.1186/s12864-019-6263-3) characterize the S. pastorianus hybrid genome, with S. cerevisiae and S. eubayanus ancestry.

Both catalogue rows contain `assembly_type=haploid`. That assembly representation
field must not be used as a biological ploidy measurement or proof of nonhybrid
ancestry. The literature evidence is curated in `config/curated_hybrid_evidence.json`;
parental assignment of individual sampled proteins remains unresolved. Hybrid
homeologs could otherwise resemble ordinary gene duplications or discordant
species-tree markers. These two records are not a comprehensive hybrid screen.

## Prepared exclusions

The builder produces explicit membership tables and exact sequence subsets for
both the 49,027-column profile alignment and 63,750-column MAFFT alignment:

| Sensitivity | Fungal entries | Outgroups | Total |
| --- | ---: | ---: | ---: |
| Exclude two curated hybrids | 499 | 25 | 524 |
| Also exclude 21 uncertain labels | 478 | 25 | 503 |

These are sensitivity inputs, not a replacement sampling target or a validated
unique-species count. All original taxa, sequences and frozen analyses remain
available. Every retained sequence is unchanged, including alignment columns;
FASTA readback and complete taxon-set checks pass. All 25 outgroups remain in
both policies. The original name-only audit remains frozen and separate from
this literature-based evidence layer; consumers of that older audit do not
implicitly acquire the additional S. pastorianus flag.

```bash
python scripts/prepare_taxon_identity_sensitivities.py \
  --output results/phylogeny/taxon-identity-inputs-v1
```

Hash-pinned source receipts, curated evidence, membership tables and matrix
checksums are versioned. Large matrices remain outside Git. Trees have not yet
been inferred from these subsets, so no topology robustness is claimed. The
21 uncertain labels still need assembly-linked taxonomic work; subgenome-aware
analyses of the hybrids remain a separate requirement before ordinary
bifurcating reconciliation or selection analyses can include their copies.


## Guide sensitivity execution

All four prepared matrices are now submitted to `run_taxon_sensitivity_trees.py`:
two concurrent trees, eight threads and 32 GB maximum memory each, LG+F+G4 and
seed 20260913. These match the full guides' model/seed; thread counts differ.
The resource envelope is 24–144 hours for the batch, based on observed ongoing
full-guide searches, with 5 GB output headroom on the existing host. No extra
support analysis is included here. IQ-TREE checkpoints support interrupted-run
continuation with an unchanged configuration; completed artifact hashes are
checked before reuse. Final tree validation requires exact expected tips and
finite nonnegative branch lengths.

The lineage-retention table shows no entire role/major-lineage group is removed
by either policy. This does not imply unchanged within-lineage coverage or
adequate marker information. The four guide fits are running, so their topology
comparisons and robustness conclusions remain pending.

```bash
python scripts/run_taxon_sensitivity_trees.py \
  --inputs results/phylogeny/taxon-identity-inputs-v1 \
  --output results/phylogeny/taxon-identity-trees-v1
```


## Source-specific paired marker coverage

Reapplying each exclusion policy to the frozen paired inputs, while retaining
per-taxon confidence/coverage criteria and rechecking the four-taxon family gate,
gives the following availability counts:

| Source and policy | Eligible markers | Taxa with an eligible marker | Taxon-marker cells |
| --- | ---: | ---: | ---: |
| AlphaFold, full | 124 | 322 | 13,565 |
| AlphaFold, hybrids excluded | 124 | 320 | 13,510 |
| AlphaFold, hybrids and uncertain labels excluded | 124 | 320 | 13,510 |
| ESMFold, full | 72 | 199 | 4,669 |
| ESMFold, hybrids excluded | 72 | 199 | 4,669 |
| ESMFold, hybrids and uncertain labels excluded | 71 | 184 | 4,286 |

The stricter ESMFold policy removes marker 776280at2759 from eligibility: its
four original taxa are Tilletia caries, Ceratobasidium sp. AG-Ba, Zalaria obscura
and Trypethelium subeluteriae. Excluding the uncertain Ceratobasidium label leaves
three. Neither hybrid has eligible coverage in this local snapshot. Thus taxon
policy changes the available source-specific data unevenly; topology or rate
sensitivity cannot be interpreted as a pure taxonomic effect without examining
these coverage changes. No new paired fits or effect estimates are claimed.

```bash
python scripts/assess_identity_policy_marker_coverage.py \
  --output results/phylogeny/taxon-identity-paired-coverage-v1
```

Source receipts and actual ready-marker FASTA taxon identities were checked.
Outputs include marker-level and role/lineage-level counts, preserving groups
with zero coverage. The full-panel baseline reproduces the earlier source
coverage counts. Exclusion policies cannot make a previously ineligible marker
eligible, and no new confidence values are imputed.

## September 27: Neocallimastix species-complex sensitivity

A primary phylogenetic study supports potential conspecificity of
*N. californiae*, *N. lanati* and *N. cameroonii*, while emphasizing the need for
stronger phylogenetic resolution ([Stabel et al., 2021](https://doi.org/10.3390/microorganisms9081655)).
The genome study for the selected `Neocon1` assembly identifies its isolate as
*N. cameroonii* var. *constans* and places it in the cameroonii/californiae
species complex ([genome study, 2025](https://doi.org/10.1093/g3journal/jkaf137)).
Its BioProject PRJNA1052201 and WGS project JBODTK000000000 exactly match the
frozen assembly catalogue for GCA_050613775.1. This supplies an assembly-linked
literature name, but does not by itself establish formal synonymy or delimit
all selected genomes as one biological species.

The selected records are:

| Taxon ID | Frozen label | Assembly | Strain |
| --- | --- | --- | --- |
| F1754190 | Neocallimastix californiae | GCA_002104975.1 | G1 |
| F2767002 | Neocallimastix lanati (nom. inval.) | GCA_016946835.1 | sp3 |
| F3108450 | Neocallimastix sp. 'constans' | GCA_050613775.1 | G3 |

`prepare_neocallimastix_identity_sensitivity.py` prepared all three possible
single-representative policies for both profile and MAFFT alignments: six exact
sequence subsets. Every subset contains 499 fungal entries plus all 25 outgroups;
retained sequences and columns are unchanged and fully checked. These are entry
counts, not validated unique-species counts. Other uncertain labels and the two
curated hybrids remain in these particular subsets. Original manifests and
running analyses are unchanged.

The two primary articles are cached as XML with identifiers, source URLs and
hashes. All six matrix hashes, membership and assembly evidence are in
`results/phylogeny/neocallimastix-identity-inputs-20260927-v1`; the versioned
closure is `metadata/neocallimastix_identity_inputs_completed_20260927.json`.
No new sensitivity tree has been inferred yet. Final species counting still
requires broader taxonomic review and genome-wide assessment of this complex.

### Neocallimastix guide searches launched

All six prepared matrices are now submitted to
`run_neocallimastix_sensitivity_trees.py`. Two searches run concurrently with
eight threads and 32 GiB maximum memory each, LG+F+G4 and seed 20260913, matching
the prior homogeneous guide model/seed. The service has a 16-CPU quota, 72 GiB
aggregate memory cap and no swap. Planning allowance is 36–336 hours and 5 GiB
output. The host had 192 logical CPUs, load near 63 and about 525 GiB available
memory at launch. These estimates are not a completion guarantee.

The first startup failed before creating tree outputs because its explicit PATH
omitted IQ-TREE's installed directory. That failed unit is retained; a new unit
with the correct executable directory now has two verified running IQ-TREE
children. No inference was restarted or existing output overwritten.

Launch identity and resource estimates are recorded in
`metadata/neocallimastix_guide_launch_20260927.json`. Results will be under
`results/phylogeny/neocallimastix-identity-guides-20260927-v1`. Checkpoint reuse
requires identical configuration and input hashes; successful outputs must have
exact expected tips and finite, nonnegative branch lengths. Full tree readback
and topology comparisons remain downstream. No bootstrap support or species
boundary conclusion is claimed by this stage.

## Neocallimastix tree verification queued (2026-09-27)

`scripts/audit_neocallimastix_guides.py` now waits for the six-tree producer to
reach inactive/success/exit-zero. The waiting service checks the live producer's
PID, creation time and command; pinned scripts and inputs must remain unchanged.
After completion it will verify every saved artifact, exact matrix and tree
identities, 524 tips including all 25 outgroups, matrix lengths, model/seed and
finite nonnegative branch lengths.

All 15 guide pairs will be compared using unrooted nontrivial splits restricted
to common tips. Projected splits are checked against explicit tree pruning.
Comparisons retaining different representatives have 523 shared tips and omit
all three disputed taxa; comparisons retaining the same representative have
524. RF distance is reported with both split counts and RF divided by their
sum, without treating guide differences as branch support or species delimitation.
The pruning check passed known six-tip examples before launch.

The checker uses at most one CPU and 2 GiB RAM, no swap, and a 0.1 GiB output
allowance; planned runtime after tree completion is 0.1–2 hours. No GPU or new
phylogenetic inference is launched. Plan and exact launch identity:
`metadata/neocallimastix_guide_audit_plan_20260927.json` and
`metadata/neocallimastix_guide_audit_launch_20260927.json`.
Full tree verification and topology results remain pending. Comparisons with
full-data guides and downstream structural analyses remain additional work.
