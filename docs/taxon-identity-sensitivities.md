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
