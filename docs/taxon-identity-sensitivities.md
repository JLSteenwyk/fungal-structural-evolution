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
