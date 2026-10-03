# Full selected-taxon identity audit

All **526 working entries (501 fungal entries and 25 outgroups)** now have an
assembly-linked audit against the preserved NCBI catalogue and taxonomy snapshot.
All effective IDs are species-rank records with complete paths to the root.
Independent names/ranked-lineage lookup agrees with node traversal for every
entry. No species-rank ID is repeated, including across study roles. No exact
selected assembly accession is repeated or paired with another selected accession.

These checks do **not verify 501 biological fungal species**. Twenty-one
uncertain fungal labels, two curated fungal hybrids, the three-entry
Neocallimastix species complex and four prior strain-linked cautions remain
explicit. The hybrid screen is not comprehensive. Database ranks and accession
uniqueness cannot resolve conspecificity, contamination or homeolog ancestry.

The audit joins 521 exact accessions to frozen raw assembly catalogues. Five
outgroup Figshare depositions retain separate provenance without fabricated
NCBI assembly records. The existing verified Chromosphaera overlay supplies the
one missing manifest identifier. Snapshot hashes are checked; no current-taxonomy
claim is made. Two Sanchytriomycota additions, F2109901/GCA_018985225.1 and
F2020955/GCA_018985235.1, are absent from the earlier species candidate inventory
but present in the complete catalogue. The initial audit's incorrect membership
assertion and source are preserved. Version 2 checks the full inventory and
reports earlier membership separately, retaining both entries.

## Amphioxus finding and interpretation

**O2700040** has the frozen label *Branchiostoma floridae*, but its selected
assembly **GCF_019207075.1** is paired identically with GCA_019207075.1, Bb-hap.
The [NCBI submitter project](https://www.ncbi.nlm.nih.gov/bioproject/706157)
assigns this principal haplotype to **B. belcheri**, derived from a
B. floridae × B. belcheri specimen. BioSample SAMN13907882, both assembly
reports and original cached statistics agree. The BioSample project
PRJNA602496, assembly project PRJNA706157 and RefSeq annotation project
PRJNA1402437 remain distinct linked records.

The [primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10013865/)
describes partitioning parental reads to assemble species-specific haploid
genomes. Hybrid specimen origin therefore does not establish mixed parental
sequence in the selected assembly. The submitter assignment is verified as
deposition evidence; independent genome-wide/locus-specific purity remains open.

All **125 marker slots** are retained: 100 selected, 25 not selected. Every
selected protein matches its cached sequence hash and maps through all annotated
CDS segments to deposited sequence accessions. Full GFF coordinates are checked
against all 79 assembly sequences. The annotation and normalized proteome each
contain 48,498 protein IDs; this is not retranslation or orthology validation.

The [versioned overlay](../config/assembly_identity_reviews_20261003_v2.json)
supplies “Branchiostoma belcheri (hybrid-derived principal haplotype)” as a display
name, retaining the stable ID and deposition TaxID 2700040. Consumers must adopt
it explicitly. All 25 outgroups, frozen manifests, query IDs, sequences and live
jobs remain unchanged. Species-level interpretation and trait joins must account
for this identity; relevant exclusion/reference sensitivities remain downstream.

## Artifacts and field meanings

Tables and primary archives remain outside Git; receipts supply exact paths/hashes.

| Artifact | Meaning |
| --- | --- |
| `results/taxonomy/selected-identity-snapshot-20261003-v2/taxon_identity.tsv` | All 526 frozen labels, effective/canonical IDs, node paths, names, exact assembly/strain/sample/project/publication links and overlapping review flags |
| `source_crosschecks.tsv` in the same directory | Root-path, independently ranked-name and raw-catalogue checks for every entry |
| `species_rank_collisions.json` in the same directory | Collision, repeated-accession and selected-paired-accession lists, including empty lists |
| `results/taxonomy/branchiostoma-assembly-identity-20261003-v2/all_marker_provenance.tsv` | All 125 marker slots, selected protein hashes and complete CDS-segment coordinates |
| `assembly_sequences.tsv` in the same directory | All 79 RefSeq/GenBank sequence pairs, roles and lengths |

`frozen_*` reproduces historical inputs. `effective_species_taxid` uses only
the previously verified identifier override. `snapshot_*` describes the pinned
database, not accepted fungal classification. `assembly_type_representation`
is not biological ploidy. `independent_locus_parental_purity=unverified` means
the coordinate join does not determine ancestry. Unknowns are not imputed.

## Reproduction and resources

Python 3.10.13 and standard-library modules suffice. Use fresh destinations:

```bash
python scripts/audit_selected_taxon_identity_snapshot_v2.py \
  --output NEW_SNAPSHOT_OUTPUT --receipt NEW_SNAPSHOT_RECEIPT.json
python scripts/review_branchiostoma_assembly_identity_v2.py \
  --sources data/taxonomy/branchiostoma-assembly-identity-20261003-v1 \
  --output NEW_IDENTITY_OUTPUT --receipt NEW_IDENTITY_RECEIPT.json \
  --overlay NEW_IDENTITY_OVERLAY.json
```

The second command consumes the first audit's pinned original destinations;
adapt those explicit paths in a new source version for relocated inputs.
Source URLs, retrieval times and hashes accompany the cached XML/reports and
identity receipt. A live taxonomy download may differ from the required snapshot;
do not silently replace it. Failed source versions are preserved, not entrypoints.

Both executions used a single-threaded process, configured 4-GiB address space
and 16-MiB per-file bound. Snapshot CPU/wall limits were 300/600 seconds, with
10–180 seconds planned; actual wall/CPU time was 34.84/34.83 seconds and maximum
RSS 1,277,435,904 bytes. Identity CPU/wall limits were 180/300 seconds, with
1–60 seconds planned; actual wall/CPU time was 9.34/9.32 seconds and maximum RSS
548,929,536 bytes. GNU time reports process RSS, not charged cgroup peak; no
observed CPU quota is claimed. No native inference, GPU work or charges followed.

Evidence: [snapshot completion](../metadata/selected_taxon_identity_snapshot_completed_20261003_v2.json),
[identity/provenance review](../metadata/branchiostoma_assembly_identity_review_20261003_v2.json),
and [execution/resource records](../metadata/taxon_identity_audit_execution_20261003_v2.json).
All table values, artifact/source hashes and overlay binding passed readback.
The [separate serialized evidence check](../metadata/selected_taxon_identity_snapshot_serialized_readback_20261003_v2.json)
verifies 43 bindings, all 526 table rows and all 125 marker slots, including
species-rank identifier distinctness across study roles.
Species boundaries, the accepted phylogeny, reconciliation, dating, calibration
and all eight biological aims remain incomplete.
