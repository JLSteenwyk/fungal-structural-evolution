# Ancestral reconstruction case inputs

Status: ancestral sequences and structures have not been inferred. This stage
identifies candidate nodes for the 13 whole/domain contrast cases; it does not
replace the requirement to reconstruct ancestral ensembles for selected families.

All 26 resolved gene trees (13 families under each of the profile and MAFFT
species guides) were extracted from the fully checked reconciliation outputs.
For each case we retained the focal duplicate MRCA and its first three ancestral
levels. All 104 neighborhoods, 3,722 descendant records and 52 guide comparisons
passed independent parent-path reconstruction and exact source-tree checks.
All 52 comparisons have identical descendant sets across guides. The two
reconciliations share underlying gene-tree information, so agreement is a
sensitivity result rather than independent support.

The largest neighborhood examined per case is shown below. These are candidate
sampling neighborhoods, not selected final alignment sets or reconstructed nodes.

| Family | Case taxon | Full family proteins | Level-three proteins | Level-three taxa |
| --- | --- | ---: | ---: | ---: |
| OG0000054 | Heliocybe sulcata | 5775 | 11 | 3 |
| OG0000095 | Rhodonia placenta | 4610 | 83 | 24 |
| OG0000107 | Cryoendolithus antarcticus | 3458 | 9 | 8 |
| OG0000152 | Scedosporium apiospermum | 4020 | 7 | 4 |
| OG0000230 | Furculomyces boomerangus | 2215 | 24 | 3 |
| OG0000294 | Jaapia argillacea | 1992 | 9 | 8 |
| OG0000336 | Synchytrium microbalum | 1735 | 10 | 8 |
| OG0000972 | Neolecta irregularis | 1548 | 622 | 389 |
| OG0001082 | Leucosporidium creatinivorum | 940 | 111 | 96 |
| OG0001200 | Smittium simulii | 762 | 29 | 9 |
| OG0001203 | Leucosporidium creatinivorum | 1297 | 21 | 13 |
| OG0002650 | Piloderma croceum | 842 | 77 | 24 |
| OG0002812 | Phycomyces blakesleeanus | 1695 | 12 | 11 |

Node labels remain guide-specific; comparisons use descendant identities. A
matched ancestral level alone would not establish node homology if descendants
differed. All focal MRCAs contain exactly the two nominated duplicate proteins.
No bootstrap support, rooting confidence or branch timing follows from this
membership check.

Next required work is to recover the exact extant sequences, assess full-family
and local alignment quality and domain architecture, and choose supported
reconstruction nodes with appropriate outside sequences. Ancestral inference
must assess model, alignment and topology sensitivity and propagate sitewise
uncertainty into sequence ensembles. Indel uncertainty needs separate treatment.
Selecting only reference-agreeing cases would bias the study; experimental
reference disagreements and missing coverage remain explicit selection factors.
GPU structure prediction remains paused.

Reproduce with `python scripts/prepare_case_ancestral_neighborhoods.py`, then
`python scripts/readback_case_ancestral_neighborhoods.py`. Both preserve existing
outputs. Local source trees, descendant lists and tables are in
`results/ancestral/case-neighborhoods-20260927-v1/`; source and artifact checksums
are in its receipt. Verified closure:
`metadata/case_ancestral_neighborhoods_completed_20260927.json`.

## Exact extant sequences recovered

Recovered all 1,025 candidate-clade proteins across 417 taxa and 13 families
(482,007 residues). Each exact string matched both the checksum-bound staged
reconciliation FASTA and its independently identified QC proteome record. A
separate export readback checked all expected family/gene identities and every
sequence hash and length. Duplicate sequences and all copies are retained.

| Family | Proteins | Minimum length | Maximum length | Unique sequences | Noncanonical proteins |
| --- | ---: | ---: | ---: | ---: | ---: |
| OG0000054 | 11 | 218 | 248 | 11 | 0 |
| OG0000095 | 83 | 213 | 749 | 82 | 0 |
| OG0000107 | 9 | 942 | 975 | 9 | 0 |
| OG0000152 | 7 | 159 | 496 | 7 | 0 |
| OG0000230 | 24 | 220 | 762 | 19 | 0 |
| OG0000294 | 9 | 460 | 1091 | 9 | 0 |
| OG0000336 | 10 | 245 | 433 | 10 | 0 |
| OG0000972 | 622 | 50 | 1009 | 616 | 1 |
| OG0001082 | 111 | 205 | 921 | 111 | 0 |
| OG0001200 | 29 | 148 | 595 | 28 | 0 |
| OG0001203 | 21 | 110 | 308 | 21 | 0 |
| OG0002650 | 77 | 55 | 483 | 77 | 0 |
| OG0002812 | 12 | 770 | 817 | 12 | 0 |

Noncanonical record: `F157183_CAD6936802.1` in OG0000972, symbols/counts `{"X": 3}`. The original sequence is preserved; its alignment and ancestral-state treatment must be explicit.

Inputs are in `results/ancestral/case-sequences-20260927-v1/`. Reproduce with
`python scripts/recover_ancestral_case_sequences.py`; it refuses to overwrite
existing outputs. Closure: `metadata/case_ancestral_sequences_completed_20260927.json`.
This recovery used one CPU, streamed existing local FASTAs, and required no GPU
or network access. Alignment, outside-clade context and ancestral ensembles
remain pending.

## Two-method alignment run started

All 13 candidate sets are queued through MAFFT 7.525 L-INS-i (`--amino
--localpair --maxiterate 1000 --thread 4 --threadit 0`) and FAMSA 2.2.3
(`-t 4 -keep-duplicates`). Runs are sequential with four threads each, an 8 GiB
service memory limit and no swap. The 26-alignment planning envelope is
0.25–24 hours and 1 GiB output; no GPU or paid resources are used.

Every output must contain exactly the input identifiers and reproduce every
original sequence after gap removal, including the three X residues. No copies
are removed, including identical sequences. Column occupancies at 50%, 70% and
90% are descriptive exports; no trimming is applied and occupancy is not a
measure of alignment homology correctness. The first family passed preservation
checks under both methods at launch verification; remaining runs are pending.

`scripts/run_ancestral_case_alignments.py` uses the frozen
`metadata/ancestral_case_alignment_plan_20260927.json`; exact process identity
is in `metadata/ancestral_case_alignment_launch_20260927.json`. Completed runs
can be reused only with matching settings and checksums. An unfinished run
directory is preserved and requires explicit recovery rather than silent
overwriting. Full independent output audit, residue-correspondence comparison,
domain/fragment assessment and reconstruction sampling decisions remain pending.

## Full alignment readback and correspondence comparison queued

The downstream checker waits for successful terminal completion of all 26
alignments, validating the producer identity while it runs. It will independently
reconstruct every original sequence and residue position, every coverage column
and all threshold counts. All IDs, source receipts and artifact checksums must
match the frozen input and execution records.

For each family, exact full-column correspondence between MAFFT and FAMSA uses
the complete vector of residue positions (or gaps) across all proteins. This
criterion is deliberately stringent: one changed protein can break a full-column
match. All columns, occupancies and matched alternative column numbers are
retained. Focal duplicate residue-pair correspondences are compared separately
so that full-clade and focal-pair sensitivity can be distinguished. Neither
agreement fraction is a probability that a column is homologous. No trimming
is applied based on these diagnostics.

The checker uses one CPU and 4 GiB RAM, no swap, with a 0.5 GiB output allowance
and a planned 0.1–2 hours after alignment completion. Script:
`scripts/audit_ancestral_case_alignments.py`. Plan and exact live identity:
`metadata/ancestral_case_alignment_audit_plan_20260927.json` and
`metadata/ancestral_case_alignment_audit_launch_20260927.json`.
Full validation and alignment-sensitivity results remain pending.

## Domain context reviewed for all candidate proteins

All 1,025 proteins were joined by taxon/protein identifier and exact sequence
hash to the fully audited Pfam candidate architecture database. All four overlap
policies were retained (4,100 protein-policy rows and 52 family-policy rows).
The policies agree for every selected protein in this subset.

| Family | Proteins | With focal Pfam hit | Without focal hit | Partial focal HMM hit | Ordered annotation patterns |
| --- | ---: | ---: | ---: | ---: | ---: |
| OG0000054 | 11 | 11 | 0 | 0 | 1 |
| OG0000095 | 83 | 83 | 0 | 34 | 3 |
| OG0000107 | 9 | 9 | 0 | 0 | 1 |
| OG0000152 | 7 | 6 | 1 | 0 | 3 |
| OG0000230 | 24 | 24 | 0 | 0 | 3 |
| OG0000294 | 9 | 9 | 0 | 0 | 2 |
| OG0000336 | 10 | 10 | 0 | 0 | 1 |
| OG0000972 | 622 | 563 | 59 | 4 | 8 |
| OG0001082 | 111 | 85 | 26 | 5 | 14 |
| OG0001200 | 29 | 13 | 16 | 1 | 4 |
| OG0001203 | 21 | 14 | 7 | 0 | 2 |
| OG0002650 | 77 | 19 | 58 | 2 | 2 |
| OG0002812 | 12 | 12 | 0 | 0 | 1 |

Across these sets, 167 proteins lack a retained focal Pfam annotation, 85 have
no retained annotation, and 46 have a focal hit covering less than 70% of the
Pfam HMM. Partial model coverage does not prove a truncated protein. Likewise,
a missing retained hit is not proof of domain absence or evolutionary loss.

The focal pairs in OG0000230 (Furculomyces) and OG0001082 (Leucosporidium) have
different ordered candidate annotations. These cases need explicit separation
of domain-content differences from changes inside homologous domains. No
proteins were removed, and no ancestor or final reconstruction eligibility was
assigned from these summaries. Annotation patterns retain Pfam types and repeats
and do not establish functional equivalence.

Full annotations and policy-specific coordinates are retained under
`results/ancestral/case-domain-context-20260927-v1/`. Reproduce with
`python scripts/annotate_ancestral_case_domains.py`. Exact source database hash,
prior full readback, protein-sequence hashes and output hashes are recorded.
Closure: `metadata/ancestral_case_domain_context_completed_20260927.json`.
These results motivate domain-aware alignment assessment before ancestral
reconstruction, alongside the ongoing two-method full-sequence comparison.

## Domain-specific sequence inputs completed

Prepared 26 FASTAs: all 13 cases under both Pfam alignment-span and envelope
boundaries. The 858 proteins with exactly one retained focal hit contribute
1,716 domain records and 318,670 residues across the two definitions. Both focal
duplicates are present in every set. All 46 proteins with partial focal HMM
coverage remain included and flagged, and identical sequences are not collapsed.

Every exported residue maps to the corresponding original protein position.
Serialized coordinate readback reconstructed every domain string, and a separate
check verified all FASTA identities, lengths and sequence hashes. The complete
2,050-row protein/boundary disposition grid includes the 167 proteins without a
retained focal hit under each boundary. No multiple-hit ambiguity occurred in
this input set; the script would retain that disposition without silently
choosing a domain copy.

These are domain-sensitivity inputs. They are not a claim that missing annotations
represent true losses, nor a replacement for whole-protein context or ancestral
indel uncertainty. Domain alignment, gene-tree context after restricting tips,
and ancestral reconstruction eligibility remain pending.

Reproduce: `python scripts/prepare_ancestral_domain_sequences.py`. Artifacts:
`results/ancestral/case-domain-sequences-20260927-v1/`, including full source
coordinates and checksums. Closure:
`metadata/ancestral_domain_sequences_completed_20260927.json`.

## Domain alignments launched under both boundaries and methods

All 52 combinations (13 cases × two Pfam boundaries × MAFFT/FAMSA) are queued
with the same alignment settings as the whole-protein sensitivity run. The
largest input contains 563 proteins and the longest extracted domain is 331
residues. Each output must preserve every original domain residue and protein
copy, including partial HMM hits, and records all column occupancies. No
sequence trimming, deduplication or residue replacement is applied.

The separate service uses four CPUs and at most 8 GiB RAM, no swap and no GPU.
The planning envelope is 0.25–24 hours and 1 GiB output. Exact script/settings
and source hashes are frozen in `metadata/ancestral_domain_alignment_plan_20260927.json`;
process identity is in `metadata/ancestral_domain_alignment_launch_20260927.json`.
Run script: `scripts/run_ancestral_domain_alignments.py`.

At launch verification, 28 alignments had passed initial preservation checks.
Production is not complete; independent readback, comparison of residue
correspondence across both boundaries/methods, and comparison with whole-protein
alignments remain pending. No ancestral sequence or structure has been inferred.

## Full domain alignment verification and protein-coordinate maps queued

A separate checker waits for successful terminal completion of all 52 domain
alignments, checking the producer PID, creation time and command while it runs.
It then checks all input identifiers, original sequence strings, residue
positions, column occupancies, threshold counts and output hashes. Both methods
must preserve every domain copy and input residue.

The checker reconstructs cumulative domain positions and independently enumerates
nongap sequence positions. Every domain position is then translated through its
recorded boundary offset to an original protein position, with the residue
checked against the full-protein sequence. All protein-position maps are
exported, providing common coordinates for downstream comparisons.

Within each of the 26 family/boundary sets, full-column position vectors and
focal duplicate residue pairs are compared between MAFFT and FAMSA. The current
comparison exports identify sets as `family-boundary`; residue-pair coordinates
are absolute protein coordinates. A whole-column match requires agreement for
every protein and is not a probability of alignment correctness. Comparisons
across boundary definitions or against whole-protein alignments remain separate
downstream work.

Checker: `scripts/audit_ancestral_domain_alignments.py`. Frozen plan and exact
launch identity: `metadata/ancestral_domain_alignment_audit_plan_20260927.json`
and `metadata/ancestral_domain_alignment_audit_launch_20260927.json`. Resources:
one CPU, 4 GiB RAM, no swap, 0.5 GiB output allowance, planned 0.1–2 hours after
production finishes. No completed readback or ancestral inference is claimed yet.

## All 52 domain alignments and full readback completed

Both the domain producer and checker terminated inactive/success/exit-zero.
All 52 alignments passed full sequence/copy preservation, column occupancy and
protein-coordinate checks: 637,340 mapped residue records, 9,349 columns and
3,828 focal-pair union records. The latter counts include both boundaries and
methods and must not be interpreted as independent observations.

MAFFT and FAMSA give identical focal duplicate pair maps in 21 of 26
family/boundary sets. They differ for Rhodonia alignment-span, Cryoendolithus
alignment-span, Jaapia envelope and both Piloderma boundaries. Even when focal
pair maps agree, other sequences can yield substantial full-clade correspondence
differences. For Neolecta, all 213 alignment-span and 225 envelope focal pairs
agree, while only 45.8–46.4% and 39.7–40.6% of full-clade columns respectively
match exactly between methods. The stringent full-column measure depends on
clade size and is not directly comparable as an error rate across families.

| Family/boundary | MAFFT columns | FAMSA columns | Shared full columns | Shared focal pairs / union |
| --- | ---: | ---: | ---: | ---: |
| OG0000054-alignment | 90 | 90 | 88 | 74 / 74 |
| OG0000054-envelope | 95 | 95 | 93 | 94 / 94 |
| OG0000095-alignment | 147 | 160 | 90 | 108 / 122 |
| OG0000095-envelope | 151 | 161 | 85 | 125 / 125 |
| OG0000107-alignment | 112 | 112 | 106 | 91 / 99 |
| OG0000107-envelope | 112 | 112 | 108 | 100 / 100 |
| OG0000152-alignment | 130 | 130 | 130 | 114 / 114 |
| OG0000152-envelope | 133 | 133 | 133 | 119 / 119 |
| OG0000230-alignment | 125 | 124 | 121 | 108 / 108 |
| OG0000230-envelope | 125 | 124 | 121 | 116 / 116 |
| OG0000294-alignment | 231 | 231 | 231 | 224 / 224 |
| OG0000294-envelope | 232 | 232 | 226 | 226 / 236 |
| OG0000336-alignment | 147 | 147 | 141 | 136 / 136 |
| OG0000336-envelope | 148 | 148 | 142 | 137 / 137 |
| OG0000972-alignment | 397 | 392 | 182 | 213 / 213 |
| OG0000972-envelope | 453 | 443 | 180 | 225 / 225 |
| OG0001082-alignment | 105 | 106 | 88 | 72 / 72 |
| OG0001082-envelope | 114 | 117 | 80 | 73 / 73 |
| OG0001200-alignment | 168 | 167 | 114 | 123 / 123 |
| OG0001200-envelope | 200 | 204 | 107 | 164 / 164 |
| OG0001203-alignment | 89 | 90 | 81 | 72 / 72 |
| OG0001203-envelope | 94 | 95 | 86 | 74 / 74 |
| OG0002650-alignment | 194 | 192 | 178 | 174 / 180 |
| OG0002650-envelope | 218 | 216 | 202 | 186 / 192 |
| OG0002812-alignment | 327 | 328 | 317 | 316 / 316 |
| OG0002812-envelope | 331 | 332 | 318 | 320 / 320 |

Every displayed shared-column count/fraction and focal shared/union count was
reaggregated from the serialized detailed comparison tables. These are
correspondence-sensitivity diagnostics, not evidence that one alignment is
correct or that an ancestral state is certain. Both methods remain in scope.

Results: `results/ancestral/case-domain-alignment-readback-20260927-v1/`.
Closure: `metadata/ancestral_domain_alignments_completed_20260927.json`.
Whole-protein alignments, cross-boundary comparison, topology/node eligibility
and ancestral sequence ensembles remain incomplete.

## Local full-protein and domain tree inputs completed

Prepared 52 local trees: 13 cases × two reconciliation guides × full/domain
tip sets. Every tree retains its original level-three local root and the exact
focal duplicate MRCA. All 208 candidate ancestral mappings are unique, and no
two of the four examined ancestral levels collapse to the same output node.
All 4,100 original tip dispositions are explicit.

Every retained rooted descendant set and edge length was checked against the
source tree after serialization. Pruned paths sum their constituent source
lengths; all 7,420 exported edge records match, with maximum difference
5.55e-17. The local root stem is zero because it lies
outside all within-clade paths.

Across guides, 26 of 26 corresponding full/domain tree pairs have identical
rooted edge sets. Maximum shared-edge length difference is 0. These
are sensitivity comparisons of reconciliations sharing original gene trees,
not independent estimates or replication. Full pairwise results are retained
in the completion record.

The trees provide candidate topologies and node identities, not finished ASR
models. Inherited source branch lengths must be re-estimated against the
corresponding full-protein or domain alignment. Model choice, alignment and
topology uncertainty, informative-site coverage, indel treatment, and additional
sampling/root sensitivity remain required before interpreting ancestors.

Reproduce: `python scripts/prepare_ancestral_case_trees.py`. Local artifacts:
`results/ancestral/case-local-trees-20260927-v1/`. Completion record:
`metadata/ancestral_case_local_trees_completed_20260927.json`.

## Domain substitution-model and branch-length fitting launched

All 52 domain alignments are being fitted under LG+F+G4, WAG+F+G4 and JTT+F+G4
with IQ-TREE 3.0.1, producing 156 model/alignment combinations. The empirical
exchangeability models share observed amino-acid frequencies and four discrete
gamma rate categories. All models remain retained; no best model, adequacy or
ancestral-state conclusion is assigned at this stage.

Each fit fixes the local unrooted topology (`-t ... --tree-fix`) while estimating
model parameters and branch lengths from that alignment. `-keep-ident` retains
identical sequence copies. Output tips, nonnegative finite branch lengths and
unrooted splits are checked against the input; the reported likelihood must
be finite and the model label must match. The full likelihood/parameter audit
and convergence assessment remain downstream.

The profile-domain tree is used once because every corresponding MAFFT-guide
local tree has identical rooted edge sets and source lengths. This saves
duplicate computation while retaining both guides' provenance; it does not
create independent replicate estimates or resolve topology/root uncertainty.
Fits under different boundary definitions or alignment methods have different
data and their raw likelihoods must not be compared as a model-selection test.

Two fits run concurrently, four threads and a 4 GiB IQ-TREE memory allowance
each. The service has an eight-CPU quota, 12 GiB memory limit and no swap.
Planning envelope: 1–48 hours, 2 GiB output, existing local CPU only. First-family
production checks passed; the full run is still active at launch verification.

Script: `scripts/run_ancestral_domain_model_fits.py`. Frozen plan and exact
launch identity: `metadata/ancestral_domain_model_fit_plan_20260927.json` and
`metadata/ancestral_domain_model_fit_launch_20260927.json`. No ancestral
sequences, ancestral structures or indel states have been inferred.

## All 156 model fits and report/checkpoint readback completed

The producer terminated successfully. The full checker verified all 156 model
labels, input/output tips and unrooted edge sets; checkpoint/report likelihood
agreement; checkpoint/tree branch lengths; empirical frequencies; gamma category
means; and AIC/BIC arithmetic. Model parameters and all flagged diagnostic lines
are retained under `results/ancestral/case-domain-model-readback-20260927-v1/`.

The frequency check initially assumed normalized nongap counts. Inspection of
the [tagged IQ-TREE 3.0.1 source](https://github.com/iqtree/iqtree3/blob/v3.0.1/alignment/alignment.cpp)
shows eight ambiguity-allocation iterations from uniform frequencies. Recreating
that calculation, including gaps, matches all reported frequencies to rounding
precision. The default `keep_zero_freq=true` is documented in the tagged
[parameter initialization](https://github.com/iqtree/iqtree3/blob/v3.0.1/utils/tools.cpp).
Fetched sources and hashes are retained locally; this is version-specific
numerical agreement, not proof of the executable's full build provenance.
Raw nongap count frequencies and reconstructed IQ-TREE frequencies are exported
separately. The original fitted outputs were not changed.

Twelve fits have unresolved input topologies, so parameter counts were checked
using the actual retained branch count plus 19 frequencies and one gamma shape
parameter, rather than assuming a fully bifurcating tree. All parameter-count
and information-criterion checks passed.

Important remaining diagnostics:

- 135 fits report near-zero internal branches; 140 contain an edge below 1e-5.
- 45 fits report sequences with more than 50% gaps/ambiguity.
- 12 report rare amino-acid states.
- Gamma shape estimates range from 0.0200255 to 1.89305, including estimates near
  the lower search boundary that need further assessment.

The diagnostic export also retains composition-test lines containing “failed,”
including summaries reporting zero failures. Its 624 lines and 156 fits with
matching diagnostic text must not be read as 156 failed fits. Composition-test
failures and warning types require their own interpretation.

This closes report/serialization checks only. It does not independently recompute
the phylogenetic likelihood, establish convergence or model adequacy, or provide
ancestral sequences. Likelihood replay, boundary/optimization sensitivity and
node-specific information assessment remain before interpreting ancestors.
Reproduce: `python scripts/audit_ancestral_domain_model_fits.py`. Closure:
`metadata/ancestral_domain_model_reports_completed_20260927.json`.

## Independent likelihood calculation passed all 156 domain fits

A separately implemented scaled-pruning calculation reproduced every saved
IQ-TREE domain likelihood. Maximum absolute difference was 1.25591e-5 log
likelihood units, below the declared 0.001 tolerance for finite checkpoint
parameter/branch serialization. All 156 signed differences remain retained.
The service terminated inactive/success/exit-zero.

The implementation constructs reversible rate matrices from the versioned
[LG, WAG and JTT exchangeabilities](https://github.com/iqtree/iqtree3/blob/v3.0.1/model/modelprotein.cpp),
recreates version-specific empirical frequencies and gamma category means,
and computes transition matrices using SciPy. Partial likelihoods are scaled
through each child contribution to avoid underflow; rate categories are mixed
with log-sum-exp. Very short branch/rate products use direct matrix exponentials
to avoid eigendecomposition cancellation. An analytic two-tip mixture identity
with an unknown/gap state passed before the full calculation.

This implementation shares input alignments, empirical rate tables, serialized
trees and parameter definitions with IQ-TREE, but not its likelihood engine.
Gaps are marginalized as unknown states; no ancestral deletion state is inferred.
Numerical agreement validates evaluation at the saved parameters, not optimum
convergence, model adequacy, topology certainty or biological ancestral states.

Reproduce with `python scripts/replay_ancestral_domain_likelihoods.py` after
materializing the sources in `metadata/ancestral_domain_likelihood_replay_plan_20260927.json`.
The run used one CPU and at most 4 GiB RAM, no swap or GPU, against a 0.1–4 hour
planning envelope. Output:
`results/ancestral/case-domain-likelihood-replay-20260927-v1/`. Closure:
`metadata/ancestral_domain_likelihood_replay_completed_20260927.json`.

## Candidate nodes mapped into fitted unrooted trees

All candidate node identities were checked using the complete partition of tip
identities around each incident branch, independently of node labels. Across
156 model/alignment fits, the focal duplicate MRCA and its next two ancestral
levels map uniquely: 468 distinct fit/node combinations. Both reconciliation
guides give the same partitions and mappings; the 936 matched guide/node rows
are aliases, not independent reconstructions.

The source local root has two incident branches. Its position is suppressed in
the fitted unrooted tree and is not identifiable under these reversible models.
All 156 such fit/root candidates (312 guide rows) are explicitly marked
`degree_two_root_position_not_identified`. No arbitrary branch midpoint or
nearby vertex substitutes for that root. This limits which nodes can receive
conditional ancestral-state calculations with the current inputs.

The complete 1,248-row mapping retains branch-partition signatures, source/fitted
labels, original levels and incident branch lengths. Short adjacent branches
remain flagged by their numeric lengths, without collapsing nodes or treating
small distances as support for a historical ancestor.

Reproduce: `python scripts/map_ancestral_fitted_nodes.py`. Local output:
`results/ancestral/case-fitted-node-mapping-20260927-v1/`. Closure:
`metadata/ancestral_fitted_node_mapping_completed_20260927.json`.
This stage establishes correspondences only. Conditional ancestral states,
optimization and model sensitivity, indel uncertainty and predicted ancestral
structures remain incomplete.

## Conditional ancestral amino-acid marginals

The complete 156-fit grid now has a CPU calculation for the three mapped
internal vertices per fit (468 fit/node combinations). The calculation retains
all 20 amino-acid probabilities at every alignment column, including columns
with incomplete descendant coverage. Gaps supply unknown-state evidence, not
evidence of an ancestral deletion. Descendant and outside-clade observed
residue counts accompany every site. No final ancestral sequence or structure
is selected at this stage.

The implementation reroots the likelihood traversal at each requested vertex,
uses scaled pruning, and integrates four gamma rate categories jointly with
amino-acid states. Analytic three-tip and all-unknown checks precede production.
Every node calculation must reproduce the saved fit likelihood within 0.001
log units, and site likelihoods must agree across the three traversals within
1e-8. These checks establish internal consistency, not independent verification
of every posterior or adequacy/convergence of the fitted evolutionary models.

Run `python scripts/infer_conditional_domain_ancestors.py` using the pinned
`metadata/conditional_domain_ancestor_plan_20260927.json`. Outputs are in
`results/ancestral/conditional-domain-ancestors-20260927-v1/`: one compressed
probability array per fit, a complete site summary, and a node summary. The
resource plan is one CPU, 4 GiB RAM, no swap/GPU, up to 2 GiB output and
0.1–4 hours. Launch identity is recorded in the corresponding launch metadata.
Independent posterior verification, optimization sensitivity, cross-alignment
site comparisons, indel uncertainty and ancestral structure prediction remain
required downstream work.

The producer terminated successfully: all 468 nodes and 84,141 node/site
distributions completed. Maximum saved-likelihood difference was 1.25591e-5
log units. Complete artifact hashes, probability normalization/bounds, and every
serialized site/MAP correspondence passed readback. Completion record:
`metadata/conditional_domain_ancestors_completed_20260927.json`. This does not
close the independent posterior-validation or model-sensitivity requirements.

## Independent posterior calculation completed

All 1,682,820 amino-acid probabilities across 84,141 node/sites and 468 mapped
vertices passed a separate fixed-root inside/outside calculation. It evaluates
full columns without pattern compression and uses direct matrix exponentials
rather than the producer's eigendecomposition and traversal rerooting. The
largest probability difference was 4.39669e-12; the largest per-site
log-likelihood difference was 2.44825e-10, both below the prespecified 1e-8
tolerances. A two-internal-node analytic enumeration with missing observations
also passed. All 156 fits were checked without exclusions.

Both implementations share input alignments, trees, saved parameters and
empirical exchangeabilities. The audit therefore checks numerical calculation,
not model convergence, adequacy, ancestral polarity or indel history.
Reproduce with `python scripts/audit_conditional_domain_posteriors.py`; pinned
plan, launch and completion metadata use the prefix
`conditional_domain_posterior_audit_` and date `20260927`. The terminal run used
one CPU, 4 GiB memory limit and no swap/GPU.

## Whole-protein alignment comparison completed

All 26 whole-protein alignments and their independent readback terminated
successfully. The audit checked every sequence/copy and residue coordinate,
24,191 alignment columns and 6,287 focal-pair union records. Only four families
have identical focal-pair mappings between MAFFT and FAMSA: OG0000054,
OG0000230, OG0001200 and OG0002812. All nine remaining families retain
alignment-dependent mappings. Both alignments remain available for downstream
sensitivity analysis; matching focal pairs does not establish full-clade
column equivalence.

The largest family, OG0000972, has 2,818 MAFFT versus 2,212 FAMSA columns, with
474 exact shared full-clade columns. These stringent identity counts are not
probabilities of homology. Whole-protein ancestral fitting and comparisons
between whole-protein and extracted-domain correspondences remain pending.
Completion: `metadata/ancestral_case_alignments_completed_20260927.json`.

## Sensitivity of conditional ancestral probabilities

Compared all 66 pairs of the 12 model/boundary/aligner settings per family:
858 context pairs, 2,574 node comparisons and 551,997 site-union records.
Matching requires the exact same original protein coordinate (or gap) for
every member of the common clade. The complete signature map is retained.
There are 373,554 matched site comparisons; unmatched records are explicit,
not counted as agreement and not assigned a posterior distance.

For comparisons changing one factor only:

| Changed factor | Matched node/site comparisons | Different MAP residues | Unmatched column records |
| --- | ---: | ---: | ---: |
| Substitution model | 84,141 | 1,578 | 0 |
| Aligner | 33,732 | 146 | 16,677 |
| Domain boundary | 32,724 | 278 | 18,693 |

No comparison in the complete grid assigns at least 0.90 probability to both
of two different MAP residues. This threshold is descriptive, not a validation
criterion. Posterior differences can nevertheless be substantial: maximum
total-variation distance reaches 0.6493 for model-only comparisons. Full
sitewise values and summaries are retained, including all combinations of
changed factors. Comparisons reuse nodes/sites and are not independent
replicates or counts of unique evolutionary changes. MAP switches near ties
also need not imply a meaningful biological difference.

Full saved-table readback reproduces all summary counts and distances.
Reproduce: `python scripts/compare_ancestral_domain_probabilities.py`; plan,
launch and completion metadata use `ancestral_domain_probability_comparison_`
and date `20260927`. Local output:
`results/ancestral/ancestral-domain-probability-comparison-20260927-v1/`.
One CPU, 4 GiB memory limit, no swap/GPU. These comparisons remain conditional
on the current fixed trees and optimized parameters. Optimization sensitivity,
whole-protein context, indel uncertainty and ancestral structural ensembles
remain unfinished.

## Finite likelihood-neighborhood diagnostic running

The full 156-fit set now has a CPU diagnostic evaluating 15 parameter
combinations each (2,340 points): fitted gamma shape multiplied by
0.5, 0.9, 1, 1.1 or 2 and all branch lengths multiplied by 0.9, 1 or 1.1.
Frequencies and topology remain fixed. The original point must reproduce its
checkpoint likelihood within 0.001 log units before any perturbations are
accepted. The same threshold flags improvements. All signed changes are
retained, not just improving points.

Some gamma perturbations extend below IQ-TREE's default lower bound of 0.02,
confirmed in the cached version-specific source (`tools.cpp`, line 6399;
checksum recorded with the launch). Those points are labeled separately:
improvement outside the original feasible domain does not itself demonstrate
a failure of the original constrained optimization. Improvement inside it
would identify fits needing further optimization. Failure to find improvement
on this finite grid cannot establish convergence or global optimality.

Run `python scripts/diagnose_ancestral_domain_fit_neighborhoods.py`; plan and
launch use prefix `ancestral_domain_fit_neighborhood_`, date `20260927`.
Output: `results/ancestral/domain-fit-neighborhoods-20260927-v1/`.
One CPU, 4 GiB memory limit, no swap/GPU; estimate 0.1–6 hours and at most
0.2 GiB output. Verified live at launch; full results remain pending.
Independent branch optimization, alternate starts and posterior propagation
after refitting remain required.

## Alternate-start domain model refits running

All 156 baseline fits are being refitted from three gamma/branch starts:
(alpha 0.05, branch scale 0.5), (alpha 0.5, scale 1), and (alpha 2, scale 2).
Each runs with gamma lower bound 0.02 and 0.005, totaling 936 fits. Branch
scaling is applied to each baseline fitted tree; unrooted topology and all
tip identities are retained. The nonexistent root stem is explicitly zero.
All exported tree labels and lengths are checked after serialization.

IQ-TREE 3.0.1 receives `-a <start> -optfromgiven --alpha-min <bound>` and
`--epsilon 0.000001`, with branch lengths and gamma shape free to optimize.
The cached version-specific `model/rategamma.cpp` confirms that positive
starting alpha is optimized when `optimize_from_given_params` is enabled;
its hash and source URL are pinned. The first completed fit moves alpha from
0.05 to about 1.2053 and retains 39 free parameters, confirming runtime
behavior for that job. Full output qualification remains pending.

Preparation: `python scripts/prepare_ancestral_domain_multistarts.py`.
Execution: `python scripts/run_ancestral_domain_multistarts.py --plan
metadata/ancestral_domain_multistart_plan_20260927.json`. Inputs are under
`results/ancestral/domain-multistart-inputs-20260927-v2/`; outputs under
`results/ancestral/domain-multistarts-20260927-v1/`. Two concurrent four-thread
fits, aggregate eight CPUs and 12 GiB memory, no swap/GPU. Planning envelope:
1–48 hours and up to 8 GiB output. Launch identity and an early failed
preparation/launch attempt are preserved in launch metadata. Original fits
and posterior calculations are unchanged.

Completion will require every fit's input/tip/topology/parameter audit,
likelihood replay and comparison of solutions across starts and bounds.
Posterior sensitivity must then be propagated from qualified solutions.
Alternate-start agreement alone cannot establish model adequacy or global
optimality.

## Complete alternate-start audit queued

A separately launched one-CPU checker waits on the exact refit producer
PID, creation time and command line, then requires terminal success. It
checks all 936 input/output hashes, tip and edge identities, gamma parameters,
eight-iteration empirical frequencies, gamma-category means, report/checkpoint
likelihoods, free-parameter counts and information-criterion arithmetic.
Every saved likelihood is independently recomputed with the previously
validated scaled-pruning engine (absolute tolerance 0.001 log units).

All three starts are compared within each of the 312 baseline/bound groups,
retaining likelihood ranges, best-fit identifiers and changes from baseline.
A second table compares bounds for every baseline fit. The checker does not
require agreement between starts or an improved fit: discrepancies remain
results for further optimization and posterior sensitivity, not omitted cases.

Reproduce with `python scripts/audit_ancestral_domain_multistarts.py`; plan
and launch metadata use prefix `ancestral_domain_multistart_audit_` and date
`20260927`. Output: `results/ancestral/domain-multistart-readback-20260927-v1/`.
One CPU, 4 GiB memory limit, no swap/GPU; estimated 0.1–6 hours after refits,
up to 0.2 GiB output. Producer and checker verified live; audit results pending.

The finite-neighborhood diagnostic subsequently completed and passed full
2,340-row readback. None of the 156 fits improved by more than 0.001 log units
at the specified grid points (maximum increase 2.14938e-6). This result only
concerns the tested gamma/uniform-branch directions; the full alternate-start
branch refits remain active and may find improvements outside those directions.
Closure: `metadata/ancestral_domain_fit_neighborhood_completed_20260927.json`.

## Whole-protein model fitting launched

All 26 independently checked whole-protein alignments are now fitted under
LG+F+G4, WAG+F+G4 and JTT+F+G4 (78 fits), with branch lengths and gamma shape
optimized at epsilon 1e-6 on the complete local whole-protein topologies.
The two guide aliases have identical tree bytes for every family and are
fitted once. All 1,025 original proteins remain represented. The three X
residues in one OG0000972 protein are retained as unknown amino-acid evidence;
no sequence was discarded or silently changed. The downstream frequency and
likelihood audit must explicitly handle X as well as alignment gaps.

These whole-protein fits include 167 proteins without a retained focal domain
annotation, unlike the domain-only fits. Therefore whole/domain comparisons
must account for both sequence context and taxon/copy membership differences;
we cannot attribute every contrast to domain boundaries alone. Alignment
correspondence sensitivity and insertion/deletion uncertainty remain relevant.

Prepare with `python scripts/prepare_ancestral_whole_model_fits.py`; execute
`python scripts/run_ancestral_whole_model_fits.py --plan
metadata/ancestral_whole_model_fit_plan_20260927.json`. Output:
`results/ancestral/whole-protein-model-fits-20260927-v1/`. Plan and launch
metadata freeze all inputs and record the live process. Two concurrent
four-thread fits, aggregate eight CPUs and 12 GiB memory, no swap/GPU;
planning estimate 1–48 hours and up to 4 GiB output. All final model reports,
likelihoods, ancestor mappings and probabilities remain to be audited or
inferred before any ancestral sequence is selected.

## Whole-protein report and likelihood audit queued

The all78-fit checker is now waiting on the exact live whole-protein producer,
requiring successful terminal exit before reading its complete receipt. Every
model label, tip identity, unrooted edge, checkpoint/report likelihood, branch
length, gamma category, empirical frequency and AIC/BIC calculation is checked.
It then independently recomputes every likelihood with the scaled-pruning
engine, requiring absolute agreement within 0.001 log units. Warnings and
short branches are retained.

Protein X and alignment gaps are both unknown evidence in IQ-TREE3.0.1
(`alignment.cpp`, protein-state conversion at lines1712–1728 and gap handling
at1656–1657). Both enter the eight-iteration frequency reconstruction as
unknown cells. For likelihood replay only, X is mapped to the helper's unknown
marker; source FASTAs and residue coordinates are unchanged. Every fit's X
count must match its frozen input job. Other ambiguity symbols are rejected
because they are absent from this input set and need different evidence masks.

Run `python scripts/audit_ancestral_whole_model_fits.py`; plan and launch
metadata use prefix `ancestral_whole_model_audit_` and date `20260927`. Output:
`results/ancestral/whole-protein-model-readback-20260927-v1/`. One CPU,4GiB RAM,
no swap/GPU; planning0.1–6 hours after production and0.2GiB output. Full results
remain pending. Passing numerical checks will not establish optimization
convergence, model adequacy or ancestral states.

## Whole-protein/domain alignment correspondence completed

Projected each whole-protein alignment onto the exact annotated domain
intervals and protein membership of each domain input. All original domain
residues were recovered with exact identity and coordinate checks. Every
retained column is a full-clade vector of original residue positions or gaps;
columns without any retained domain residue are counted separately. Compared
all13 families ×2 boundaries ×2 whole aligners ×2 domain aligners =104
combinations, retaining23,037 matched/unmatched column-union records.

Only26 of104 comparisons have identical projected column sets. OG0000152
is identical in all eight settings. For OG0000972, the fraction of domain
columns with an exact projected whole-protein match ranges from0.339 to0.615;
for OG0000095 it ranges from0.497 to0.687. These are coordinate-correspondence
fractions, not homology probabilities, and comparisons reuse data. All
settings remain retained, including low-correspondence cases.

Projection controls the members and residue intervals being compared; it
does not remove the influence of extra proteins or flanking sequences on the
original whole-protein alignment. Consequently it does not isolate sequence
context effects from membership effects. Matched ancestral-column comparisons
must retain that limitation, and unmatched positions require separate
alignment/indel uncertainty handling.

Reproduce: `python scripts/compare_whole_domain_ancestral_alignments.py`.
Plan and completion metadata use `whole_domain_ancestral_alignment_comparison_`
and date `20260927`. Outputs:
`results/ancestral/whole-domain-alignment-correspondence-20260927-v1/`.
One CPU,4GiB memory limit,noGPU; terminal success and complete saved-record
identity/count readback confirmed.

## All-refit ancestral probability propagation queued

All936 domain alternate-start/bound fits have finished production; their
independent audit is actively checking reports and recomputing likelihoods.
A downstream CPU calculation now waits on the exact audit process and
requires terminal success, the expected936-fit audit receipt, its frozen
plan hash and all artifact hashes before using fitted parameters.

It will compute probabilities at the three mapped internal vertices for all
six refits per baseline model/alignment, totaling2,808 fit/node combinations.
Every setting remains retained; no only-best-fit selection hides differences
between starts or bounds. Target identities are recovered from complete
incident-tip partition signatures, independently of changed internal labels.
Every traversal must reproduce its fitted likelihood; independent posterior
validation and comparisons to baseline remain downstream requirements.

Run `python scripts/infer_refitted_domain_ancestors.py`. Frozen plan/launch
metadata use `refitted_domain_ancestor_`, date `20260927`; output:
`results/ancestral/refitted-domain-ancestors-20260927-v1/`. One CPU,4GiB memory
limit,no swap/GPU; estimate0.2–8 hours after audit and up to4GiB output.
These are conditional residue probabilities, not final ancestral sequences
or evidence for ancestral deletion states. Model adequacy, indel uncertainty
and structural prediction remain open.

## Independent verification of all-refit probabilities queued

Queued the fixed-root inside/outside checker for every refitted domain
probability array. The expected scope is936 fits,2,808 mapped nodes,504,846
node/sites and10,096,920 probabilities. It waits for verified terminal success
of the refitted-probability producer, then checks its receipt against the
fit audit and fitted parameters before evaluating any array. Every node is
matched by complete incident-tip partitions, independently of node labels.
The second evaluator uses direct matrix exponentials and uncompressed
alignment columns; its analytic two-internal-node enumeration precedes
production. Prespecified absolute tolerances are1e-8 for both probabilities
and site log likelihoods. Every fit remains in scope.

Run `python scripts/audit_refitted_domain_posteriors.py`; plan/launch metadata
use `refitted_domain_posterior_audit_`, date `20260927`. Output:
`results/ancestral/refitted-domain-posterior-audit-20260927-v1/`. One CPU,4GiB
RAM,no swap/GPU; estimate0.2–8 hours after probability production,0.2GiB output.
Results remain pending. Numerical agreement does not establish model adequacy,
convergence, indel history or biological support for ancestral sequences.

## Alternate-start audit completed: optimization differences remain

All936 report/checkpoint checks and independent likelihood replays passed,
maximum absolute error1.23417e-5 log units. Complete saved group extrema and
all artifact hashes were checked; producer and auditor terminated successfully.
At each gamma bound,48 of156 model/alignment groups have across-start
likelihood range above0.001; maximum range14.2882. Best refits improve49 of156
baseline fits by more than0.001, maximum gain0.383662, at either bound.
No lower-bound best fit improves over the original-bound best fit by0.001
(maximum difference1.98e-5). These results do not establish global convergence.

Thus the earlier finite gamma/uniform-branch grid missed optimization
differences involving individual branches or other optimization paths.
Poorer local solutions must not be interpreted as equally plausible biological
uncertainty. All-refit posterior calculations remain useful diagnostics, but
selection or weighting for ancestral sequence ensembles requires convergence
assessment and comparison to the best qualified solutions. Short edges remain
common (840 of936 fits); warning exports also include informational zero-failure
lines and are not counts of failed jobs. Completion:
`metadata/ancestral_domain_multistart_audit_completed_20260927.json`.

## Best-solution refinement launched

The alternate-start audit found remaining optimization differences. For each
of156 baseline models under each of two gamma lower bounds, selected the
highest checkpoint likelihood among the three audited same-bound refits and
the original audited baseline when feasible. This produces312 starts. Twenty
selected starts are original baselines; small likelihood differences are not
treated as evidence of biological preference. Selection is deterministic and
all source job IDs, input-tree hashes, starting alphas and likelihoods remain
recorded. Every prior solution is retained.

Each selected tree/gamma is now freely refitted with IQ-TREE epsilon1e-8,
compared with1e-6 for the alternate starts. Topology, frequencies, gamma lower
bound and alignment stay fixed. This is a stability test of the best available
solutions, not a proof of global optimization. Full report/likelihood audits
and signed changes from selected starts are still required before qualifying
refined estimates. Ancestor probabilities currently running from the earlier
936 fits remain diagnostics of optimization sensitivity.

Prepare: `python scripts/prepare_ancestral_domain_refinements.py`.
Run: `python scripts/run_ancestral_domain_refinements.py --plan
metadata/ancestral_domain_refinement_plan_20260927.json`. Output:
`results/ancestral/domain-refinements-20260927-v1/`. Two concurrent four-thread
fits,8CPUs aggregate,12GiB memory,no swap/GPU; planning0.5–24 hours,4GiB output.
Plan/launch metadata freeze all selected starts and verified live identity.

## Refinement audit queued

The312 best-start refinements have finished production. A separate checker
requires successful producer termination and checks every report/checkpoint,
input hash, topology, tip identity, frequency and gamma parameter, information
criterion and independently recomputed likelihood (tolerance0.001).
It verifies each selected starting likelihood directly from its original
checkpoint, then retains signed likelihood changes, gamma-shape changes and
maximum branch-length changes for every refinement. Improving and worsening
solutions both remain represented. Stability is not global optimality.

Run `python scripts/audit_ancestral_domain_refinements.py`; plan/launch
metadata use `ancestral_domain_refinement_audit_`, date `20260927`. Output:
`results/ancestral/domain-refinement-readback-20260927-v1/`. One CPU,4GiB
memory,no swap/GPU; estimate0.1–6 hours after production,0.2GiB output.
Full audit results and refined ancestral probabilities remain pending.

## Refinement audit passed; refined ancestral probabilities launched

All312 refined fits passed full report/checkpoint checks and independent
likelihood replay, maximum error1.21704e-5 log units. None changed its selected
starting likelihood by more than0.001 in either direction: signed range
-3.00006e-8 to+0.0002558. Complete serialized differences and artifact hashes
were checked. This supports stability under the specified best-start restart
and tighter epsilon; it does not establish global optimality or model adequacy.
Closure: `metadata/ancestral_domain_refinement_audit_completed_20260927.json`.

Launched conditional amino-acid probabilities for all312 refined fits at three
signature-matched internal vertices each (936 fit/node combinations). Both
gamma bounds and every alignment/model combination are retained. Prior
936 alternate-start posteriors remain separate diagnostics, including poorer
solutions, rather than equally weighted biological alternatives.

Run `python scripts/infer_refined_domain_ancestors.py`; plan/launch metadata
use `refined_domain_ancestor_`, date `20260927`; output:
`results/ancestral/refined-domain-ancestors-20260927-v1/`. One CPU,4GiB RAM,
no swap/GPU; estimate0.1–4 hours,2GiB output. The complete audit receipt and
parameter file are frozen before launch. Full probability verification,
model/aligner/bound sensitivity, indel handling and sequence ensembles remain
incomplete.

## Alternate-start probabilities completed; refined verification queued

All936 alternate-start/bound probability calculations finished:2,808 mapped
fit/nodes and504,846 node/site distributions, retaining10,096,920 amino-acid
probabilities. Analytic and within-fit likelihood checks passed (maximum
likelihood difference1.23417e-5); full artifact hashes and terminal success
confirmed. Independent fixed-root posterior verification is running. Closure
for production only: `metadata/refitted_domain_ancestors_completed_20260927.json`.

Separately queued independent verification of all312 refined fits'936 node
marginals, expected168,282 node/sites and3,365,640 probabilities. It uses
direct matrix exponentials and fixed-root inside/outside messages, requires
the frozen completed fit audit, and compares every probability and site
likelihood at1e-8 tolerance after successful production. All model/aligner/
boundary/bound settings remain represented. No ancestral sequence selection
is made by either calculation.

Run `python scripts/audit_refined_domain_posteriors.py`; plan/launch metadata
use `refined_domain_posterior_audit_`, date `20260927`. Output:
`results/ancestral/refined-domain-posterior-audit-20260927-v1/`. One CPU,4GiB
RAM,no swap/GPU; estimate0.1–4 hours after production,0.2GiB output.
Refined probabilities and their verification remain pending.

## All posterior audits and optimization comparisons completed

The independent evaluator checked all10,096,920 alternate-start probabilities
and all3,365,640 refined probabilities. Maximum absolute probability errors
were4.40484e-12 and4.40528e-12, respectively; site-likelihood errors were below
2.47e-10. Both auditors terminated successfully and full artifact hashes
were checked. Completion records use prefixes `refitted_domain_posterior_audit`,
`refined_domain_ancestors` and `refined_domain_posterior_audit`, date `20260927`.

Compared the optimization stages at identical aligned sites and mapped nodes,
retaining925,551 site comparisons across1,716 fit pairs. Saved site records
independently reproduce the aggregate counts and probability distances.

| Comparison | Node/site comparisons | Different MAP residues | Maximum total variation |
| --- | ---: | ---: | ---: |
| Alternate starts versus baseline | 504,846 | 220 | 0.83055 |
| Refined versus baseline | 168,282 | 38 | 0.23188 |
| Refined versus selected start | 168,282 | 0 | 0.01113 |
| Refined lower versus original gamma bound | 84,141 | 0 | 0.00017584 |

No conflicting MAP pair in these comparisons has both probabilities>=0.90.
These are dependent comparison counts, not unique evolutionary events.
Unchanged MAP residues can still have different probability distributions.
Poorer alternate-start fits are retained as diagnostics, not assigned equal
biological weight. Refinement and bound stability do not resolve model
adequacy, alignment correspondence, ancestral node support or indel history.

Reproduce `python scripts/compare_ancestral_optimization_probabilities.py`;
plan/launch/completion metadata use `ancestral_optimization_probability_comparison`,
date `20260927`. Output:
`results/ancestral/optimization-probability-comparisons-20260927-v1/`.
One CPU,4GiB memory limit,noGPU. Both independent audits were required before
comparison; their receipt hashes are recorded with the output.

## Indel reconstruction implementation started

Amino-acid posteriors treat gaps as unknown observations and cannot establish
whether an ancestral residue existed. Moving toward sequence ensembles
therefore requires explicit gap/indel inference and uncertainty. Retrieved
the published FastML3.11 source archive from the authors and started a local
build of its `indelCoder` and `gainLoss` components. Archive SHA256:
`2d1ec87116eae1163177631ad39083db01cdd78ffb8d119bebbdbb6a43e99451`.
A full source-file checksum manifest, compiler version and build log are
recorded under `results/software/fastml-indel-build-20260927-v1/`.

The wrapper selects simple indel coding (`SIC`) and by default includes
leading/trailing gaps. Its coder also supports treating terminal gaps as
unknown; that distinction matters for incomplete proteins and extracted
domain boundaries. The source treats X as unknown for gap coding, distinct
from a known residue identity; this behavior must be checked against our
three whole-protein X positions before interpreting any indel call. We will
retain rather than silently convert ambiguity. Overlapping/nested gap
semantics, encoded coordinate mapping, and the advertised200-sequence limit
must be validated before full project inference. No indel result is claimed
from downloading or compiling software.

Build: `python scripts/build_fastml_indel_tools.py`; frozen plan/launch:
`metadata/fastml_indel_build_{plan,launch}_20260927.json`. Uses four CPUs,8GiB
RAM,no swap/GPU; estimate0.1–2 hours,2GiB output. GNU++11 mode is passed to
the legacy source build without modifying source. No remote sequence upload,
paid infrastructure, ancestral FASTA selection or GPU prediction occurs.
See `docs/bibliography.md` for primary methods/source references.

The FastML indel-tool build subsequently terminated successfully; both binary
hashes verified. Closure: `metadata/fastml_indel_build_completed_20260927.json`.
Functional/capacity validation and indel inference remain pending.

## Whole-protein models independently checked

All78 whole-protein fits and their full report/checkpoint audit terminated
successfully. Independent likelihood replay passed throughout, maximum
absolute error3.46471e-5 log units, including correct X/unknown-state handling.
All artifact hashes were checked. Sixty-three fits have edges shorter than
1e-5; warnings remain recorded and are not synonymous with failed fits.
Completion: `metadata/ancestral_whole_model_audit_completed_20260927.json`.
Whole-protein alternate-start stability, ancestral node mapping and
probabilities remain unfinished.

## Complete indel coding underway

Started SIC coding of all78 whole/domain alignments under both explicit
terminal-gap policies (156 encodings). Synthetic tests cover exact, nested,
partially overlapping, terminal and unknown-flanked gap characters before
production. A separate Python implementation independently reconstructs every
character interval and every0/1/? matrix cell from the source alignment.
Coordinates are exported as1-based inclusive start/end with original
0-based half-open tool output verified.

The622-protein whole MAFFT alignment has completed both policies with every
cell verified (1,516 characters with terminal gaps included;1,354 with terminal
gaps unknown). This demonstrates local indelCoder capacity for that full case;
it does not establish gainLoss inference capacity or accuracy. The complete
156-encoding run remains active. Here1 means an exact coded gap,0 means the
absence of that exact gap under SIC, and? denotes a containing gap or qualifying
unknown span/flank. A partial overlap may be0; this is not a direct binary
residue-presence matrix. Dependent overlapping characters cannot simply be
counted as independent insertion/deletion events.

Run `python scripts/encode_ancestral_indels.py`; plan/launch metadata use
`ancestral_indel_coding_`, date `20260927`. Output:
`results/ancestral/indel-coding-20260927-v1/`. One CPU,4GiB RAM,no swap/GPU;
estimate0.1–4 hours,2GiB output. Ancestral gap histories and final sequences
remain pending.

All156 indel encodings subsequently completed:12,957 character records and
4,042,247 coded cells independently verified, maximum622 proteins. All156
job receipts and artifacts checked after successful termination. Closure:
`metadata/ancestral_indel_coding_completed_20260927.json`. These are encoded
gap characters, not inferred historical indel events or ancestral sequences.

## Refined whole-protein probability propagation (September 27)

The 156 whole-protein refinements retain all 26 alignments, three substitution
models and two gamma-shape bounds. A queued propagation stage waits for the
independent refinement auditor to finish successfully, verifies its exact
process identity while waiting, and checks the final fit/parameter artifacts.
It then computes all 20 amino-acid probabilities at the three identifiable
candidate vertices for every fit (468 node/model combinations), preserving
all input proteins. X and gaps remain unknown observations; the degree-two
root position remains unidentifiable under these unrooted fitted models.

The stage uses the previously independently checked pruning implementation,
replays each fitted likelihood while computing probabilities, and retains full
probability arrays. Its outputs will still require an independent posterior
readback before downstream interpretation. No refined probability result is
claimed while the producer is waiting. The original estimates remain separate.

Reproducibility: `scripts/infer_refined_whole_ancestors.py`,
`metadata/refined_whole_ancestor_plan_20260927.json` and the matching launch
record. Resource allowance: one CPU, 8 GiB RAM, no swap, 4 GiB output and
0.5–12 hours after the refinement audit; no paid infrastructure or GPU work.

The matching independent posterior auditor is also queued. It checks every
probability for all 156 fits and 468 candidate vertices using fixed-root
inside/outside messages and direct matrix exponentials, rather than the
producer's rerooted, pattern-compressed spectral calculation. The analytic
two-internal-node enumeration check passed before launch. Complete array
shapes, finite values, normalization and saved site likelihoods are checked.
Tolerances remain 1e-8 for probabilities and 1e-7 for site log likelihoods.
The auditor waits for verified successful producer termination and checks
its artifact hashes; the full-data comparison has not yet run. Resource
allowance is one CPU, 8 GiB RAM, no swap, 1 GiB output and 0.5–12 hours after
the producer. See `scripts/audit_refined_whole_posteriors.py` and its
`metadata/refined_whole_posterior_audit_*_20260927.json` records.

All 468 whole-protein multistart fits and their independent likelihood audit
are now complete. Maximum replay error was 3.63e-5 log-likelihood units.
Among 156 model/bound groups, 48 had across-start differences above 0.001
(maximum 0.138223); 36 improved over the baseline by more than 0.001
(maximum 0.014044). Best-start refinements are running; these numerical
results do not establish global convergence or biological significance.

A pre-execution mapping check found that selected fitted trees suppress the
original degree-two root. The v1 probability producer and auditor were stopped
while still waiting, before probability output began. Separate v2 scripts,
plans and launch records supersede those two waiting stages. The v2 producer
uses the pinned original labeled trees to identify candidate signatures and
then maps them to the fitted tree. All 156 refined job inputs were checked:
all three identifiable signatures are retained and the root remains absent,
as required. Existing refinement jobs and earlier result collections are
unchanged. The v2 probability producer and independent auditor are verified
live and waiting for their required predecessors.

A full-coordinate comparison is queued behind the independent refined
posterior audit. It compares all 156 refined fits against their original
baselines and all 78 pairs of gamma-shape bounds, retaining 702 node
comparisons and 653,157 node/site comparisons. Each pair must share the
alignment hash and ancestor partition signature. It reports changes in
most probable amino acids, opposing calls supported at least 0.9 in both
fits, and total variation across all 20 amino-acid probabilities. These are
dependent optimization diagnostics, not biological replicates. Intermediate
alternate-start probabilities are not part of this comparison.

The script is `scripts/compare_whole_optimization_probabilities.py`; matching
plan/launch records use `whole_optimization_probability_comparison`. It has
one CPU, 4 GiB RAM, no swap, 2 GiB output allowance and a 0.1–2 hour estimate
after its required audit finishes. The live process is waiting; no comparison
result is claimed yet.

All 156 best-start whole-protein refinements and their independent likelihood
audit are now terminal with successful exit status. Every bound output
artifact was rechecked. Maximum independent likelihood replay error is
3.5902e-5. Relative to selected starts, the largest gain is 7.7e-6 and the
smallest change is approximately -1.0e-8; no absolute change exceeds 0.001.
These results support stability of this additional refinement step, without
proving a global optimum or model adequacy. The v2 ancestral probability
producer has begun writing results; complete posterior validation and
site-by-site comparisons remain pending. The completed audit is recorded in
`metadata/ancestral_whole_refinement_audit_completed_20260927.json`.

Full alternate-start probability propagation is now running for all 468
audited whole-protein fits, including poorer starts. It retains 1,404
identifiable ancestor/model combinations and all 20 amino-acid probabilities.
This permits assessment of probability sensitivity across optimizer outcomes
beyond the original-versus-refined comparison. It preserves all proteins,
alignment columns, models, starts and bounds; optimizer outcomes are not
independent biological replicates or weights for ancestral ensembles.

The producer and independent full-column posterior auditor use separate
`alternate_whole_ancestor` and `alternate_whole_posterior_audit` plan/launch
records. Producer allowance: one CPU, 8 GiB RAM, no swap, 6 GiB output and
1–24 hours. Audit allowance: one CPU, 8 GiB RAM, no swap, 1 GiB output and
1–24 hours after production. No paid resources or GPU predictions are used.
Production has started and the auditor is waiting; full completion and
optimization comparisons incorporating these alternate probabilities remain
pending. The existing refined comparison runs unchanged.

The expanded comparison is queued in a separate output collection. It covers
468 alternate-versus-baseline pairs, 156 refined-versus-baseline pairs, 156
refined-versus-selected-start pairs and 78 bound pairs: 858 fit pairs, 2,574
node comparisons and 2,394,909 node/site comparisons. Selected-start identities,
alignment hashes and expected site counts were checked across the complete
input plans before launch. Runtime checks also require identical ancestor
signatures and completed independent posterior audits for all collections.

The `compare_whole_optimization_probabilities_full.py` process is verified
live and waiting. Its allowance is one CPU, 8 GiB RAM, no swap, 4 GiB output
and 0.1–4 hours after the audits. All comparison categories retain full
coordinates; they do not supply independent replicates or ensemble weights.
Earlier comparison outputs remain separate.

All 8,708,760 refined whole-protein probabilities passed independent
fixed-root/direct-exponential readback (maximum absolute difference
1.672e-11; site log-likelihood difference 3.656e-10). Both producer and audit
are terminal with successful status and all artifacts verified.

The initial refined/baseline and bound comparison also completed: 14 changes
in the most probable amino acid among 435,438 refined-versus-baseline node/site
comparisons, with zero opposing calls having at least 0.9 support in both
fits. Maximum total variation was 0.003031. Across 217,719 bound comparisons,
there were zero most-probable-state changes (maximum total variation 0.0002595).
A separate readback of all underlying probability arrays reproduced the
site, state-change and high-support-disagreement counts. These results support
stability for these specific optimization contrasts; alternate-start comparisons,
model adequacy and insertion/deletion uncertainty remain separate. Evidence is
in `metadata/refined_whole_posterior_audit_completed_20260927.json` and
`metadata/whole_optimization_probability_comparison_completed_20260927.json`.

The [refinement sensitivity figure](figures/whole_refinement_sensitivity_20260927.pdf)
summarizes all 13 families. The 14 most-probable-state differences are
concentrated in OG0000336 (2), OG0001082 (8) and OG0002650 (4), counting
dependent model/bound/node comparisons rather than independent events. PNG,
PDF and editable SVG copies are provided. The plot script verifies its
complete 702-row input table and aggregate totals. Its final rendering was
visually checked; artifact hashes and scope are recorded in
`metadata/whole_refinement_sensitivity_figure_completed_20260927.json`.

The alternate-start producer has now finished all 468 fits, writing 1,404
candidate-node arrays with 26,126,280 probabilities. Complete artifact hashes,
array shapes, finite values, ranges and normalization were rechecked after
successful producer termination. This is recorded in
`metadata/alternate_whole_probabilities_produced_20260927.json`. The independent
probability replay is verified live; it and the full comparison remain pending.

All 26,126,280 alternate-start probabilities passed independent readback
(maximum absolute difference 1.674e-11). The full comparison finished all
858 fit pairs, 2,574 node comparisons and 2,394,909 node/site comparisons.
A separate readback of the underlying arrays reproduced every aggregate
count and total-variation summary. Alternate-versus-baseline comparisons
have 128 most-probable-state changes across 1,306,314 dependent node/sites,
with no opposing calls supported at least 0.9 in both fits. Refined versus
selected-start comparisons have zero state changes across 435,438 node/sites;
maximum total variation is 0.0002352.

The maximum alternate-versus-baseline total variation is 0.298073, at
OG0000972 whole-FAMSA WAG, start 2, alpha minimum 0.02, candidate level 2,
alignment column 1650. Thus low state-change counts do not establish uniformly
small probability differences. Largest site shifts are retained for follow-up.
These comparisons assess optimizer sensitivity conditional on the models and
alignments; they do not resolve joint indel histories or qualify final sequences.
Evidence: `metadata/alternate_whole_posterior_audit_completed_20260927.json`
and `metadata/whole_optimization_probability_comparison_full_completed_20260927.json`.


### Descriptive investigation of extreme whole-protein support shifts

The reproducible `scripts/diagnose_whole_ancestral_extremes.py` joins the ten
retained per-fit maxima to the original alignment, mapped descendant sets,
audited likelihoods and baseline/alternate/refined probability arrays. All
arrays are checksum-checked against producer receipts; the full probability
audits remain the numerical validation. Output and input hashes are recorded
in `metadata/whole_ancestral_extreme_diagnostics_completed_20260927.json`.
These ten records are not the ten largest individual sites globally.

The largest shift is OG0000972, FAMSA/WAG, candidate level 2, alignment column
1650 (one-based). Of 622 tips, 595 contain canonical amino acids and 27 have
gaps at this column; there are no unknown X characters there. The candidate
has 621 descendants; its one outside tip carries proline. This documents
coverage and phylogenetic context, not alignment correctness or reliable
root placement.

Both baseline and alternate fits favor asparagine. Its conditional probability
changes from 0.602051 to 0.900123, while the alternate fit has a lower log
likelihood by 0.126265. The selected refined fit has probability 0.602082 and
is only 0.00003576 total variation from the baseline at this site. Gamma
shape and counts of near-zero edges also differ, but these observations do
not isolate the cause of the probability shift. Higher apparent confidence
in the poorer fit is not stronger historical evidence. All alternative
results remain retained; no likelihood-based ensemble weights or biological
substitution counts are inferred from these diagnostics.


### Refined whole-protein versus domain conditional probabilities

All 104 alignment-coordinate comparisons now connect to the audited refined
probability arrays: 624 fit pairs across three substitution models and two
parameter bounds, with 1,872 non-root candidate comparisons. Source-node
identity is checked and whole-tree descendant sets must restrict exactly to
the domain-tree sets. Only exact projected residue-coordinate signatures are
compared; all unmatched counts remain in each node summary and the complete
coordinate union remains in the source artifact. This includes both matching
and differing alignment methods and both Pfam boundary definitions.

Of 269,154 matched node/site comparisons, 10,320 change their most probable
amino acid. In 104 comparisons, different residues each receive at least 0.9
support in the corresponding fit. The same-alignment-method subset contains
137,988 comparisons, 5,528 state disagreements and 52 opposing high-support
calls; cross-method comparisons contain 131,166, 4,792 and 52, respectively.
A separate full serialized-table readback reproduces every site probability,
state call, total variation and all 1,872 node summaries from saved arrays.

These are dependent comparisons, not unique substitutions or independent
replicates. Differences jointly reflect membership, flanking sequence,
alignment and fitted model parameters; they do not isolate a causal effect of
protein context. Even with the same alignment method, separate alignments and
different protein sets remain. The results are conditional on residue presence
and assumed local trees. They expose uncertainty beyond optimizer sensitivity
and must accompany any future ancestral structural interpretation.

Reproducibility: `scripts/compare_refined_whole_domain_probabilities.py` and
`scripts/readback_refined_whole_domain_probabilities.py`; completed results
and provenance are recorded in
`metadata/refined_whole_domain_probability_comparison_completed_20260927.json`
and `metadata/refined_whole_domain_probability_readback_completed_20260927.json`.


### Localization and dependence of opposing context calls

The 104 opposing high-support comparisons reduce to **three exact extant
residue-coordinate signatures across five source-tree nodes**, not 104
independent sites. The grouped results retain every available comparison at
each signature, including those without opposing high support:

| Family | Source node | Whole → domain state | Opposing / available comparisons |
| --- | --- | --- | --- |
| OG0001203 | n36 | H → T | 8 / 48 |
| OG0001203 | n37 | H → T | 28 / 48 |
| OG0002650 | n619 | Y → M | 20 / 24 |
| OG0002650 | n620 | A → V | 24 / 24 |
| OG0002650 | n622 | A → V | 24 / 24 |

Arrows indicate a change in the estimate between analysis contexts, not an
evolutionary substitution. The two H/T nodes share one coordinate signature;
the two A/V nodes share another. All seven observed residues in the retained
A/V coordinate set are valine, but the whole-protein analysis also includes
other proteins and sequence context. This does not establish the cause of
the conflicting alanine estimate. The H/T and Y/M signatures contain 14 and
13 observed residues, respectively. Missing/gapped proteins are not represented
by those coordinate counts.

`scripts/localize_ancestral_context_conflicts.py` verifies each signature
against original protein positions, stores every gene/position/residue, and
retains all model/bound/method comparisons for these groups. The output
`results/ancestral/ancestral-context-conflict-localization-20260927-v1/conflict_dossiers.json`
is checksum-bound by
`metadata/ancestral_context_conflict_localization_completed_20260927.json`.
Grouping removes duplicate settings from descriptive site counts, but shared
nodes, data and parameters still preclude treating the groups as independent.
These sites require contextual review before mechanistic interpretation.
