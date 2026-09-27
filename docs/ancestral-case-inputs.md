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
