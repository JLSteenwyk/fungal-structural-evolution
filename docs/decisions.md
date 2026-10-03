# Decisions and unresolved questions

- User-directed scope: approximately 500 unique fungal species and 25 additional non-fungal outgroups, no separate pilot.
- Repository name: fungal-structural-evolution. User confirmed public repository JLSteenwyk/fungal-structural-evolution; origin configured and remote verified empty before initial push.
- NCBI RefSeq and GenBank fungal assembly catalogs are initial discovery sources, not sufficient coverage by themselves. Augment sparse lineages with published datasets and repositories.
- Fungal circumscription is unresolved. Szánthó et al. (2025) use a narrower osmotrophic Fungi definition, treating Aphelida, Rozellida and Microsporidia as close relatives. Record NCBI classification separately from study role; do not silently count boundary taxa as non-fungal outgroups. Assess alternative definitions before freezing roles.
- Shared host: do not infer exclusive resource allocation from physical capacity. Use conservative metadata jobs initially; no GPU jobs launched.
- Existing `gh` executable is not the GitHub CLI (it rejects `auth status`); a GitHub configuration file exists. Do not expose credentials. Establish an appropriate client before remote creation.

Open: outgroup availability and boundary definition; trait sources/replicate transitions; predictor/database availability; actual reusable structure coverage; sustained shared-machine allocation.

- Source audit: Szánthó et al. Supplementary Table 1 assigns the Rozella allomycis JGI URL to Parvularia atlantis as well. Retain the original source row but do not download that URL as Parvularia. NCBI independently lists Parvularia assembly GCA_943704415.1, annotation availability pending.
- All 35 originally missing taxonomy IDs are in NCBI delnodes.dmp; none are merged IDs. Keep their records and flag deleted; do not infer replacements from similar names.
- Availability validation exposed trailing slashes in NCBI ftp_path. The initial malformed URL checks are invalid for taxon exclusion. Corrected URL construction strips trailing slashes; check cache is keyed by full URL, so corrected URLs are checked independently. HEAD availability remains distinct from validated FASTA content.

- Initial 500-genome QC draft uses the broad NCBI Fungi circumscription, including Microsporidia, Rozellomycota and Aphelidiomycota as ingroup. These cannot also count toward the 25 non-fungal outgroups. A narrower osmotrophic clade can be tested with explicit re-rooting/pruning and separate counts.
- Draft selection prioritizes unrepresented orders/families/genera and up to 12 annotated species from smaller phyla, plus focal genus sampling. This is a transparent taxonomic proxy, not branch-length-based phylogenetic optimization. Focal genera provide potential contrasts, not verified ecological transitions. Sanchytriomycota requires an external annotation source and is a known coverage gap.
- Full-scale proteome QC began directly for the 500-species draft. QC failures may require replacement; final sampling is not frozen.

## Taxon identity and hybrid sensitivity

Name-based review flags 22 of the 526 working entries (21 incompletely identified fungal labels and one explicit Saccharomyces hybrid). Distinct taxon IDs must not be equated with established distinct species. Preserve current data and exploratory runs, but resolve names/assembly identity and assess exclusion of flagged taxa before final species-level inference. In particular, hybrid ancestry can violate a single bifurcating species-tree model. The hybrid is involved in several large exploratory global RMSDs; this is a review priority, not evidence that hybridization caused structural change. See `metadata/taxon_label_review.tsv` and reproduce with `python scripts/audit_taxon_labels.py`.

## Paired structural branch inference

Use identical observed cells and sequence-based residue correspondence for AA and native coordinate-derived 3Di fits. Preserve invariant sites and missing-data reasons. Fit structural branch lengths on the sequence marker topology using the published AF empirical model, with frequency/training-source sensitivities. Preserve unrooted splits and full taxon-set identity; do not compare edges from different taxon subsets as though they are the same branch. Treat near-zero lengths and rare states as uncertainty/model concerns, and withhold branch ratios and acceleration rankings until those concerns are addressed.

## Dependence and conditional resampling

Resample identical columns in AA and 3Di alignments across all taxa. Use 200 draws each at block lengths 1, 10 and 30 as a conditional sensitivity analysis, retaining invariant sites and recording unestimable draws rather than silently redrawing. Do not interpret percentile intervals as overall calibrated confidence: the audited inputs contain many spatially linked feature pairs farther apart than the chosen blocks. Preserve paired sampling covariance for later error-aware coupling analysis; it is not an evolutionary correlation.

## Source-stratified structural snapshots

Use an explicitly selected provider and prediction tool for the primary evolutionary mapping; rank candidates only within that source and leave unavailable sources missing. Retain all alternative models for same-sequence source comparisons. Higher pLDDT across pipelines is not evidence of greater accuracy and must not silently switch the production source. The first available GDM/ColabFold pair is a control with limited scope, not a calibration of pipeline effects across fungi.

## Cloud availability and prediction strategy

The user confirmed that no existing Google Cloud account is available and asked whether ESMFold would be faster. Cloud authentication is therefore not an available acceleration route for the previously attempted proteome archives. Continue the already authorized individual-model downloads and local ESMFold execution; no account creation or paid cloud provisioning is planned.

At this check, the active ESMFold run had approximately 3,400 completed per-sequence prediction receipts, with median measured inference approximately 1.69 seconds and mean 2.43 seconds for the completed subset. The current production queue is limited to canonical proteins of at most 512 residues; this is a local execution limit, not an ESMFold hard limit. Completed short-sequence timings do not forecast long proteins or the entire multi-million-protein atlas. Individual AlphaFold acquisition and PAE acquisition are also progressing, so replacing already available models with new predictions would duplicate work. Continue complementary retrieval and prediction, with explicit predictor strata and same-sequence controls to avoid confounding clade effects with prediction methods. EBI also documents public FTP bulk subsets; the earlier authentication failure concerned the specific Google-hosted archives attempted, not all AlphaFold access.

## Integrating disjoint ESMFold acquisition cohorts

Combine completed cohorts only after checking identical inference settings apart
from input receipts and physical GPU identifiers, disjoint model/sequence IDs,
and unchanged coordinate provenance. Remap the complete union against the fixed
full-panel sequence matrix. Before sharing qualified encodings, require every
combined model-provenance, marker-link and residue-mapping row to equal the
disjoint source union. Retain original encoding paths/checksums and native,
coordinate, PAE and qualification receipts; the integration receipt explicitly
describes a derived union, not a new native audit. Reassess marker eligibility
using the combined full-design grid. Matching settings do not establish absence
of batch effects, ecological replication, or sequence-independent evidence.

## Ancestral uncertainty and numerical checks — September 27

Retain the complete family inputs, alternative fits and failed capacity
outcomes when advancing ancestral work. Independent amino-acid probability
agreement verifies numerical implementation conditional on the specified
model; it does not qualify insertion/deletion histories or final sequences.
The observed incompatibility of separate gap characters and largest-family
Historian branch-floor sensitivity require joint-history uncertainty to remain
an explicit unfinished component. Neither compatibility conditioning nor the
highest-scoring selected reconstruction substitutes for that component.

Exact BAli-Phy input equivalences retain every original label, model and tree
choice. Computation reuse is allowed only for identical effective inputs and
seed/settings; aliases are not independent posterior chains. Completed
initialization and short diagnostic chains do not establish mixing. Future
ancestral structure predictions must use qualified alternatives and the
current GPU authorization. See [ancestral evidence](ancestral-case-inputs.md),
[Historian](historian-method-assessment-20260927.md) and
[BAli-Phy](baliphy-method-assessment-20260927.md).

## Sequence-first reference sensitivity — September 30

Preserve the full available-reference design and add a sequence-first choice
sensitivity for every modeled duplicate target. Choose nearest nonfocal sister
genes and a fixed lexical representative before checking structural availability;
retain every tied gene/model assignment and missing lexical representative.
Do not silently substitute a modeled tie or farther gene, promote ineligible
parent contexts, switch predictors, or rematch fixed background controls.
Keep gene-choice availability sensitivity distinct from sequence-locked residue
correspondence. Test resulting structural contrasts only after complete input/
result reuse, geometry, coverage/confidence/PAE/domain and orthology checks.
The completed native check finds a farther modeled reference in approximately
25.9% of provisionally available targets; this requires missingness/reference-
choice treatment but does not itself establish structural bias or asymmetry.


### September 30: full reference measurement partition preserves pending reuse

Keep the complete expanded reference ledger while measuring every pair without
matching catalog sources. Existing-source matches are pending actual
coordinate/input/mask/residue/executable/settings/checkpoint and numerical
qualification; a catalog match never authorizes importing results or clearing
flags. No structural outcome selects computation or fixed background matches.
Check both masks/orders of all new pairs and retain errors and short inputs.
The final reference/background union must also account for overlapping pairs
before sharing verified results. CPU native measurements remain authorized;
GPU predictions stay paused.


### September 30: reuse compatible results with all original exclusions intact

Qualify actual old/current raw and masked coordinate bytes, original residue
positions, sequence/length, executable/options, directed checkpoints and
complete numerical/geometry provenance. Prefer completed expanded primary
results, otherwise old reference results, without choosing on response or
usability. Compatible excluded/short/error results remain excluded; never
replace them with zero distances or clear flags using a better-looking
alternative source. Different source directions map by ordered model/version
endpoints. Archive large raw-byte/checkpoint provenance outside Git with
versioned locations/checksums. Independent full readback and journal closure
precede downstream result union and inference.

### September 30: query native reference orthology before structural inference

Retain all source target contexts under both native guides and both nearest
reference designs, including unmodeled references, all ties and parent
exclusions. Deduplicate physical gene queries only; preserve every logical
context/design/tie/duplicate-side link. Bind source-line ordinals through full
native identity audits and resolved-tree readbacks. Compare both duplicate
memberships under both native guides; absence or disagreement remains explicit.
Native coorthology does not override parent exclusions or prove biological
orthology. Missing reference genes are unqueried contexts, not negative
membership. Full independent source/link reconstruction and exact process
journals precede completion. Preserve failed attempts and use fresh versioned
scripts/plans/outputs for recovery.

### September 30: complete the whole reference measurement grid before union

Require the complete new-native disposition grid, unchanged numerical/geometry
proof lineage and all four exact completion journals before joining new
measurements to qualified old results. Preserve every old/new state and
original native checkpoint/metric/flag/source order. Map current directions
by ordered model/version endpoints, without editing original source records.
Independent full field reconstruction and both union process journals precede
downstream qualification. Keep large byte indexes outside Git with versioned
locators. A numerical union does not qualify biological coverage, residue
correspondence, prediction uncertainty, orthology or calibrated effects.

### September 30: require both orders and full-protein reference coverage

Apply the existing six 30/50-residue by 50/70/90-percent coverage thresholds
to the entire closed reference measurement union. Preserve current ledger
endpoint directions separately from original checkpoint orders. Keep missing,
failed and numerically excluded states with blank usable metrics and explicit
raw aligned lengths/flags. Both native orders must pass; do not choose or
average favorable measurements. Confidence masking never changes the full-
protein denominator. Label retained-input coverage separately from original
coverage. Independently reconstruct every source/order/metric/blank/exclusion
field and exact threshold decision; use decimal ceilings as a check on integer
fraction screening. Require full input preflight and source/journal-closed
union before production, and full independent readback plus both original
process journals before accepting numerical coverage completion. Pair coverage
does not establish common residue correspondence, biological reference
orthology, independent events or calibrated duplication asymmetry.

### September 30: fix complete phylogenetic context-to-measurement roles

Join every closed native target context to the frozen duplicate gene/model/
version queue and complete reference pair/availability-side ledger before
projecting coverage or effects. Match genes by identity rather than positional
a/b labels. Preserve both designs, every tie and original lexical choice,
all native source fields and parent exclusions. Classify missing structures,
identical-model comparisons, pairs in the physical measurement design and
unmeasured distinct pairs separately. A pair measured for a different event
does not make an excluded parent context eligible. Require complete independent
SQL field/role reconstruction, exhaustion of both full context and original
availability-link universes, source hashes and both original process journals.
Work-count exports label context counts separately from logical side-link
counts; do not pool counting units, guides, reference designs or ties.

### September 30: retain the full reference-context coverage denominator

Project the source/journal-closed full original-protein coverage matrix onto
every fixed native context/design/tie and both duplicate sides. Require both
reference-side comparisons to pass under both native input orders; carry
all unavailable/error/numerical flags and explicit null unmeasured order
fields. Retain missing genes/models, identical-model and outside-design work
states. Parent eligibility is mandatory for each contextual policy even when
a physical pair is measured in another event. Report five separate lexical
coverage, lexical own-guide assignment, lexical both-guide assignment, any-tie
and all-tie both-guide diagnostics. Any/all sensitivities do not change the
original lexical reference; empty all-tie sets are excluded. Preserve full
source, parent, guide and model availability denominators. Independently
reconstruct each exported field and derive all contextual policies using SQL
boolean/count aggregates. Full producer/reader journal closure precedes
completion. Reference-side coverage does not qualify the primary duplicate
AB comparison, an independent common-residue triad, prediction accuracy,
biological orthology, accepted phylogeny or a calibrated asymmetry effect.

### September 30: require all three duplicate/reference edges and explicit identity screens

Carry primary AB work alongside AR/BR for every original context/reference tie.
Preserve the complete closed source record, gene-based A/B orientation, parent
exclusions, lexical choices and modeled/unmodeled reference records. Deduplicate
only ordered physical A/B/reference model/version triples and retain every
logical occurrence. Require both three distinct versioned models and three
distinct model IDs for source-only correspondence work; two versions of one
model do not count as independent proteins. Explicitly retain identical-model
and outside-design edge states with null endpoint/order fields. Source readiness
also requires the original eligible parent, queued primary AB and all three
edges in their designated frozen catalogs. Record native own/both-guide
assignment sensitivities separately without changing lexical choices or ties.

Require full independent SQL reconstruction, complete physical/logical exhaustion
and both exact original producer/reader completion journals before publishing
counts. Label source contexts, logical reference records, ordered physical
triples and future mask/order states separately. Source readiness does not
qualify native coverage, common residue correspondence, biological orthology,
prediction error or asymmetry. The full future mapping design retains both
masks and all eight edge-order combinations, with reference-common and
AB-cycle-consistent residue sets assessed separately before physical fitting.


### October 1: close original correspondence jobs before geometric fitting

Interactive access failure did not stop the original CPU jobs. Use their closed
full-scope outputs rather than relaunching them. Preserve preflight-pinned
prototype sources and production v2 plans unchanged. Record full-input/grid,
independent reconstruction, artifact hashes and exact original journals before
claiming completion. Common-reference and cycle-consistent cores remain separate;
failed sources and short cores remain explicit, and favorable native orders are
not selected. Subsequent three-edge fits must use identical triples and recheck
proper-rotation geometry, inherited exclusions and original coverage. Raw signed
RMSD differences are descriptive contrasts, not directional evolutionary rates.


### October 1: separate shared-core and inherited three-edge qualification

Use exact same-residue rigid fits for all three distances under both mapped
cores and every native order. Reuse computations only for identical physical
triple/mask/residue sets; preserve every output disposition and logical source
link. Keep inherited both-order three-edge coverage separate from shared-core
geometry/coverage and require both for combined screening. Blank uncomputed
metrics distinguish failed/short states from computed degenerate geometry.
Require the exhaustive closed work preflight before production, then an
independent quaternion reader and both exact original completion journals.
The completed fourth PMSF audit does not establish an accepted species framework;
retain its gap warnings and require all-four topology/consensus and broader
phylogenetic sensitivities.


### October 1: preserve tree-specific support and sparse-taxon uncertainty

Keep all four alignment/guide combinations and both ML/consensus views. Report
ML SH-aLRT80/UFB95 criteria separately from consensus UFB95, leaving unavailable
SH-aLRT blank. Exhaust all split/missing/conflict/witness records and require
independent raw-tree reconstruction and original journal closure. Compare full
character coverage on both matrices and retain every nearest declared-role
boundary without selecting a root or changing taxon roles. The sparse
Amoeboaphelidium protococcarum placement motivates taxon/marker sensitivities;
its observation does not prove causal missing-data bias. Pruned-topology
checks, if used next, cannot substitute for native sensitivity refits.

Same-residue numerical completion does not establish calibrated asymmetry.
Preserve all order/mask/core alternatives and perform complete logical-context
linkage and robustness qualification before comparative inference. Full
additional-background input closure precedes any actual alignment reuse or
new native measurement acceptance. Keep all eight scientific aims open.


### October 1: qualify complete physical grids before original-context policies

Retain all eight native alignment orders for each ordered model triple, mask
and common-core definition. Publish all-order and any-order coverage separately;
the latter cannot qualify a strict all-order cohort. Numeric ranges retain null
uncomputed values and original A/B sign. Both-mask qualification intersects
the same order bits; it does not average distances across different residue
sets. Join physical results to every unchanged original context/tie with its
original parent/model-identity and native-guide gates. A measured overlap cannot
promote an excluded context, and a missing lexical model cannot be replaced
after observing geometry. Empty all-tie sets remain excluded.

The full independently closed contextual design supplies future comparative
cohorts, not an accepted biological orthology assignment or calibrated effect.
Guide/design/core policies overlap and their counts must not be pooled as
independent events. Preserve all source choices in models of availability,
missingness, family/taxon effects and phylogenetic uncertainty.


### October 1: retain legacy-only active models and qualify full background inputs

The modern three input collections omit 224 active background models. Include
the complete older reference supplementary collection after full written-input
readback rather than silently dropping these models. Require every full
semantic overlap and actual written byte hash to agree, and preserve all source
origins. Use declared collection order independently of readiness or geometry;
restore original residue positions through exact source-row hashes before
numeric/geometry work. Keep all 146,172 physical background pairs and all
305,434 active models. The 74,712 genuinely new native pairs come from full
catalog screening, not the inventory's older additional-pair flag. All 71,460
matching catalog pairs remain pending real input/checkpoint/numeric qualification.

Queue eight-CPU native work only behind full independent input-union closure.
Use separate full vectorized/quaternion reconstruction after fused numeric/
geometry assessment; retain all short/failure/parse/timeout/RMSD/rank states.
Require three original journals and full artifact/source hashes for new
measurement completion. Software fixtures are synthetic implementation checks,
not a pilot or evidence that production data passed. Preserve failed legacy
completion v1; corrected v2 adds explicit evidence pins and reuses the passed
reader. Keep GPU inference paused, existing CPU jobs untouched and every
scientific aim open.


### October 1: qualify all actual old background states before result union

Keep all 146,172 pairs and both masks/directions in the reuse qualification
ledger, including 74,712 new pairs as pending. Use old reference first and then
old background according to an outcome-independent fixed order; preserve all
matching alternatives and source endpoint order. An input mismatch requires
new measurement, not an alternative source chosen after observing a fit.
Require actual raw/PDB/checkpoint hashes, full original residue/sequence/length/
mask/status identity, native settings and complete original proof lineage.
Reconstruct all reused aligned numeric values on current byte-identical inputs,
then independently audit all exports, quaternion RMSDs and rotation geometry.
Preserve numerical flags, unavailable and error states. Final reuse closure
requires two new and eight original native/numeric/geometry/reader journals.
The full new/reused result union and evolutionary controls remain separate
requirements; source or journal matches alone do not qualify a measurement.


### October 1: preserve every old/new measurement state in the full background union

Require the full new-native three-journal and actual-reuse ten-journal closures
before merging all 584,688 states. Keep fixed old-reference/background source
preference, original source directions, model roles, raw numeric/geometry
records and every native/parse/timeout/unavailable/numerical exclusion beside
normalized fields. A mismatched reused input blocks full union until an
additional explicitly closed native workload exists. Independently reconstruct
all export fields through SQL source lookup and require both original union
journals/full hashes for completion. Do not infer coverage, biological effects
or phylogenetic eligibility from merge or numerical eligibility alone. The
15 false-export software checks and synthetic prior proof contracts qualify
implementation only. Queue under two CPU/64 GiB per stage, preserve all original
CPU jobs and leave GPU inference paused.


### October 1: apply identical full-protein screens before fixed-match attrition

Queue the complete background coverage stage after successful full measurement
union closure. Keep the same six thresholds as duplication targets, require both
native directions, preserve all original exclusions and use original full
protein lengths for masked coverage. Distinguish all 305,434 inventory models
from 292,326 background endpoints rather than equating these universes; retain
the broader inventory and test an unused model before production. Reconstruct
all rows independently with SQL/decimal ceilings and close with complete hashes
and both original journals. Physical target/background counts are descriptive;
full 4,250,692 selected/56,965,652 unmatched fixed decisions, target/control
attrition and post-screen balance remain required without rematching. Keep
all eight aims open and GPU inference paused.


## Source-specific matched coverage and journal revalidation — October 1

Use the dedicated full graph/covariate proof for graph and covariate identity,
then tie fixed matching to that covariate receipt; selected matching exports
do not alone bind the complete graph. Check duplication targets against their
original catalogue and target physical QC, backgrounds against their original
catalogue and current background physical inventory. Keep same-model nodes
explicitly excluded. Preserve failed launched versions; corrections use new
scripts/plans/units/output paths.

Freeze all original control choices across the entire scenario design before
coverage outcomes. Preserve original unmatched and failed selected records,
complete denominators and empty strata. Attrition exports do not replace
post-screen balance, dependence control or biological calibration.

Raw `journalctl -o json` byte hashes can vary with serialization even when
parsed messages are identical. Preserve original closure hashes and revalidate
exact original process identity plus recorded terminal completion fields;
store fresh raw/canonical diagnostics separately. Do not claim original raw
byte equivalence, accept collected-unit default success fields alone, or
reinterpret journal validation as a repeat scientific-input/geometry audit.


## Full post-screen baseline and reuse definitions — October 1

Keep three baselines distinct: all original modeled target records, all
metadata-matched targets in each fixed scenario, and all target-screen-eligible
records including unmatched targets. Paired log-distance summaries exclude
nonpositive distances explicitly without an epsilon; constant vectors have
exact zero variance and standardized quantities remain nonestimable.

Report control reuse both by original node and by canonical model/version
physical pair. Gene-labelled aliases of the same pair do not create additional
physical observations. Reciprocal reuse weights and Kish concentration are
descriptive mass/concentration diagnostics; they do not qualify an inferential
weight scheme, evolutionary sample size, independent events or phylogenetic
correction. Retain original matching and all source failures/unmatched records.


## Native taxon refits and conditional guide reuse — October 1

Preserve the full 526-entry baseline and all 25 original outgroups. Test four
explicit taxon policies by native supported refitting on both alignments and
both conditioning guide alignments. Preserve every retained original sequence,
site/marker mapping and all exclusion/lineage counts. The low-occupancy policy
retains 23 outgroups; detailed lineage-bin losses are sampling sensitivities,
not a preferred filtered design or diagnosed biological loss.

Reuse exact-matrix homogeneous guide files only as frozen conditioning inputs
with full input/config/build/artifact and report/tree/tip checks. Do not claim
new supported inference or retrospective original-process completion from
those files. Refit every mixture profile and supported sensitivity tree;
topology pruning alone is insufficient. Only one mixture fit at a time and
only after currently available memory reaches its recorded prelaunch threshold.
Preserve native/checkpoint process identities and complete successful records;
full collection/journal/biological qualification remains separate.


- October 1: require a second, algorithmically separate complete reader for
  all native taxon PMSF sensitivities. DendroPy/Decimal/raw-bootstrap checks
  reconstruct numerical outputs independently of the Bio.Phylo first audit.
  Software validation uses the complete existing baseline and explicit
  synthetic fixtures; it is not a smaller biological pilot. Preserve every
  policy's retained-role boundary, including absence, without assigning a
  root. Close production only after both original batch/reader journals and
  all output/source bindings pass. Cross-policy/baseline comparisons and
  biological framework qualification remain separate required stages.


- October 1: recompute baseline bootstrap support on each exact retained
  cohort before comparisons with native taxon refits. Deduplicate projected
  splits within each original replicate; do not sum source split frequencies
  or inherit SH-aLRT labels. Record collapsed path lengths/component counts
  as full-inference reference quantities, not native subset estimates.
  Require complete independent pruning and two original completion journals.
  Keep the actual 16 supported refits and subsequent biological qualification
  as required stages.


- October 1: retain every comparison/conflicting split and quartet across all
  32 retained-reference views, rather than only selected stable branches.
  Compare views on the same cohort and report projected UFB95 criteria with
  unavailable SH-aLRT. Require full projection closure before production and
  independently reconstruct RF/compatibility from raw pruned trees. Keep
  actual native-versus-reference comparisons and biological framework
  qualification as separate required work.


- October 2 UTC: compute native gCF across every one of the 250 audited marker
  trees and all eight completed candidate species views, retaining all
  1,046,000 branch–gene states and actual decisive denominators. Join original
  support by canonical split because native reorientation can move labels.
  Require a separate complete bipartition/NEXUS reader and two original
  completion journals. Preserve the first reader's inventory failure and
  implement corrections in new v2 scripts/plans/outputs. Keep unavailable
  factors, gene-estimation uncertainty and biological causes explicit;
  these measurements do not substitute for phylogenetic qualification.


- October 2 UTC: estimate the full gene-tree species-tree alternative using
  all 125 markers from each alignment under every closed taxon cohort and
  uncontracted/SH10/SH80 settings. Preserve missing-support dispositions;
  contract before pruning and keep exact-cutoff edges. Do not equate SH-aLRT
  with the literature's bootstrap contraction thresholds or call filtered
  trees posterior draws. Require full independently reconstructed input
  splits and both original completion journals before native inference;
  qualify quartet/local-posterior/coalescent-length estimates and compare
  all candidate trees before biological integration.

- October 2 UTC: independently recompute local and global quartet evidence
  for every coalescent candidate. Preserve unavailable and unresolved
  evidence before gene-specific normalization, check default-prior posterior
  and MAP coalescent lengths, and require both original completion journals
  after all 30 cases. Compare NNI alternatives as an unordered pair; do not
  infer their biological polarity from native labels. Full named-case proof
  validates the implementation and leaves the remaining batch and biological
  model/root/reconciliation qualification open.

- October 2 UTC: compare every coalescent candidate with all same-cohort
  original/projected concatenated references and all other coalescent
  alignment/support candidates. Preserve every conflict and support metric;
  do not equate local posterior, UFB and SH-aLRT, compare lengths across
  units, assign a root from a role split, or interpret overlapping split
  pairs as independent evolutionary events. Require complete native
  numerical closure, raw-tree comparison readback and both original
  completion journals before accepting the full comparison stage.

- October 2 UTC: display the full 315-comparison grid after actual full
  source closure, retaining every cohort and candidate/reference pair.
  Use a common normalized-RF scale across the two pages; blank triangular
  cells are not additional observations. Require a separate complete
  source-to-cell and PDF numeric-label reader, full hashes and both original
  completion journals, followed by visual inspection of both actual pages.
  Software fixture previews are explicitly labeled and remain outside the
  production results. Point-tree distances do not select a preferred model
  or establish biological causes, roots or calibrated structural effects.

- October 2 UTC: preserve the failed full quartet reader and every affected
  dependency attempt. Implement ASTRAL's exact effective-N substitution rule
  in new v3 readers, keeping raw fractional evidence, available genes and
  native-effective denominator distinct. Verify actual installed-jar boundary
  behavior and retain the strict 2e-8 comparison tolerance; no broad tolerance
  increase or native inference rerun. Use new v2 plans/output directories
  throughout the full audit/comparison/figure chain. Dependency handoffs
  require invocation-linked original completion/resource records rather than
  collected unit default fields. Retain exact failure exit records when a
  short failed invocation has no CPU summary; never infer successful completion
  from those missing resource records.

- October 2 UTC: publish actual full coalescent/reference figures only after
  all 30 native candidates pass numerical proof, all 315 comparisons pass
  independent raw-tree readback, and every figure cell/PDF value, full archive,
  original journal and both actual rendered pages are verified. These stages
  are now complete. Preserve the earlier failed attempts and the separate
  effective-N evidence; no inferred tree is rerun or silently altered.

- October 2 UTC: compare all actual native taxon-subset fits against every
  matched projected reference/coalescent candidate and each other, covering
  all 88 views/560 pairs after full native collection closure. Preserve native
  consensus multifurcations and report observed split counts/missing slots,
  maximum-resolved RF normalization and observed-split normalization separately.
  The latter is unavailable at a zero/zero denominator. Verify all exports
  independently from raw trees; require both original producer/reader journals
  and full hashes. Fixed source cohorts, every missing state and separate
  support/length semantics remain explicit. The full synthetic workflow passed
  with all 15 corrupted exports rejected; actual production remains dependent
  on the original 16 native fits. No GPU or paid resources are added.

- October 2 UTC: apply sequence-derived correspondence controls to every
  source-ready expanded physical triple, while retaining all original
  missing/excluded/unscheduled context and role records. Align full sequences
  before confidence filtering with both installed MAFFT/FAMSA under all six
  input orders. Do not select alignments from structural outcomes, substitute
  the previous 48 identical-domain cases for the full design, or count methods/
  orders as biological replicates. Require independent raw sequence/position
  reconstruction and both original journals/full hashes before downstream
  same-residue geometry under both masks. Preserve source-parent eligibility,
  original-protein coverage and explicit native failures. Sequence-derived
  correspondence is a sensitivity control; predicted coordinates still derive
  from sequences and biological homology/asymmetry require further evidence.

- October 2 UTC: require exhaustive independent work preflight before the full
  649,344 sequence-derived geometry grid. Use authoritative original confidence
  mask positions; PDB rounding must never create eligibility. Export rounded
  PDB confidence summaries separately from authoritative mask fractions. Fit
  all three distances on the identical residue triples, with coverage relative
  to full original protein lengths. Retain six core screens, inherited
  both-order pair gates and their intersection separately. The complete
  432,896-view source gate-constancy check permits reuse of those inherited
  gates within an ordered triple/mask; it does not qualify excluded parents or
  logical contexts. Preserve all native failure/short/rejected/degenerate
  dispositions. Full raw-MSA/quaternion/SQLite readback and original journals
  plus hashes must close before structural/phylogenetic interpretation. Resume
  only with source/plan-bound checkpoints reconstructed and rechecked; never
  duplicate a live producer. Full software tests passed and all six subsequent
  jobs are queued after the original native closure. Production completion and
  biological integration remain outstanding; GPU prediction stays paused.

- October 2 UTC: compare every sequence method/order with every structural
  core/order for the same original triple and mask. Retain all 10,389,504
  pair states and 216,448 joint groups alongside the 108,224 all-sequence-order
  groups. Preserve exact correspondence overlap, original-length coverage,
  statuses and exclusions, all 18 shared numeric differences and each source's
  six screening decisions. Empty-union overlap is unavailable. Strict numerical
  sign agreement is descriptive; different residue cores, masks and prediction
  sources cannot establish evolutionary polarity or extra independent events.
  Require independent raw-MSA/sorted-tuple/SQL reconstruction, deterministic
  checkpoint/export checks and both original journals/full hashes. Source
  acceptance requires both full geometry closures; original parent/context
  gates remain required after physical comparison. Full software contracts
  and interrupted recovery passed, with all 17 false pair/group exports
  rejected. Three jobs are queued under pre-estimated CPU-only resources;
  production and calibrated biological interpretation remain pending.

- October 2 UTC: preserve every original context and tied reference when
  integrating correspondence controls. Add all 27 mask/method/core scenarios
  using intersections of exact 48-order-pair bitmaps; use joint both-mask,
  both-method and both-core requirements without favorable combinations.
  Physical availability cannot override original parent/distinct-model/
  three-pair/native coorthology gates. Missing physical results stay NULL,
  distinct from measured zero-pass bitmaps. Keep lexical choices and tie order
  fixed; exclude empty sets from all-tie policies. Require independent original
  model SHA/source-gate reconstruction, complement-union bitmap arithmetic,
  SQL policy/denominator checks and all source/checkpoint exports, followed
  by original two-journal/full-hash closure. Regenerate and check deterministic
  1,000-context checkpoint chunks under an exclusive lock on interrupted
  restart; never duplicate a live job or restart completed output. Full
  software contracts and recovery passed, rejecting all 18 false exports.
  Complete production is queued after original comparison closure. Repeated
  eligibility decisions are dependent screens, not events, tests, accepted
  duplication or calibrated effects. Biological/phylogenetic/prediction
  uncertainty remains; GPU prediction stays paused.

## October 2 UTC: preserve complete signed contrast sensitivity

Retain all nine mask/core combinations, including joint scenarios, and every
selected structural alignment order for each of the 27,056 source-ready triples.
Report original-role AR-minus-BR extrema and missing/nonunique states before any
contextual or biological interpretation. Export exact-zero and 1e-9 Å numerical
direction classes; the numerical boundary is not an effect-size cutoff, and
sensitivity ranges are not confidence intervals. Strict all-selected-order
quality qualification preserves existing inherited exclusions. Never select a
favorable mask/core/order or promote an excluded result because its sign is
consistent. The independent reader must reconstruct the full raw-fit grid and
every exported field/count/checkpoint. Repeated scenario decisions remain
dependent. Full source context/tie/native-parent gates, accepted phylogeny,
prediction/domain/PAE controls and calibrated inference are still required.
[Complete scope](full-triad-contrast-sensitivity-20261002.md).

## October 2 UTC: separate reference-pool direction from quality eligibility

Preserve all original contexts, lexical references and tied choices when linking
closed physical contrasts. A physically measured shared triple cannot override
parent, model-identity, three-pair or native-guide exclusions. Null gated
directions do not mean zero contrasts. Summarize every eligible native-both tie
in the any-eligible pool, reporting category disagreement and direction support;
uniform agreement of an incomplete eligible subset does not establish complete
original-tie agreement. All-tie pools require all original references eligible
and nonempty; missing lexical references are never replaced. A quality-passing
context may have disagreeing reference directions. Retain quality flags and
direction states separately; do not choose a favorable reference/mask/core/order
or claim biological polarity/significance from these conditional descriptions.
[Complete full-context workflow](full-triad-context-contrasts-20261002.md).

## October 2 UTC: require joint direction across correspondence sources

Retain separate sequence, structural and joint signed contrast envelopes across
every selected mask, sequence method, core and input order. Deduplicate actual
geometry keys rather than repeating structural fits for each aligner. The full
joint envelope has 56 actual fits; source-category agreement involving sign
uncertainty is not stable sign or effect concordance. Missing/nonunique grids
cannot become complete. Require intersection of all 48 original order-pair
screen bits across every selected leaf; complementary failures and favorable
mask/method diagonals must not promote qualification. Keep numerical direction
boundaries distinct from biological effect sizes and uncertainty intervals.
Independent original raw-fit SQL extrema/grid/status reconstruction and separate
Cartesian order-flag/complement-union checks precede full hashes and both original
journals. Sequence-derived mappings still use predicted coordinates and do not
remove predictor circularity. Original parent/native/model/reference/tie gates
remain mandatory for subsequent contextual interpretation. Full production is
queued; all eight scientific aims remain incomplete.
[Complete joint-direction workflow](full-triad-joint-directions-20261002.md).


## October 2 UTC: carry joint directions through complete original reference pools

Apply all 27 joint mask/method/core direction scenarios to every original
context and tied reference. Preserve explicit scenario axes, source identities,
lexical order and parent/distinct-model/three-pair/native gates; physical overlap
cannot override source exclusions. Keep gated null direction distinct from
near-zero measurements and geometry qualification separate from category
agreement. Any-eligible pool direction uses all eligible ties and retains all
original missing/ineligible references; complete original-tie agreement requires
a nonempty fully eligible pool. Reject favorable method/reference selection,
incomplete axes and altered support/count/type fields using independent original
model-role SHA/gate and SQL pool/count reconstruction. Require full hashes and
both original completion/resource journals before production acceptance.
Full production is queued, not completed. These dependent context descriptions
do not establish biological polarity, significance or a duplication effect.
[Full contextual joint-direction workflow](full-triad-context-joint-directions-20261002.md).


## October 2 UTC: close complete whole-protein numerical evidence before integration

Require full Cartesian original fit/audit manifests, every input/fit checksum,
every five-fit audit shard and all flagged follow-ups before integrating complete
optimization results into a new comparison version. Preserve passes, flags and
original/follow-up errors with exact source/readback identities. Neither a
checksum nor a completed process promotes an unverified numerical error. Bind
all original fit/audit/follow-up/reader completion/resource journals and complete
source hashes; preserve the old five-candidate export as its own frozen version.
Numerical closure is not global optimality, identifiable covariance components,
accepted phylogeny or calibrated biological inference. Also retain older launch
schemas in the live inventory rather than excluding those original jobs or
rewriting their provenance.
[Full handoff and live inventory scope](full-whole-protein-optimization-closure-20261002.md).

## October 2 UTC: integrate all optimization follow-ups with explicit selected parameters

Require complete four-journal optimization closure before integrating every
original flag into a new full comparison version. Preserve the old five-case
version. Select closed refinement/recovery results while retaining unresolved
status; on follow-up error retain the original audited flagged fit and record
the failed source explicitly. Export selected coefficients and conditional
covariance separately from millions of repeated scalar comparisons; verify
source selection and unit transformations independently. Conditional covariance
does not establish calibrated uncertainty. Bind identical-plan checkpoint
recovery to deterministic parameter hashes and require full readback and both
original journals before completion. No favorable review/error promotion,
negative-gain clamping, transfer to expanded measurement cohorts or automatic
biological/model acceptance.
[Full workflow, validation and resources](full-whole-protein-comparisons-20261002.md).

## October 2 UTC: preserve correspondence sensitivity and target-selection shifts in inference

Publish every closed physical mask/method/core direction and retain full
eligible/reference uncertainty when moving into original contexts. Sign
uncertainty over alternative measurements is not a confidence interval or
ancestral polarity. Do not infer counts of sign reversals by subtracting
different quality cohorts. Publish every fixed-matching scenario's post-screen
balance and all three target-selection baselines. Aggregate scenarios only as
descriptive extrema/medians with explicit estimability, never independent
events or an optimally selected matching specification. Retained target–control
balance does not establish representative duplication-event coverage; observed
selection shifts remain requirements for missingness/sampling sensitivity and
calibrated expanded models. Preserve failed publication contracts and prior
figures while correcting immutable versioned output paths/labels.
[Complete background results](full-screened-background-results-20261002.md),
[closed joint directions](full-triad-joint-directions-20261002.md#complete-production-and-published-results-october-2-utc).


## October 2: retain gene contexts through the expanded measurement handoff

Index the complete expanded fixed matching by ordered target/background node
IDs. Keep physical model-pair identity as a separate reuse descriptor; genes
sharing predictions are not interchangeable independent observations. Preserve
all selections, endpoint mappings, scores, ties, unmatched scenarios and
mask/screens. All four target/control alignment-order comparisons remain
explicit; means require both orders and complete envelopes require all four
cells. Native TM dissimilarity remains labelled as native, with its length and
prediction-source limitations. Raw excluded diagnostics are quarantined by
numerical usability, and target coverage's intentional blanking is preserved.
Do not reuse older fitted inputs as evidence of fitting this expanded cohort.
Independent full SQL/source/Decimal readers and original journal/full-hash
closures gate all handoffs. New integration failures retain their original
versions and proof records; corrected source-lineage and quarantine contracts
use fresh versions. Expanded modelling, phylogenetic/reconciliation uncertainty,
ascertainment/dependence/prediction controls and calibration remain required.
[Full workflow and schemas](full-expanded-matched-measurements-20261002.md).


## October 2: recompute dependence inputs for the full expanded cohort

Do not reuse the earlier 4,568-pattern factors as coverage of the expanded
9,812-pattern design. Preserve all gene and versioned model endpoint occurrences
and target/control reuse, across guides. Export signed and unsigned incidence
separately; they represent different prospective models, with cancellation
explicit. Derive shared-entity family components and species patterns from
every closed logical case, including excluded/same-model/zero-pattern records.
Independent graph and raw-tree-edge checks must cover the entire export.
The five current kernel trees remain working substitution-distance alternatives,
not accepted dated phylogenies or a complete dependence correction. Expanded
designs, fitting, calibration and the other biological controls remain required.
[Complete resources, sources and interpretation](full-expanded-covariance-20261002.md).


## October 2: preserve ancestral failures and separate integrity from mixing

Close full 1,620-chain accounting and both original journals before using the
new terminal ancestral handoff. Retain all three memory failures and all five
allocation-warning chains, including integrity-passing outputs. No quartet
passes every scalar screen at this horizon, so individual scalar passes cannot
qualify ancestral posteriors. Recover failed chains with isolated same-seed
attempts under higher memory allowance, retaining the original model, horizon,
inputs and timeouts and running serially. Do not join failed samples or change
scientific thresholds. Full quartet overlays, warning review and length/
alignment/category/scalar mixing remain required before interpretation.
[Full accounting and recovery evidence](baliphy-method-assessment-20260927.md#full-initial-horizon-completed-and-source-closed-october-2).


## October 2: full ancestral diagnostic publication and matching-version review

Retain all 405 original quartets and both burn-in cuts, including unresolved
groups with blank diagnostics. Separate integrity, scalar/length/category
mixing and joint posterior qualification. The corrected coordinate-schema
closer must compare all original contributing chain coordinate files and
complete original artifact/journal closure. Publish complete count partitions
and inspect actual plots. Audit original scalar logs before extending sampling
or inferring that memory alone caused native allocation failures. Nonfinite
gamma values confined to discarded initialization and explicitly handled by
matching software are not automatically labelled software defects. Preserve
the repeated higher-memory failure and source/version evidence; do not trim
taxa, truncate priors, repair logged values or change diagnostic thresholds.
[Complete publication and remaining requirements](baliphy-full-first-horizon-diagnostics-20261002.md).


## October 2: preserve the complete expanded model-input grid

Retain every logical case under both masks and all five order contrasts.
Keep distinct gene endpoints despite shared predictions; link full expanded
covariance indices rather than old-subset factors. Transform identity and
distance powers before averaging orders. Missing structural measurements stay
null and raw observed zero stays zero. Preserve original matching decisions,
all 54 scenarios, all fixed screens and both eligibility gates. Check every
numeric input with separate Decimal calculations and every cohort/reuse count
with original-membership SQLite reconstruction. The full 622,080 future
model-setting records and 3,110,400 nominal tree-setting fits require rank/
recipe/resource qualification before fitting. No GPU authorization changed.
[Complete workflow](full-expanded-model-inputs-20261002.md).


## October 2: exact cohort reuse and explicit design qualification

Retain all 622,080 expanded model settings, including empty and unestimable
ones. Share computation only for exact ordered original case identities,
mask/order/axis/degree/response and source contracts; preserve distinct genes
sharing models. Record every declared term and deactivate only exactly zero
columns, with no approximate cutoff. Keep nonzero collinearity, absent sequence
information and constant response explicit. Check rank stability across the
recorded tolerance band using independent long-double normalization and
pivoted QR/gesvd. Complete original-membership SQL reconstruction, every
source/artifact hash and both original process journals gate production
acceptance. No expanded variance fit is launched until the full census
informs a separate compute/resource plan.
[Full design inventory](full-expanded-model-designs-20261002.md).

For version2, verify condition-number exports through their singular-value
ratio and compare reciprocal condition using a dimension/epsilon error bound.
A direct rank-boundary contract demonstrated that raw condition numbers
amplify otherwise acceptable last-singular-value roundoff. Rank thresholds
and review dispositions remain unchanged; the original version1 queue was
stopped with exact original handles before production and journal-preserved.


## October 2: shared genes/models and explicit variance-basis identities

Preserve each original entity namespace and distinct genes despite shared
coordinates. Retain signed endpoint cancellation rather than inventing
noise for a cancelled contrast. A family effect on the target-control
contrast is a separate one-per-case definition; unsigned family incidence
and this contrast intercept have proportional covariance kernels and cannot
be independently estimated together. Check complete kernel dependencies
before final covariance model identities/fitting, and repeat qualification
within each actual cohort. Fixed full-case numeric benchmarks are software/
resource evidence, not estimated variances, biological effects or a pilot.
[Full bank and numerical backend](full-entity-operator-bank-20261002.md).

## October 2: qualify covariance bases within each setting and its fixed design

All 4,250,692 original selection records establish one control per target
within each of 432 nonempty matching strata. Under a uniform residual diagonal,
the setting-specific target-node covariance equals residual identity. Do not
fit these as separately identified variances. The global pooled-case kernel
differs because targets recur across settings; its full Gram rank cannot
replace cohort-specific qualification. Nonuniform weighting requires its own
audit. Signed endpoint-family variance remains zero, while unsigned endpoint
family and contrast-family intercept variances are proportional.

Qualify the residual-space bases `H K H` as well as raw covariance bases before
REML fitting. Explicitly retain unresolved cancellation and rank boundaries.
No generic independent-column selection is authorized to change the
nonnegative variance cone. The new numerical audit passed all 24 dense
error-contrast tests and 13 invalid-input rejections; complete production
cohort/design qualification still awaits source closure.

## October 2: retain the complete uniform covariance qualification grid

Queue every original expanded design under both loading modes and all five
working trees, with all 6,220,800 setting/mode/tree links retained. Represent
target-node plus uniform residual variance as their exact combined variance;
represent unsigned endpoint-family plus contrast-intercept variance as
`variance_intercept + 4*variance_endpoint_family`. Preserve signed-family
zero in the audit. These specific exact combinations preserve the covariance
cone; do not drop any other dependent basis automatically. Nonuniform weights
remain a separate required qualification, not implied by this uniform stage.

Reuse component products across fixed designs. The independent reader caches
sparse global entity overlaps and reassociates species contractions to avoid
dense latent-by-species caches for every tree. Independent producer/reader
qualification and exact SQLite source-setting links gate numerical acceptance;
full source/artifact/two-original-journal closure gates the handoff. Source
review, constant-response, empty, dependency, rank-boundary and numerical-error
states remain explicit. No fit or calibrated inference is established.
[Full stage contracts, resources and scope](full-uniform-covariance-qualification-20261002.md).

## October 2: independent spectral candidate checks and complete fit accounting

Verify supplied ML/REML candidates with an eigen/SVD whitening implementation
independent of production Cholesky algebra. Compare coefficients in balanced
units, preserve conditional uncertainty labels, and require explicit score
boundary checks. Local curvature must be positive on a conservative superset
of the critical cone; weak boundaries stay tested. Flat/negative/asymmetric or
resolution-dependent curvature remains review. The empirical stability screen
is not a global-optimum proof or calibrated inferential uncertainty.

Preserve all original settings and exact source/input identities in prospective
uniform ML/REML fitting. Full arithmetic permits at most 5,208,000 unique
candidates and 12,441,600 setting links; exact source sharing reduces repeated
computation without shrinking the study. Retain all nonqualified source rows.
Close original design/qualification inputs, specify restartable full exports
and independent global checks, and measure full-scope timing before production
launch. No biological pilot, new GPU schedule or paid resource is introduced.
[Validation, limits and resource inventory](shared-entity-likelihood-20261002.md).

## October 2: checkpoint full fits and independently qualify every finite candidate

Bind full-grid candidate identity to closed qualification, original response/
design/audit identities and all backend pins/settings. Finalize each cohort
atomically and verify it on resume. Keep source exclusions, exact zero columns,
constant responses, dependencies and optimizer failures in the setting links.
Refuse completed reruns, including alternate reader output names. Rebuild only
derived reader scratch under its exclusive lock and unchanged input state.

Use separate spectral SLSQP starts, objective/KKT checks and conservative local
curvature before calling a working candidate independently numerically audited.
Replay every finite production start and reproduce declared failures, so error
records cannot conceal source cases. Agreement supports numerical qualification
but proves neither global optimality nor calibrated biological effects. The
full 71-pin draft remains unlaunched while qualification/timing are incomplete;
all eight aims and nonuniform/control/calibration requirements stay open.


## October 2: preserve whole recovery attempts and finish full diagnostic accounting

All three higher-memory BAli-Phy attempts are terminal. Accept output integrity
for the one fully parsed successful attempt; retain the other two allocation
failures alongside their original failures. Verify the full original grid and
all source/artifact hashes, then select whole attempts without concatenation.
Full recovery closure checked 22,486 bindings and both original journals.
Recompute only the newly complete quartet and new chain state trace, reuse
unchanged closed results, and retain the full 405-quartet scope. Scalar/length/
category diagnostics remain convergence screens, not posterior qualification.
Preserve the failed first audit and correct legacy launch metadata through a
new version with command/plan agreement checks, not edits to captured launches.

## October 2: keep independent spectral contractions within verified components

Use component-local entity trace/energy arithmetic with explicit global
whitening fallback under cancellation. Retain the unchanged streamed checker,
original production optimizer, all numerical thresholds and full study scope.
Dense/80-digit and complete fitting-grid contracts passed. The 80-pin draft
supersedes an unlaunched 71-pin draft; neither is production fitting. Software
speedup is not a project ETA. Close full qualification and measure actual
qualified-input runtime before installing fitting resources and launching.

## October 2 evening: publish the recovered diagnostic overlay separately

The corrected full diagnostic closure checked 50,019 bindings and both original
producer/readback journals. Publish all 405 quartets at both unchanged cutoffs,
including 403 complete and two failed quartets. Preserve the earlier 402-quartet
publication, every failed attempt and all original input/alias/seed identities.
Derive allocation notices from each selected whole attempt rather than carrying
old failed-attempt warnings into replacement chains. Inspect native PNG and
rendered PDF and check every plotted total/percentage against serialized counts.
No complete quartet passes every scalar or length screen; increased integrity
coverage is not posterior convergence or accepted ancestral structures.

## October 2 evening: measure every eligible full-scope fitting group

Queue timing behind the original full qualification closure, without launching
variance fitting. Retain all potential 5,208,000 candidate identities/4,340 cohorts
and every original source review. Choose representatives within every eligible
cohort/loading/tree/method/response group by greatest active design dimension,
then condition, then candidate ID. Original covariance guards and both numerical
backends are measured at scaled variances zero and one; reviews stay unmeasured.
Preserve all search budgets and numerical thresholds in conditional extrapolation.
Tests cover the complete synthetic grid, genuine six-variance probes and separate
full eligible selection. Timing agreement does not qualify fits, prove runtime
bounds or supply a project ETA. Original pinned sources and queues remain fixed.
