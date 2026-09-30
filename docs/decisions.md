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
