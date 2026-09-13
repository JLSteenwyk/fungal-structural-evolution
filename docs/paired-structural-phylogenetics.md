# Matched sequence and structural branch estimates

## Executed checkpoint

The full 526-taxon, 125-marker input space was audited (65,750 taxon-marker rows). The immutable 422-model AFDB acquisition snapshot supplies sufficient confident observations for 52 markers: 4–12 taxa and 91–1,181 retained columns per marker, with 105 taxa represented across these alignments. The other 73 markers remain explicitly listed as coverage-limited. This is a checkpoint within the full sampling design; ongoing acquisitions and predictions will expand coverage.

All 208 requested fits completed. They provide 416 matched unrooted branches, each with one amino-acid estimate and three structural-model estimates. These are point estimates conditional on marker sequence topologies, confidence-selected sites and substitution models. They are not the supported species phylogeny, change per year, physical displacement, or evidence of structural acceleration.

## Observation and alignment rules

`prepare_paired_phylogenetic_inputs.py` verifies the audited encoding, structural mapping and original profile-matrix receipts, model coordinate checksums, encoding sequence hashes, and residue correspondence. Each retained observation requires a canonical aligned amino acid, exact physical residue correspondence, a valid native 3Di feature, all six feature-context residues at pLDDT ≥70, and maximum directional PAE within that context ≤10 Å. Invalid terminal features are missing even when their printed state is the ordinary coil symbol. Nonfinite confidence is missing.

Both alphabets receive exactly the same missing-data mask (`?`). Taxa require at least max(50, ceil(0.3 × original marker columns)) observations; markers require four eligible taxa. Only columns that are entirely missing among eligible taxa are removed. Invariant columns remain, and every output column maps to the original marker and profile alignment. The emitted FASTAs are read back to verify dimensions, taxa and masks. All excluded taxa and markers remain in the coverage audit. These eligibility thresholds do not prove sufficient phylogenetic signal.

## Models and fitting

The authors' [Edmond deposit](https://doi.org/10.17617/3.1MJJBH), version 3.0, supplies the [Garg and Hochberg model](https://doi.org/10.1093/molbev/msaf124) files Q.3Di.AF and Q.3Di.LLM. Retrieval verifies publisher MD5, size and local SHA256. Each file has 190 positive finite lower-triangular exchangeabilities and 20 positive frequencies whose sum is unity within six-decimal publication rounding. Published bytes are retained unchanged.

The [IQ-TREE custom-model format](https://iqtree.github.io/doc/Substitution-Models#user-defined-empirical-protein-models) uses state order ARNDCQEGHILKMFPSTWYV. Here those letters denote structural states, not amino acids; the files contain exchangeabilities plus stationary frequencies, not a preconstructed instantaneous rate matrix. `-st AA` selects the 20-letter parser while the explicit custom file defines structural evolution.

IQ-TREE 3.0.1 searches an unrooted LG+F+G4 amino-acid tree for each matched marker. On this fixed topology (`-te`), it fits structural branches under Q.3Di.AF+G4, Q.3Di.AF+F+G4 and Q.3Di.LLM+G4. The AF matrix with published frequencies is the reference; empirical structural frequencies and the LLM-trained matrix are sensitivities. LLM here describes the published matrix's training source; our observed structural states were extracted from coordinates, not predicted directly from amino acids by ProstT5. The current sequence model also requires adequacy/sensitivity assessment.

Identical sequences remain (`-keep-ident`). Four independent workers use one CPU thread and a requested 2 GB memory allocation per fit. The resource plan preceded execution: 208 fits, provisional 0.1–10 minutes each, 0.1–9 hours at four workers, and 2 GB output headroom. Actual summed elapsed time across individual fits was 160.4 seconds. This short runtime does not forecast full-cohort or uncertainty-analysis costs.

Unrooted branches are identified by canonical taxon bipartitions, with a hash of the complete marker taxon set. A degree-two display root is suppressed by summing its two edge segments. Different marker taxon sets do not produce interchangeable branch identities. Structural fits must preserve the sequence splits exactly. No arbitrary display root is interpreted as an evolutionary ancestor.

## Audit findings and limits

The post-run audit verifies all input and result checksums, the model recorded in each IQ-TREE report, paired masks, complete marker coverage, tree taxa/topologies, branch-table values and total tree lengths. All 208 fits passed those integrity checks. There were warnings in 169 fits, and 85 fits had at least one branch ≤10⁻⁵ substitutions/site. No fitted branch reached 10 substitutions/site. The last threshold is a descriptive flag, not a saturation test or confidence bound.

Warnings include rare states, near-zero internal branches, high missingness, some tree-search convergence diagnostics and small-tree memory messages. They are preserved verbatim in `metadata/paired_marker_fit_warnings.tsv`. An audit of the pinned IQ-TREE source found that its memory-slot minimum can exceed its leaf-count cap for very small trees, triggering the “Too low -mem” message even with a large requested allocation; the source URL, checksum and reasoning are in `metadata/iqtree_small_tree_memory_audit.json`. This allocation diagnostic does not resolve the substantive rare-state or near-zero-branch concerns.

Do not divide structural by sequence branch length when the denominator may be near zero. Do not compare likelihoods across different observed alphabets. Comparisons among models fitted to the same structural alignment require parameter penalties and model checks. Site resampling must be paired between alphabets, and uncertainty work must consider dependence among overlapping and spatially linked 3Di features, topology uncertainty and confidence-driven ascertainment. Broader structure coverage, supported genealogies, model adequacy, alignment/indel sensitivity, prediction-source controls, direct-geometry benchmarking of branch estimates and shared-ancestry statistical analysis remain required before ranking accelerated clades or testing sequence–structure coupling.

## Reproduction and evidence

Run from the repository root, using new output directories if artifacts already exist. The fitting runner is restartable with unchanged configuration and verified completed outputs; an incomplete IQ-TREE job retains its native checkpoint.

```bash
python scripts/prepare_3di_models.py --output data/structural_models/garg-hochberg-v3
python scripts/prepare_paired_phylogenetic_inputs.py \
  --encodings results/structural_alphabet/audited-marker-v2 \
  --snapshot results/structural_markers/snapshot-v2 \
  --matrix results/phylogeny/profile-matrix-50-v1 \
  --output results/phylogeny/paired-structural-inputs-v1
python scripts/run_paired_marker_fits.py \
  --inputs results/phylogeny/paired-structural-inputs-v1 \
  --models data/structural_models/garg-hochberg-v3 \
  --output results/phylogeny/paired-marker-fits-v1
python scripts/summarize_paired_marker_fits.py \
  --inputs results/phylogeny/paired-structural-inputs-v1 \
  --models data/structural_models/garg-hochberg-v3 \
  --fits results/phylogeny/paired-marker-fits-v1 \
  --output results/phylogeny/paired-fit-audit-v1
python -m unittest discover -s tests -p 'test_paired*.py'
```

Seven focused tests cover invalid coil masking, physical correspondence rejection, nonfinite confidence, missing-state site statistics, invalid substitution files, display-root invariance and changed/invalid trees. Full-data validation is provided by the preparation and completed-fit audit, not by those tests alone. Compact receipts, marker coverage, fit summaries, warnings and 416 matched branch rows are versioned under `metadata/paired_*`. FASTAs, coordinate/state arrays, reports, downloaded files and detailed coverage remain outside Git with checksums and paths in receipts.

## Paired resampling and feature dependence

The next executed stage requests 200 paired draws for each of 52 markers under column/block lengths 1, 10 and 30: 31,200 paired draws and 62,400 fixed-topology fits. A draw samples circular blocks uniformly by starting column, concatenates them and truncates to the original alignment length. Both alphabets and all taxa use exactly the same indices. Circular wraparound is a resampling convention, not a claim that protein termini are biologically adjacent. Branches and gamma shape are re-estimated, with AA empirical frequencies recomputed for each draw. Structural fits use the published AF frequencies. Unestimable draws (an entirely missing taxon or no variable columns in either alphabet) are recorded without replacement; unexpected executable failures stop the affected work for inspection.

The resource plan reserves eight CPU workers, one thread and a requested 2 GB per fit, with 20 GB output headroom. A provisional 0.1–5 seconds per fit gives approximately 0.2–11 hours before initialization/I/O. Each request, deterministic PCG64 seed, source columns, alignment, command, tree, diagnostic log and receipt is retained. Per-draw checksums and complete-batch receipts support audit and continuation. The full run's status is established by its terminal receipt, not this planned count.

A feature-overlap audit of the actual paired inputs found 83,455 observed features across 286 taxon-marker alignments and 474,044 within-model pairs of features sharing at least one physical residue. Of these pairs, 223,510 (47.1%) have circular alignment separation ≥10 columns and 138,415 (29.2%) separation ≥30. Such pairs cannot co-occur within one sampled block of the corresponding length. These are potential dependencies, not measured covariances, independent observations or effective sample sizes. The finding shows why block-length sensitivity alone cannot establish valid uncertainty for all nonlocal structural dependencies.

The summarizer verifies every draw's source columns, FASTAs, checksums, model identity, topology and recorded branch lengths before calculating 2.5/97.5 percentiles. It requires at least 90% of requested draws to be estimable for a branch interval. These are conditional sampling-sensitivity intervals, not overall calibrated confidence bounds: protein-site stationarity, topology uncertainty, alignment uncertainty, model fit, prediction errors and nonlocal feature sharing remain unresolved. Paired resampling covariance can later inform an errors-in-variables analysis, but it does not itself measure evolutionary sequence–structure coupling.

```bash
python scripts/run_paired_branch_resampling.py \
  --inputs results/phylogeny/paired-structural-inputs-v1 \
  --models data/structural_models/garg-hochberg-v3 \
  --fits results/phylogeny/paired-marker-fits-v1 \
  --audit results/phylogeny/paired-fit-audit-v1 \
  --output results/phylogeny/paired-resampling-v1
python scripts/audit_3di_feature_overlap.py \
  --inputs results/phylogeny/paired-structural-inputs-v1 \
  --encodings results/structural_alphabet/audited-marker-v2 \
  --snapshot results/structural_markers/snapshot-v2 \
  --output results/structural_alphabet/paired-feature-overlap-v1
python scripts/summarize_paired_resampling.py \
  --inputs results/phylogeny/paired-structural-inputs-v1 \
  --fits results/phylogeny/paired-marker-fits-v1 \
  --resampling results/phylogeny/paired-resampling-v1 \
  --output results/phylogeny/paired-resampling-audit-v1
```

The [original phylogenetic bootstrap paper](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x) motivates character resampling under independence assumptions. Our fixed-topology paired resampling addresses branch-estimation sensitivity rather than bootstrap clade support, and the native feature-overlap audit directly limits that independence assumption here.

### Completed resampling audit

All 156 marker/block batches finished. Of 31,200 paired draws, two were unestimable: marker 4977163at2759, block 30, draw 19 had an entirely missing taxon; marker 5001734at2759, block 30, draw 134 had no variable structural columns. Both remain recorded and were not replaced. The independent audit validated 62,396 fits and produced 2,496 interval rows (416 branches × three block lengths × two alphabets). Warnings were present in 47,590 fits and remain in the original reports. All interval batches had at least 199 estimable draws. Outputs occupy approximately 2.2 GB.

Median AA interval widths were 0.1447, 0.1552 and 0.1516 substitutions/site for blocks 1, 10 and 30; the corresponding structural widths were 0.06374, 0.06793 and 0.06431. Thirty-column blocks widened intervals for 203/416 AA branches and 202/416 structural branches compared with single columns. Thus larger blocks do not uniformly widen intervals, and similar medians do not validate independent sites or resolve nonlocal feature dependence. No acceleration ranking follows from these results.

![Conditional interval-width sensitivity](figures/paired_resampling.svg)

The figure is generated with `python scripts/plot_paired_resampling.py --audit results/phylogeny/paired-resampling-audit-v1 --output results/phylogeny/paired-resampling-figures-v1`. Five focused tests cover resampling length/bounds/contiguity, paired missingness, deterministic seeds, remote feature overlap and graph components. All passed. Metadata receipts, conditional intervals and batch summaries are versioned; complete replicate artifacts remain outside Git.

## Expanded paired input execution

`run_expanded_paired_inputs.py` is running with control output
`results/phylogeny/paired-expanded-control-v1`. It waits on the recorded live
confidence-controller process identity, requires its complete verified receipt,
and then invokes the existing paired-input preparation on
`audited-gdm-expanded-v1`, `gdm-expanded-v1` and the unchanged
`profile-matrix-50-v1` matrix. Output is
`results/phylogeny/paired-inputs-gdm-expanded-v1`; main execution output is
`logs/expanded_paired_inputs_v1.log` and the child log is `preparation.log` in
the control directory. Configuration and resource estimates are versioned in
`metadata/expanded_paired_input_config.json` and
`metadata/expanded_paired_input_resource_plan.json`.

The planned coverage comparison requires identical matrix, confidence-mask and
eligibility definitions across snapshots and records eligible taxa and retained
columns for every marker. Availability and chosen predictions can both change.
These are acquisition/selection diagnostics, not biological rate estimates.
The job uses a new immutable paired output; interrupted partial preparation
requires inspection and a new output configuration. Expanded input completion,
model fits, support, geometric calibration and uncertainty remain pending.

## Support-aware expanded fitting

The fitting runner now records marker/taxon/column dimensions from its verified
input summary instead of retaining the first snapshot's 52-marker resource
assumptions. An empty eligible set fails before fitting. Expanded runs require
a separate resource estimate based on their completed inputs.

Optional `--alrt 1000 --bootstrap 1000` adds sequence-topology SH-aLRT and
NNI-refined ultrafast-bootstrap analysis. All bootstrap trees are retained;
every replicate must contain the full original taxon identities and the total
must match the requested count. Structural models still estimate branches on
the same fixed AA tree. Default support counts remain zero for the earlier
point-estimate behavior. These options do not provide structural branch-length
intervals or propagate sequence-topology uncertainty into structural fits.
Support-aware expanded fitting has not yet been launched.

Three focused tests cover support only on the searched sequence tree,
rejection of invalid replication counts, and detection of missing, repeated
or changed bootstrap taxon identities and incomplete replicate counts. The
two existing unrooted-edge correspondence tests also pass. Existing outputs
retain their original producer hashes; use a new output for the updated runner.

## Lineage coverage of the earlier paired snapshot

A full-grid coverage summary independently checked the actual 52 paired AA/3Di
FASTA taxon sets against every per-taxon usable-marker count. All 105 retained
taxa are Ascomycota; only two have at least ten usable markers. There are 286
usable taxon/marker combinations. The source snapshot has models for 118 taxa,
also entirely Ascomycota. Its conditional branch analyses therefore provide no
evidence for structural-rate differences across fungal phyla or outgroups.

Across 526 × 49,027 matrix cells, the sequential mask assigns 135,903 observed,
3,009,758 noncanonical/missing sequence, 22,618,607 absent structural mapping,
84 invalid native feature, 23,736 low feature pLDDT and 114 high feature PAE.
These counts include all markers before per-taxon and per-marker eligibility;
they are not the dimensions of the final paired alignments. Missing matrix
sequence includes gaps and noncanonical symbols, and missing structures do not
establish biological absence. Sequential reason counts are not independent
causal effects of each filter.

Reproduce with `summarize_paired_lineage_coverage.py --inputs
results/phylogeny/paired-structural-inputs-v1 --output
results/phylogeny/paired-lineage-coverage-v1` using a new output for reruns. The
same script will assess the expanded paired snapshot after it completes. All
27 manifest lineage/role groups remain in the tables, including zero-coverage
groups. Compact tables and readback receipts are versioned under
`metadata/paired_*coverage*`.

## Expanded paired inputs completed; supported fitting started

The expanded acquisition snapshot passes the unchanged paired masks and
coverage rules for 124 markers, up from 52. Marker 4999134at2759 remains
coverage-limited. Emitted alignments have 5–135 taxa and 73–1,184 columns per
marker; their union includes 322 taxa (304 fungal entries, 18 outgroups).
Independent FASTA readback checked every AA/3Di identity set, dimension,
missingness mask and per-taxon usable-marker count.

Twenty-three of 27 manifest lineage/role groups have usable markers.
Aphelidiomycota, Calcarisporiellomycota, Sanchytriomycota and the Corallochytrea
outgroup group have none in this snapshot. All 526 entries remain in the
coverage audit and ongoing acquisition/prediction design. Group presence does
not imply dense coverage or sufficient support for every branch. The complete
matrix-level sequential mask retains 4,160,922 observations before per-marker
eligibility/column removal; this is distinct from whole-model state counts.

The 496-fit batch has started: each of 124 markers receives an AA LG+F+G4 tree
search with 1,000 SH-aLRT and 1,000 NNI-refined ultrafast-bootstrap replicates,
then AF+G4, AF+F+G4 and LLM+G4 fits on that exact AA topology. Four concurrent
one-thread jobs use up to 2 GB memory each. Input dimensions, a conservative
1–96-hour runtime/10-GB output planning envelope, software and commands are
pinned in `metadata/expanded_paired_fit_resource_plan.json` and
`metadata/expanded_paired_fit_config.json`. Actual optimization has begun.

```bash
python scripts/run_paired_marker_fits.py \
  --inputs results/phylogeny/paired-inputs-gdm-expanded-v1 \
  --models data/structural_models/garg-hochberg-v3 \
  --output results/phylogeny/paired-marker-fits-gdm-expanded-v1 \
  --alrt 1000 --bootstrap 1000
```

These are candidate-marker branch fits conditional on selected sites and
sequence topology. Gene-copy issues, model adequacy, topology and branch-length
uncertainty, phylogenetic placement, geometry calibration and predictor
circularity remain to be addressed before acceleration or selection claims.

### Coverage depth across lineages

![Paired lineage coverage and marker depth](figures/paired_lineage_coverage.svg)

The expansion increases represented taxa from 105 to 322, but depth remains
uneven. Of all 526 entries, 204 have no usable paired marker, 200 have 1–9,
17 have 10–49 and 105 have at least 50. All members of the last category are
Ascomycota. Thus the presence of 23 lineage/role groups does not establish
balanced coverage or readiness for comparisons of rates across fungal phyla.
The figure uses every lineage's complete manifest count as its denominator and
retains zero-coverage groups. Thresholds summarize coverage, not statistical
power or independent replication.

Reproduce with `plot_paired_lineage_coverage.py --previous
results/phylogeny/paired-lineage-coverage-v1 --expanded
results/phylogeny/paired-lineage-coverage-expanded-v1 --output
results/phylogeny/paired-lineage-figure-v1`, choosing a new output for reruns.
SVG, PNG and PDF artifacts are available there; the PNG was visually inspected.
Source table hashes, group/denominator equality and nonnegative depth partitions
are checked. The figure receipt and an SVG copy are versioned.
