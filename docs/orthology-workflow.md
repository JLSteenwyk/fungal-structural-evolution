# Full-cohort orthology workflow

The corrected representative input contains 526 taxa and 5,815,847 proteins under `results/gene_representatives/full-v2`. All source proteins and alternative-product decisions remain available. The receipt is copied to `metadata/gene_representatives_receipt.json`. Provisional ORFs, unresolved gene features and the multi-gene protein are explicitly flagged; their inclusion in family discovery does not authorize interpreting them as resolved gene-copy counts.

OrthoFinder 3.1.5 is installed from its official Linux Intel release archive, verified against the release API's SHA256 digest. This release contains a Python package and bundled executables. Installation uses Python 3.12 in `.cache/envs/orthofinder`; conda and pip dependency locks are under `environments/`. Reproduce with `python scripts/install_orthofinder.py`. The installed command is `.cache/envs/orthofinder/bin/orthofinder`; help and version checks are part of installation verification. Runtime dependency checks remain necessary before launching the actual analysis.

The intended computation follows the [official scalable workflow](https://orthofinder.github.io/OrthoFinder/tutorials/advanced-tutorial/): infer reference orthogroups from a phylogenetically diverse core, then assign the remaining taxa and run the combined phylogenetic analysis. All 526 taxa are part of the experimental design from the outset. Core inference is a computational dependency, not a pilot or a reduced final dataset.

Before launch, freeze the core and additional-taxon manifests, retain every distinct species, verify FASTA checksums and identifier handling, and document resource estimates. Core selection must span the represented fungal lineages and suitable outgroups; high broad-BUSCO scores alone must not exclude evolutionarily reduced lineages. Default FAMSA alignment, DIAMOND searches, FastTree gene trees and hierarchical orthogroups will be evaluated against the task's accuracy needs and confirmed installed configuration. Full gene-tree reconciliation and independent species-tree comparison remain required. A simple sequence cluster or top-level MCL group is not sufficient evidence of orthology.

## Frozen inputs and launched core inference

`python scripts/prepare_orthology_inputs.py` created immutable symlink inputs under `data/orthology_inputs/v1` and recorded checksums for all 526 proteomes. The 64-taxon core contains 54 fungi and 10 outgroups, with 597,213 proteins. The additional cohort contains 462 taxa and 5,218,634 proteins. `metadata/orthology_input_manifest.tsv` records each species' stage and selection reason; the receipt records coverage and source hashes.

Core selection first requires one representative from each of the 19 fungal groups and eight outgroup categories, then adds a second choanoflagellate and ichthyosporean. Remaining fungal slots maximize taxonomic-prefix diversity, with broad-BUSCO single-copy recovery used as a tie-break. This produces 15 ascomycetes and 15 basidiomycetes plus 24 fungi from other groups. Mandatory group coverage prevents high broad-BUSCO scores from displacing reduced lineages. This is a taxonomic proxy; selection sensitivity and the eventual inferred phylogeny remain relevant.

`python scripts/run_orthology_core.py` launched OrthoFinder with 32 search threads, eight analysis workers, DIAMOND, FAMSA and FastTree. Runtime dependency checks passed. Configurations and logs are under `results/orthology/core-control-v1`; output is under `results/orthology/core-v1/Results_Sep12`. Existing outputs are protected from automatic overwriting. A failed/interrupted stage requires explicit restart review using OrthoFinder's supported continuation mechanism.

Core computation completed successfully in 6,066.6 seconds with all 64 expected taxa. The output reports 597,213 genes, 510,318 in orthogroups (85.4%), 86,895 unassigned genes and 24,215 orthogroups. Full assignment, combined gene-tree analysis, reconciliation, quality assessment and independent species-tree comparison remain outstanding. The core output alone does not complete this project's orthology requirement.

## Full-cohort assignment

`run_orthology_assignment.py` advances the validated 64-taxon MSA core using `--assign` for all 462 remaining taxa (5,218,634 proteins), following the official scalable workflow. It rechecks every input FASTA checksum, exact core-tree taxa, pinned executable and core completion receipt. The combined run uses 32 search threads and eight analysis workers, DIAMOND, FAMSA and explicit FastTree (not the assign-mode fastest default). The software creates a new named result directory beside the core output; protected core tree, orthogroup and overall-statistic checksums are recorded and checked after completion. Internal reusable databases may be added by OrthoFinder. All 526 final tree taxa must be present before the controller records computational completion; scientific reconciliation validation and sensitivity are separate.

Pre-launch planning envelope: the completed core used 6.6 GB disk and 101 minutes. The full assignment includes 8.74 times as many additional proteins as the core but uses a different algorithm, so linear time extrapolation is inappropriate. Reserve 1 TB disk headroom, up to 512 GB aggregate host RAM and 1–14 days as a provisional planning range, not an observed estimate or enforced cap. The host has approximately 1 TB RAM and more than 10 TB free disk. These settings reuse released core CPU resources alongside ongoing tree, domain and single-GPU prediction jobs. No paid infrastructure is provisioned. Record actual performance from execution logs and monitor the large-family stages.

The first assignment attempt (full526v1) stopped during reference-profile preparation because NumPy 2.5.3 lacks the top-level `np.chararray` called by OrthoFinder 3.1.5. OrthoFinder printed the traceback yet returned exit code zero; the controller correctly withheld a completion status because no full species tree existed. The revised controller also exits nonzero for incomplete results. The isolated dependency was pinned to NumPy 2.2.6; pip dependency checking passed, and the previously failing function reconstructed a real core alignment exactly. `environments/orthofinder-core-pip.lock.txt` preserves the original core dependencies; the updated `orthofinder-pip.lock.txt` reproduces assignment dependencies. The install report records the wheel URL/hash. Failed outputs and runner source remain under assignment-control-v1 and Results_full526v1. A fresh full526v2 run uses assignment-control-v2 without deleting or overwriting the first attempt.

## Recovered family alignment failure during full assignment

The full526v2 log reported FAMSA exit 134 for `OG0001522`, followed by missing
alignment/tree errors. The main assignment remained live in its MSA/tree stage.
The family input contained 769 unique nonempty protein records (79–4,623 AA).
Isolated reruns with bundled FAMSA 2.2.3 and locally installed 2.5.2 both exited
successfully and preserved every identity and ungapped sequence. The original
failure did not reproduce; its cause is unresolved.

Recovery used the **same bundled 2.2.3 binary** as production, avoiding a method
change. Its 20,197-column alignment was trimmed with OrthoFinder's existing
`trim.main(...,10,0.1,500,0.75,False)` to 2,112 columns. Every retained column was
verified as an ordered column of the original alignment. Bundled FastTree with
the production invocation completed a 769-tip tree; identities and finite,
nonnegative branches passed readback. This gene tree is a production-equivalent
recovery, not an independent orthology validation or supported species tree.

Before installation, the assignment was verified live and still before the
species-tree stage; no family writer or Astral-Pro process was active. Both
production artifact paths were absent. The installer copied the validated
alignment and tree to independent inodes and atomically created only those two
missing files, refusing to overwrite any existing file. Final hashes match the
isolated recovery. No reference-core artifact changed. Downstream full-run
completion and scientific validation remain pending.

```bash
python scripts/recover_orthology_family_tree.py \
  --diagnostic results/orthology/famsa-OG0001522-diagnostic-v1 \
  --output results/orthology/OG0001522-recovery-v1
python scripts/install_orthology_family_recovery.py \
  --recovery results/orthology/OG0001522-recovery-v1 \
  --production-pid 841175 \
  --receipt metadata/orthology_family_recovery_installation_receipt.json
```

The PID records this execution and must be re-observed for another execution;
the installer refuses a missing/wrong process, advanced stage, active writer,
existing output or changed source. Diagnostic commands, executable hashes,
resource plans, recovery outputs and installation are recorded in
`metadata/orthology_family_*`. Raw intermediate alignments remain outside Git.


## Audited large-family repair installation handoff

The isolated `OG0000017` FastTree rebuild remains live. Its wrapper validates the tree with Bio.Phylo and leaves production untouched. The additional service `fungal-family-tree-installation-20260916.service` waits for that exact wrapper PID, creation time and command, then requires its successful completion receipt and unchanged source hashes. It independently parses the repaired tree with OrthoFinder's bundled ETE-derived reader, verifies the complete source tip universe and explicit finite nonnegative branch lengths, and records the parser checksums. This uses the actual native parser that downstream OrthoFinder consumes, independently of the repair wrapper's Bio.Phylo parser.

After validation, the handoff stages and fsyncs an exact tree copy, checks that the production target remains the audited empty file, and atomically installs the copy. The isolated repair, log and validation receipt remain available. An occupied target, changed input, absent completion receipt or failed native readback stops for review. No source alignment or running rebuild is modified. The installation controller and its validator are pinned while this service is active.

Configuration and launch identities are in `metadata/orthology_OG0000017_installation_controller_{config,launch}.json`; live state and subsequent readback are under `results/orthology/OG0000017-installation-controller-v1`. The stage uses one CPU, an 8 GiB memory planning allowance with a 16 GiB service ceiling, 1 GiB output allowance, and 0.01–1 hour after the rebuild for checks and installation. Launch requires 16 GiB available memory and 5 GiB free disk at the validation stage. No additional charges are incurred. The transient service does not survive reboot; a partial installation or failed state requires inspection before any fresh controller launch.

`python scripts/check_family_installation_gates.py` exercises the actual native parser and installer in disposable three-tip fixtures. It accepts and installs a valid tree, and rejects wrong tips, negative branches, omitted branch lengths, an occupied target and a missing completion receipt while preserving the target. The check record is `metadata/orthology_family_installation_gate_checks.json`. These fixtures verify the handoff logic, not large-tree runtime or biological topology.

Installation is a prerequisite to a separately reviewed native `--from-trees` / `-fgt` continuation. Reconciliation has not been launched by this handoff and remains required, together with species-tree/root sensitivities and independent output validation.


## Native continuation descriptor and remaining launch requirements

The interrupted full-cohort `Log.txt` contains its working-directory base but no `FN_Orthogroups` or `WorkingDirectory_Trees` records. Direct readback with the installed `PreviousFilesLocator_new` returns no cluster file and no species-tree paths. The script `prepare_orthology_continuation_descriptor.py`, run with the installed OrthoFinder Python, creates an explicitly constructed log in `results/orthology/full526-continuation-descriptor-v1`. It adds the locations of the independently audited cluster file and existing tree working directory, verifies those locations with the native reader, preserves the exact working-directory chain, and confirms the original log checksum is unchanged. No gene tree, species tree or inference output is created by this preparation.

The native reader's inferred rooted/unrooted species-tree files are both absent. Thus the descriptor is **not launch-ready**. The actual `--from-trees` path consumes the IDs species-tree file through `RootSpeciesTree`; a supplied `-s` option alone must not be assumed to create that file without checking native behavior. A reviewed species tree, explicit species-ID conversion, and isolated tree-input staging remain necessary. Homogeneous guides and root diagnostics are available, but supported species-tree and root sensitivities remain part of the scientific requirement.

The installed code also defaults `fix_files=True`; `--no-fix-files` disables it. The default continuation updates MSA/tree outputs and performs a second ortholog-inference call. Although `StartFromTrees` creates a new named results directory, it retains the previous tree working-directory pointer. Future execution must review these write paths and the consequences of either setting, rather than assuming that a new results directory protects source data. The exact command and resource plan are not yet finalized. These findings are archived in `metadata/orthology_continuation_native_source_review.json`; the descriptor receipt and constructed log are versioned alongside the reproducible preparation script.


## Executed native restart behavior check

`check_native_orthology_restart.py` executes two isolated four-taxon, three-family, 13-protein fixtures with the installed CLI, using `--from-trees`, a user species tree, `-M msa`, two threads, one analysis worker and `--no-fix-files`. One family has a known terminal duplication. These are disposable software fixtures, not a biological pilot or a replacement for full-cohort analysis. Inputs and outputs are under `results/orthology/native-restart-fixture-v1`; the resource plan caps each CLI attempt at 120 seconds and terminates the entire fixture process group on timeout.

The fixture with only the user-facing `-s` tree fails inside `RootSpeciesTree` because the expected IDs tree file is absent, yet the installed CLI returns exit code **zero**. It produces no final root-level HOG table or overall statistics. Therefore process exit code alone is not a completion criterion. The fixture that also prepares `WorkingDirectory/SpeciesTree_unrooted_ids.txt` with the same rooted topology completes. Its observed original source-file inventory and checksums are unchanged with `--no-fix-files`.

`audit_native_orthology_restart.py` independently checks all 13 genes across three HOGs and resolved trees, all 42 directed cross-species ortholog pairs, the one expected terminal duplication, the species-tree root split, overall counts, and exact source preservation in both cases. The observed zero-exit failure is explicitly retained. Output hashes, installed source hashes, the execution receipt and independent readback are archived as `metadata/orthology_native_restart_fixture_*`.

This establishes a usable restart path for the fixture. The full-cohort launch must prepare and verify matching species-name and native-ID trees, validate the repaired complete family-tree universe, isolate the tree-input directory, and freeze resources and output protections. It must check final HOG, resolved-tree, duplication, ortholog and species-tree contents; a zero process status or completion banner is insufficient. These tests do not establish biological accuracy or guarantee file behavior for other flags or datasets. The default file-fixing path remains distinct from the tested `--no-fix-files` path.


## Rooted guide inputs for full-cohort reconciliation sensitivities

`results/orthology/reconciliation-guide-inputs-v1` now contains both audited homogeneous guide alternatives (`profile` and `mafft`), each with a taxon-name tree and matching native OrthoFinder-ID tree. Both are rooted on the independently established 501-fungal-entry/25-outgroup separating edge. Equal division of that edge is an explicit file-format convention; it estimates neither the root coordinate nor elapsed time. These inputs support conditional topology sensitivities and do not replace the pending supported-mixture, marker and root analyses.

Preparation verified all 526 species IDs and preserved every original unrooted split and branch length. The installed native species-tree converter independently checked renaming/topology. Its default serialization rounds branch lengths (maximum observed edge discrepancies about 4.77e-6 and 4.58e-6), so the intended `species_tree_taxa.nwk` and `species_tree_ids.nwk` files preserve higher precision; `native_converter_check.nwk` is a diagnostic artifact only. All six output trees passed an independent undirected-graph readback: exact 526-tip identity, complete bifurcation, the required root partition, and 6,294 canonical-edge comparisons. The intended inputs preserve original edge sums within 1e-12, including the two root halves.

Reproduce preparation with `python scripts/prepare_reconciliation_guide_trees.py --plan metadata/reconciliation_guide_input_plan.json`, using a fresh plan/output directory for subsequent runs, then run `audit_reconciliation_guide_trees.py` with that plan and a fresh readback output. Input and script hashes, resources, preparation and independent-readback receipts are versioned under `metadata/reconciliation_guide_input_*`. The original guides remain unchanged. Isolated native tree-input staging, completed repair installation, full reconciliation, supported-guide replacement/sensitivity and biological output validation remain necessary.


## Measured full-family workload

`estimate_full_reconciliation_workload.py` read all 104,138 audited families (4,485,438 clustered genes; 526 taxa) and checked each native family index/count against the artifact inventory. Per-family taxon multiplicities imply 5,177,610,933 unordered cross-species candidate gene pairs, or 10,355,221,866 directed pairs. Two independent linear-time formulas agreed for every family; no pairs were materialized. These are possible within-family comparisons, not inferred orthologs. There are 355,559,436 occupied unordered family/taxon pairs. The top 25 families account for 47.94% of candidate protein pairs, but no family has been excluded on that basis.

A flat directed-pair table at 64/128/256 bytes per pair would occupy approximately 617/1,234/2,469 GiB. These are arithmetic scenarios, not predictions of OrthoFinder's grouped output, runtime or peak memory. The current tree input files occupy 136,011,525 bytes across 32,638 expected filenames, with OG0000017 still empty pending repair. The installed default keeps only O(number of species) ortholog files open concurrently; the potential directed species-pair file count is 276,150. The final resource plan should retain grouped ortholog output and avoid unnecessary Cartesian expansion, while estimating native output and scratch separately.

The workload inventory also flags a material interpretation issue: OG0000017 has 46,787 genes from Austropuccinia psidii (F181123; GCA_000469055.2), 83.58% of that family and 34.92% of its 133,973 representative proteins. Native species IDs, manifest mapping and source receipts agree. This requires annotation/repeat/homology review before treating these memberships as biological duplications; it is not evidence of a lineage-specific expansion by itself. The existing repair continues unchanged. Counts and provenance are recorded in `metadata/orthology_OG0000017_taxon_concentration.json`.

The full per-family workload table remains outside Git at `results/orthology/full-reconciliation-workload-v1/family_workload.tsv`; its checksum, compact top-family summary, resource estimate and export readback are versioned under `metadata/full_reconciliation_workload_*`.


## Full sequence-quality review of OG0000017

`profile_large_family_quality.py` profiled all 55,981 source sequences and their complete 27,697-column alignment in `results/orthology/OG0000017-quality-v1`. It retained every gene, verified source/alignment identity and ordered residue preservation, measured canonical amino-acid composition with two independent entropy formulas, counted exact sequence redundancy globally and within taxa, and accumulated non-gap/canonical column coverage. No exclusions or live tree changes were made.

The 46,787 Austropuccinia psidii members comprise 46,782 distinct sequences; only ten genes participate in within-taxon exact-duplicate groups. Thus repeated identical records do not explain most of the concentration. Median source length is 172 residues in this taxon versus 417 in other taxa, with median non-gap occupancy of 0.621% versus 1.506% of the full alignment. Median source-retention fraction is 1.0 in both groups. Median global canonical composition entropy is 4.135 versus 4.142 bits; this does not rule out local low-complexity segments, repeats or gene-model artifacts.

Of 27,697 columns, 6,903 contain any focal-taxon residues, 24,791 contain any other-taxon residues, and 3,997 contain residues from both groups. These are descriptive coverage counts, not a homology-confidence threshold or evidence that every represented gene pair shares informative columns. The heterogeneous lengths and sparse alignment reinforce the need for domain/homology and annotation review before interpreting reconstructed duplications. They do not identify a repeat class, establish biological copy-number expansion, or justify automatic gene removal.

An independent export readback checked all 55,981 row identities, exact-hash multiplicities, retained-length ratios, 27,697 column rows, the 15,381,640 retained-residue total across both axes, and all 105 group quantiles using NumPy. It did not independently reparse every raw alignment character. Resource plan, source pins, compact group summaries and receipts are under `metadata/orthology_OG0000017_quality_*`; full protein/column tables remain outside Git.


## Isolated retained-family reconciliation inputs

`scripts/stage_full_reconciliation_inputs.py` prepared independent physical copies for the profile and MAFFT guide alternatives at `results/orthology/full-reconciliation-inputs-v1/{profile,mafft}`. Each contains 33,168 copied files, approximately 3.007 GB: all 526 species FASTAs, the full native species/sequence mappings and cluster partition, 32,637 previously validated gene trees, and matching named/native-ID guide trees. The FASTAs were resolved across the full526 and original64 working directories. No symlinks or hardlinks to the original inputs are used.

`scripts/readback_staged_reconciliation_inputs.py` independently checked the exact complete file universe, every copied/source checksum and size, absence of source inode aliases, and the native sequence-counting routine. Both variants resolve all 526 species FASTAs and 5,815,847 representative proteins. This input-protein count is distinct from the 4,485,438 genes in the audited family partition; staging does not assign previously unclustered proteins to families.

The native file locator resolves only the isolated working directory in each case. Both directories deliberately contain `Log.pending.txt`, not the `Log.txt` required for a restart. OG0000017 remains absent until its rebuilt tree has passed independent validation and installation. The staged inputs are therefore incomplete and no reconciliation has launched. Final launch still requires the repaired tree in both variants, a complete final inventory, a resource plan, and output validation beyond the CLI exit status. The intended native path is `--from-trees --no-fix-files`; alignments and similarity-search intermediates are not copied. The previously executed native fixtures cover that path, but full-run source-access coverage remains to be verified during execution.

Reproduce staging with `.cache/envs/orthofinder/bin/python scripts/stage_full_reconciliation_inputs.py --plan metadata/full_reconciliation_input_staging_plan.json`, using a fresh output path if rerunning. Reproduce the independent readback with `.cache/envs/orthofinder/bin/python scripts/readback_staged_reconciliation_inputs.py --plan metadata/full_reconciliation_input_staging_plan.json --output <new-readback.json>`. Source pins, resource allowances, compact receipts and readback are versioned under `metadata/full_reconciliation_input_staging_*`; per-file manifests and all copied data remain outside Git.


## Full protein accounting and unfinished discovery stage

The full native input contains 5,815,847 proteins. `audit_full_protein_family_coverage.py` checks every SequenceIDs entry and all 526 species FASTA headers, and partitions the exact native identities into 60,423 retained singleton proteins, 22,154 proteins in two-member families, 4,402,861 proteins in tree-eligible families, and **1,330,409 proteins outside the retained partition (22.88%)**. The independent `readback_full_protein_family_coverage.py` uses OrthoFinder's own cluster parser and checks every exported outside-protein identity and label against the native complement. All per-taxon totals close.

The 64 core taxa contribute 597,213 proteins and have no proteins outside the saved partition. The 462 added taxa contribute 5,218,634 proteins, including all 1,330,409 outside-partition proteins. This is a workflow-stage asymmetry, not evidence of lineage-specific novelty, loss or missing biological homologs. Full-proteome domain annotation and atlas coverage must continue to include these proteins. Their complete identifiers are exported outside Git in `results/orthology/full-protein-family-coverage-v1/outside_retained_partition.tsv`; compact per-taxon coverage, checksums and receipts are versioned.

Review of the installed `BetweenCoreOrthogroupsWorkflow` establishes that the interrupted accelerated workflow saved its assigned partition and unassigned FASTAs **before** species-tree-dependent clade discovery. The original working directory contains only the initial cluster file and lacks its native species-tree products. The `--from-trees` entry point bypasses this discovery workflow. Thus the staged restart inputs and their successful software fixtures establish a route to **retained-family reconciliation only**. They do not establish a complete continuation of orthology discovery. This qualification supersedes any earlier implication that the repaired tree alone would make full orthology ready for final reconciliation.

The intermediate discovery FASTAs contain 1,390,832 proteins: all 1,330,409 outside-partition proteins plus exactly all 60,423 retained singleton proteins. The exact overlap was checked with the native cluster reader. `write_unassigned_fasta` excludes multi-member families, so the singleton overlap is expected at this intermediate stage. A completed combined partition must replace the old singleton membership for genes assigned to new families, retain appropriate unresolved singletons, and pass unique-membership and complete accounting checks; blindly appending groups would not satisfy those requirements.

Executing installed `get_new_species_clades` on each conditional rooted guide gives 45 disjoint clades covering all 526 taxa. The largest clade contains 355,735 discovery-input proteins; the maximum taxon count is 72 for the profile guide and 71 for MAFFT. The sum of squared clade protein counts is 199,178,446,126 versus 192,712,000,282. These are possible ordered protein-pair dimensions, including self comparisons, not executed DIAMOND alignment counts or runtime predictions. The clade lists, source hashes and exact arithmetic are recorded in `metadata/unfinished_clade_discovery_scope_review.json` and reproduced by `review_unfinished_clade_discovery.py`.

Required next work is to complete discovery with a resource estimate and isolated outputs, resolve singleton replacement, infer added family trees, evaluate guide-dependent and cross-clade missed-homology sensitivity, and reconcile the expanded validated partition. Existing retained-family trees and staged copies remain useful and unchanged. No new discovery or reconciliation inference was launched by this audit.


## Complete discovery inputs for both guide alternatives

`prepare_clade_discovery_inputs.py` stages all 1,390,832 intermediate discovery proteins without filtering, preserving their original native-ID FASTA headers. The two guides define 90 clade occurrences, of which 30 clades have identical complete species membership across guides. These reduce to 60 distinct input sets at `results/orthology/clade-discovery-inputs-v1/clades/<clade_id>/proteomes`, using deterministic species-set identifiers. Guide-dependent clades remain separate. No partial-clade, approximate or sequence-similarity-based reuse is performed.

The stage contains 676 independent FASTA copies (479,409,438 bytes), comprising 1,782,124 protein records and 452,966,656 residues across the 60 unique input sets. Copies preserve the complete audited source bytes and do not link to live inputs. `clades.json` maps each native clade index under each guide to its staged input set, and `files.tsv` records species, source checksums, counts, lengths and paths.

`audit_clade_discovery_inputs.py` independently reads every copied FASTA and verifies hashes against sources, physical isolation, protein IDs, species identity, counts and residue totals. Both guide partitions have exactly the same 1,390,832 distinct protein IDs and no repeated IDs within a partition. Eight unique clades contain only one taxon; these are retained explicitly and require a verified single-taxon discovery path, rather than omission or automatic singleton assignment. Discovery inference has not yet launched. Remaining prelaunch work includes the native workflow contract, resource plan and exact mapping/merge validation; the earlier singleton-replacement and guide-sensitivity requirements still apply.

Run preparation with `.cache/envs/orthofinder/bin/python scripts/prepare_clade_discovery_inputs.py --plan metadata/clade_discovery_input_plan.json` (fresh output directory required). Run readback with `.cache/envs/orthofinder/bin/python scripts/audit_clade_discovery_inputs.py --plan metadata/clade_discovery_input_plan.json --output <new-readback.json>`. Compact receipts and the pinned preparation plan are versioned under `metadata/clade_discovery_input_*`; sequence data and full file/clade manifests remain outside Git.


## Executed native discovery software contract

`check_native_clade_discovery.py` exercised the installed CLI on deterministic random amino-acid strings in 24 exact-copy groups, with one or two synthetic taxa. This is a software fixture, not a biological pilot. `audit_native_clade_discovery.py` independently checked original identifier mappings and every expected family membership from the native cluster export.

For two taxa, all 96 identifiers are recovered in the expected 24 four-member groups. With default file-fixing behavior, `--only-groups` nevertheless builds 24 gene trees and runs downstream inference. Adding `--no-fix-files` produces the same family partition without gene trees. The intended discovery-only command therefore includes `--scores-v2 --only-groups --no-fix-files -I 1.2`, with explicit thread limits. These tests establish control flow and exact fixture membership, not biological accuracy or full-scale normalization equivalence.

The `-f` entry point rejects a single taxon in both file-fixing modes, but returns exit code zero. The eight single-taxon clades must therefore not use that entry point or be treated as completed based on exit status. `check_single_taxon_search_restart.py` demonstrated a native alternative using `-b` with a complete single-species self-search and its exact SpeciesIDs/SequenceIDs mapping. `audit_single_taxon_search_restart.py` independently recovered all 48 input proteins in the expected 24 duplicate-pair families, without trees or changes to the input artifacts. This path uses one real input taxon and adds no artificial species.

For `-b`, the installed CLI rejects `-o` and again returns zero despite failure. Successful restart outputs are nested beneath the isolated working directory at `WorkingDirectory/OrthoFinder/Results_<name>`; original mapping files remain in the input working directory. The executed negative case and successful output layout are recorded explicitly. Full production work still requires verified self-search generation for those eight clades, a resource plan and artifact-level completion checks. Native success messages or exit status alone are insufficient.

Final fixture outputs remain outside Git at `native-clade-discovery-fixture-v2` and `single-taxon-search-restart-fixture-v3` under `results/orthology`. Compact receipts, resource estimates and readbacks are versioned under `metadata/native_clade_discovery_fixture_*` and `metadata/single_taxon_search_restart_*`. The earlier synthetic observations remain preserved outside Git. No biological discovery batch was launched by these software tests.


## Full clade discovery batch launched

The complete 60-clade queue is running under `fungal-clade-discovery-20260916.service`, with the largest 355,735-protein clade first. All 52 multi-taxon and eight single-taxon clades are included, covering both guide alternatives through the explicit shared-clade map. This is full production execution, not a biological pilot. Shared inputs account for 1,782,124 protein records across distinct runs, representing the 1,390,832-protein discovery universe under each guide.

`run_clade_discovery.py` verifies pinned source files and complete original identifiers before each clade. Multi-taxon cases use native `-f`; single-taxon cases generate a native-ID mapping and an actual DIAMOND self-search, validate every search row, and use native `-b`. Self-search settings reproduce the installed DIAMOND configuration (BLOSUM62, gap penalties 11/1, more-sensitive, E-value 0.001 and compressed tabular output), with 16 search threads; other DIAMOND defaults are retained. All discovery commands specify scores-v2, MCL inflation 1.2, only-groups and no-fix-files. A completed clade must have one native cluster export, a bijection to its original protein IDs, a unique complete family partition, unchanged source FASTAs and no gene-tree inference. Family memberships and artifact hashes are written before its completion receipt.

The exact production runner was executed against synthetic one- and two-taxon fixtures, including new self-search generation. `readback_clade_discovery.py` independently used BioPython and the native MCL parser to verify all 48/96 identifiers and 24 expected exact-copy groups. The same independent readback is queued for all production outputs after discovery, via `advance_clade_discovery.py`. These checks establish identity and software behavior, not biological accuracy or full-scale normalization equivalence.

The batch runs one clade at a time with 16 search threads and four analysis threads, a 192 GiB service memory cap, no swap allocation, a 1 TiB output allowance, and explicit available-memory/disk checks. It requires 256 GiB available memory before a new clade and stops its own active process group for review if available memory falls below 32 GiB, disk free space falls below 1 TiB, or the output allowance is exceeded. Each command has a 336-hour upper stop limit. Completed clade receipts can be reused under the same plan; partial clades require review and are not silently overwritten. Other project jobs are unchanged.

The initial total planning envelope is 4–168 hours, not a measured ETA or confidence interval. The previous 597,213-protein core run took 6,066.6 seconds with 32 search threads but included downstream work and different sequence composition. The new clade pair dimensions, composition, hit sparsity and MCL memory prevent direct runtime transfer. Production timings will refine the estimate. The resource plan uses the existing host and creates no paid resources.

Execution plans, software/input pins, launch process identity and fixture evidence are under `metadata/clade_discovery_execution_*` and `metadata/clade_discovery_runner_fixture_*`. Live state and outputs are in `results/orthology/clade-discovery-v1`; its controller schedules independent readback automatically. Scripts, the execution plan and the pinned installed OrthoFinder sources must remain unchanged while this queue runs. Guide-level merging with singleton replacement, added gene trees, reconciliation and cross-clade homology sensitivity remain required after completion.


## Guide-specific merge prepared; awaiting discovery completion

`merge_clade_discovery_partitions.py` prepares separate complete protein partitions for the two guides, gated on the entire discovery queue and its independent readback. It keeps all 43,715 original multi-protein families intact and replaces the 60,423 original singleton slots through each guide's complete discovery partition. The resulting union must contain all 5,815,847 original protein IDs exactly once. Duplicate assignments, overlap with retained multi-protein families, unknown IDs, repeated family identifiers and incomplete coverage fail before a guide's merged files are written. No real merged partition has been generated yet; the controller described below now queues merging after successful discovery and readback.

Dropping old singleton slots can leave gaps in native family indices. The merge therefore assigns new contiguous indices and writes a mandatory `family_sources.tsv` crosswalk back to every retained or discovered source family. `singleton_relocations.tsv` records the destination and size for every replaced singleton. Exact membership hashes accompany source mappings. Original retained tree names are recorded only as candidates for later validated reuse; no trees are copied or assumed complete by this stage. Downstream scripts must use the crosswalk, never assume old and new OG numbers identify the same family.

`check_discovery_merge.py` passed two alternative synthetic guide partitions and rejected ten invalid cases, including missing proteins, duplicate membership and retained/discovery overlap. `check_discovery_merge_files.py` then executed the complete file-based merge command using previously independently audited synthetic discovery results. Both synthetic guide labels share the same discovery output by design; all 99 proteins, 25 serialized families, contiguous family indices, the retained-tree source crosswalk and replacement of an old singleton into a four-protein family were checked. These software tests do not substitute for independent readback of the eventual biological merged partitions.

The pinned plan is `metadata/guide_discovery_merge_plan.json`: one CPU, an 8 GiB memory allowance, 2 GiB output allowance and a broad 5–60-minute planning range on the existing host. Run `python scripts/merge_clade_discovery_partitions.py --plan metadata/guide_discovery_merge_plan.json` only after discovery and independent readback finish. Its receipt explicitly requires an independent merged-partition readback before added-family tree inference and expanded reconciliation. Original discovery and live retained-family inference inputs remain unchanged.


## Automatic merge and independent validation handoff

`readback_guide_discovery_merge.py` independently parses original and merged cluster files with the installed native MCL reader, reconstructs every retained/discovered source membership, and verifies the full protein universe under each guide. It also checks raw serialized multiplicity (the native parser uses sets), contiguous family numbers, source memberships, tree-candidate references, all membership hashes and every singleton relocation. A complete source-to-output checksum chain is required. Tree candidates remain provenance references, not validated installed trees.

The checker passed the complete synthetic file fixture. Five additional fixture mutations deliberately changed a source family, tree candidate, singleton target, duplicated gene or omitted gene. Artifact and plan checksums were updated in each mutated fixture so acceptance could not hinge on a digest mismatch; all five failed semantic validation without producing a success readback. Results are under `metadata/guide_discovery_merge_readback_fixture.json` and `metadata/guide_merge_readback_rejections.json`.

`fungal-guide-discovery-merge-20260916.service` is live and waiting for the exact discovery controller process identity. After that process terminates, it requires the discovery controller's successful terminal state and independently checked output hashes before running the merge and its complete readback. The service has a 16 GiB memory cap, an 8 GiB planning allowance, no swap allocation, and requires at least 32 GiB available memory and 20 GiB free disk before execution. The combined merge/readback planning envelope is 0.2–3 hours after discovery completes. No added gene trees or reconciliation are launched by this handoff.

Controller configuration and launch identity are versioned as `metadata/guide_discovery_merge_controller_{config,launch}.json`; live state is `results/orthology/guide-discovery-merge-controller-v1/state.json`. Successful biological outputs will remain at `results/orthology/guide-discovery-merge-v1`, including `readback.json`. The merge plan and controller/auditor scripts are now pinned for this live handoff and must not be edited in place. This supersedes the earlier unqueued preparation status; biological merging and validation are still pending.


## Expanded family tree workload (2026-09-17)

The completed guide-specific partitions were enumerated before scheduling
expanded alignments and trees. Exact memberships, rather than renumbered OG
labels, identify reusable tasks. Across both guides there are 96,587 distinct
families with at least three proteins: 32,638 retained tree candidates and
63,949 new tasks. The profile guide requires 62,846 new trees and the MAFFT
guide 62,956; 61,853 new tasks have identical memberships in both guides.
The largest new family contains 10,072 proteins. These counts do not provide
a runtime estimate: sequence lengths, alignment dimensions and method choice
must inform the execution resource plan.

All singleton and pair families remain in the source partitions; they are
accounted separately from inferred trees. The profile partition has 500,113
singletons and 62,586 pairs (11,077 retained plus 51,509 newly discovered).
The MAFFT partition has 500,481 singletons and 62,447 pairs. No biological
novelty or lineage restriction is inferred from these counts.

`scripts/census_expanded_family_tree_workload.py` verifies merged artifact
hashes and each serialized membership against its source crosswalk, then
writes `results/orthology/expanded-family-tree-workload-v1/unique_tree_tasks.tsv`.
`scripts/readback_expanded_family_tree_workload.py` uses the installed native
MCL reader to independently check every task's size, taxon count, source,
representative index and guide multiplicity: 96,587 tasks covering 191,078
guide-specific tree occurrences. Receipts are archived under
`metadata/expanded_family_tree_workload_*`.

This completes workload enumeration only. Expanded family FASTAs, new
alignments and gene trees still need execution; retained tree candidates need
validation against their exact memberships before reuse. The original large
family repair remains running. Reconciliation is not yet complete.


## New-family sequence staging (2026-09-17)

`scripts/stage_expanded_family_sequences.py` stages all 63,949 distinct new
families under `results/orthology/expanded-family-sequences-v1/families/`.
Filenames use exact native-ID membership hashes so later jobs need not depend
on guide-specific renumbering. The input plan pins the independently verified
task census, merged partitions and the 526-proteome source manifest.

The staged collection contains 788,325 distinct native protein IDs, 824,364
family-sequence records and 240,132,391 amino-acid residues in 249,384,973 FASTA
bytes. Some proteins occur in differing families under the alternative guide
partitions; these occurrences are preserved. Full sequences and native IDs
are retained without trimming or sequence deduplication. The per-family table
records protein/taxon counts, residue totals, minimum/maximum lengths and file
checksums. Staging was planned for one CPU, 16 GiB memory and a 5 GiB output
allowance, using existing resources.

Independent validation uses `scripts/readback_expanded_family_sequences.py`
and Bio.SeqIO against the original source FASTAs, rather than the copied
staging sources. It checks every full sequence, exact family membership,
file checksum and recorded dimension. A successful `readback.json` is required
before inference. Alignments and new trees have not yet been run; large-family
memory requirements and native alignment settings must be included in the
execution plan. Retained-tree reuse and reconciliation remain separate stages.


## Expanded alignment and tree execution (2026-09-17)

All 63,949 distinct new-family tasks are running under
`fungal-expanded-family-trees-20260917.service`, using
`scripts/run_expanded_family_trees.py` and the pinned
`metadata/expanded_family_tree_execution_plan.json`. Inspection of the original
full526v2 log established FAMSA, rather than MAFFT, as its alignment method.
The new run matches that method: FAMSA 2.2.3 with one thread, installed native
OrthoFinder trimming (minimum sequence threshold10, meaningful-character
fraction0.1, minimum columns500, retained-character fraction0.75), then the
same default FastTree executable. Untrimmed alignments are retained.

Sixteen independent workers process all tasks, largest residue totals first.
Each worker has a 16 GiB address-space ceiling; the service has a 320 GiB
memory cap with swap disabled. New tasks require 64 GiB available memory and
1 TiB free disk; each command has a 72-hour timeout. The output planning
allowance is 256 GiB (not a filesystem quota), and the initial total planning
range is 4–168 hours, not a calibrated ETA. No paid resources are used.

The runner checks exact source IDs and untrimmed residues, rectangular
alignments, trimmed columns as an ordered subset of raw columns, nonempty
trimmed sequences, and exact tree tips with finite nonnegative branch lengths.
Per-family receipts and artifact hashes allow verified completed tasks to be
reused. Partial tasks require review and are not automatically overwritten.
Failures remain explicit while other families continue. After all jobs pass,
`scripts/readback_expanded_family_trees.py` independently reads every tree
using OrthoFinder's native parser and checks explicit branch serialization,
source tips and all saved artifact hashes. Completion of the job service alone
does not establish this readback or biological accuracy.

`check_expanded_family_tree_runner.py` exercised the production worker using
three- and twelve-sequence synthetic families. Identical sequences with
distinct IDs were preserved. The twelve-sequence alignment was trimmed from
900 to 600 columns, retaining 7,200 of 7,500 residues. Checkpoint reuse passed,
and wrong alignment IDs and wrong tree tips were rejected. These are software
fixtures, not a biological pilot or runtime benchmark. Fixture and live launch
receipts are versioned in metadata. This transient service requires explicit
state review after a machine restart. Retained-tree validation and full
reconciliation remain outstanding.

## Retained-tree membership catalog (2026-09-17)

`scripts/catalog_retained_family_trees.py` maps every retained tree candidate
to the profile and MAFFT expanded partitions through their audited membership
hashes and source-family crosswalks. It reads immutable staged trees from
`results/orthology/full-reconciliation-inputs-v1/profile/`, verifies their
manifest hashes, and uses the native OrthoFinder parser to check unique tips,
the exact membership fingerprint, taxon count and explicit finite nonnegative
branch lengths. An absent staged tree remains an explicit pending row; an
unmanifested or changed tree is an error. The running OG0000017 repair is not
read while it is being written.

The output is `results/orthology/retained-tree-catalog-v1/`. Its TSV records
source paths, hashes and both destination family IDs; the receipt records
validated and pending totals. This is a catalog for later assembly of complete
reconciliation inputs, not a reconciliation result. Run with:

```bash
.cache/envs/orthofinder/bin/python scripts/catalog_retained_family_trees.py \
  --census results/orthology/expanded-family-tree-workload-v1 \
  --staged results/orthology/full-reconciliation-inputs-v1/profile \
  --merged results/orthology/guide-discovery-merge-v1 \
  --output results/orthology/retained-tree-catalog-v1
```

The validator uses one CPU, reads roughly 3 GB of existing staged inputs,
and writes a small catalog without copying trees. Re-execution requires a
new output directory to preserve the earlier evidence.
