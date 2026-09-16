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
