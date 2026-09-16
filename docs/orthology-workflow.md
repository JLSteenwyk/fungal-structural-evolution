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
