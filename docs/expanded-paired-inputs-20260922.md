# Expanded paired sequence and structural-alphabet inputs

The next paired-input snapshot uses the complete 25,322-model ESMFold inventory
after combined mapping and cohort-specific confidence qualification pass.
It retains the existing source alignment, six-residue feature confidence mask
and eligibility rules: pLDDT at least 70 throughout the feature context,
maximum directional context PAE at most 10, a valid native state, canonical
amino acid and exact residue correspondence. Each taxon needs at least
max(50, ceil(0.3 × original marker columns)) observations; each marker needs
at least four eligible taxa. Invariant sites remain included.

`scripts/readback_paired_inputs_from_encodings.py` independently reconstructs
accepted residues from qualified arrays, rather than importing the paired
producer's character/mask function. It checks all taxon/marker eligibility
cells, reconstructed emitted AA and 3Di characters, identical observation
masks, retained columns, and source matrix correspondence. It does not
reassign structural states or establish the biological accuracy of the
predictions or alignments.

The full new readback passed on the previous combined ESMFold dataset and
matches the earlier emitted-file audit: 89 ready markers, 292 taxa,
9,332 marker/taxon cells, 1,699,035 observed paired cells and 22,205 retained
marker columns. All 526 × 125 eligibility cells were checked, including
ineligible combinations. The result is archived in
`metadata/esmfold_combined_paired_full_array_readback.json`.
The pre-launch validation allowance was one CPU, 16 GiB RAM, no swap, small
output and 0.1–4 hours of uncalibrated planning; the actual job used about
25 CPU seconds.

The controller `scripts/advance_complete_paired_inputs.py` is queued under
`fungal-all-completed-paired-inputs-20260922.service`. It waits for the exact
qualified-encoding integration process, requires its successful bound receipt,
runs the unchanged paired-input producer, and then runs the full independent
array readback on the expanded output. Its plan is
`metadata/esmfold_all_completed_paired_inputs_plan.json`.

The new inputs will be under
`results/phylogeny/paired-inputs-esmfold-all-completed-20260922-v1/`; controller
records and the baseline-to-expanded marker-coverage comparison will be under
`results/phylogeny/all-completed-paired-input-controller-v1/`.
Resources are one CPU equivalent, 32 GiB RAM, no swap, 10 GiB output allowance,
a 100 GiB free-disk gate, and 0.5–12 hours of uncalibrated planning after the
predecessor completes. No GPU prediction or paid resource is used.

The expanded ready-marker count is not known yet. Increased coverage is not
evidence of structural acceleration or a prediction-source effect. Supported
sequence trees, structural fits on matched genealogies, branch uncertainty,
direct-geometry checks and gene-copy/annotation sensitivities remain required
before updated evolutionary conclusions.
