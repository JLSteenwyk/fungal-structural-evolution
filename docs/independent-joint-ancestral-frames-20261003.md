# Independent joint ancestral frame reader

The new reader keeps ancestral amino-acid states, rate categories and residue
coordinates together within the same native JSON record. This supports the
uncertainty propagation needed for aim 8. It does not qualify that aim,
posterior mixing or full-grid native sampling.

`scripts/independent_joint_ancestral_frames.py` checks every runtime node,
including ancestors outside the four focal candidates. State indices use the
native alphabet `ARNDCQEGHILKMFPSTWYV`; category indices range from 0 to 3.
Every decoded state must match its nongap letter in the record's own
`alignmentLines`. Concrete observed residues must remain unchanged; source
`X` may resolve to a valid native residue. Duplicate or extra labels, invalid
lengths/support, malformed rate matrices and unexpected conditions fail.
Rates must be finite, nonnegative, constant across amino-acid states within
each category and have mean one across the four categories.

Candidate rows are sorted by source-node identifier. A production caller
must establish the source-to-runtime mapping by rooted descendant clades,
verify source/artifact provenance and require the actual runtime tip set to
match the observed tips. The standalone decoder receives that verified
mapping; it does not independently establish likelihood correctness.

For each saved frame, anchors are ordered by tip identifier and then by
ungapped observed residue offset. Two tips can anchor the same native
alignment column, so repeated columns are expected. Candidate states and
categories are projected to those anchors together. Both arrays use 255 for
a gap; it is outside both supports. Candidate residues in columns occupied
by no observed tip are retained separately with candidate index, native
column, state and category. They are not lost during projection.

The exports retain each frame's coordinates and category rates. Native
column numbers can change between draws; they are not stable homologous-site
identifiers across frames. Category indices also do not establish comparable
biological rate classes across models or samples. Separate legacy FASTA logs
remain separate draws and are never substituted for the joint sequence.

All ten NPZ arrays have exact serialized value, shape, dtype and key checks
with pickle disabled. The [data dictionary](../metadata/independent_joint_ancestral_frames_data_dictionary_20261003.tsv)
defines their axes. Variable-length unanchored arrays are per-frame outputs;
the reader does not concatenate samples or original chain attempts.

## Qualification and remaining work

The [software receipt](../metadata/independent_joint_ancestral_frames_software_validation_20261003_v1.json)
checks all 1,620 future role metadata records, 135 alignment/tree inputs and
324 aliases. It checks source-tip/candidate membership and counts observed
residue anchors. Complete rooted candidate mapping is checked against the
native fixture trees; full-grid native mappings remain a production gate.

Read-only replay of nine retained native frames across all three priors
checks 36 ancestral records and 233 ancestral residue/category pairs.
State projections agree with the existing separate amino-acid projection
oracle after explicit alphabet conversion. A separate per-cell category
oracle uses nongap rank within each candidate sequence, checking both
anchored and unanchored coordinates. Every exported array passes readback.
Sixty malformed-frame cases, fifteen altered exports and two strict JSON
cases are rejected. Eighty synthetic state/category support combinations
exercise candidate-only insertions, deletions, repeated anchors and an
unknown observed residue. These are software checks, not a fungal pilot or
full-grid native output validation. No inference was rerun.

Execution took 3.35 seconds elapsed and 3.57 seconds of self CPU, with
87,576,576 bytes self peak RSS. The command configured 8 GiB address space,
150 CPU seconds, 180 seconds wall timeout, 256 MiB per file and one-thread
BLAS. Self RSS is distinct from charged group memory; no group limit or
long-chain memory guarantee is inferred.

Across all future roles, three saved frames would require 587,929,824 bytes
for the two candidate state/category arrays alone. This excludes coordinate
arrays, unanchored residues, raw all-node JSON, legacy logs and filesystem
overhead. It is a storage component estimate, not a complete resource budget
or finish ETA.

Full native startup/readback/closure, the historical sampler/resource and
reader closures, full joint-frame native sampling qualification and longer
adequate posterior ensembles remain pending. Production stage integration
must preserve every failed role and unresolved quartet and check original
process, command, source and artifact provenance. Root/model/likelihood,
predictor, calibration, accepted phylogeny/reconciliation/dating and all
eight biological aims remain open. GPU prediction stays paused.

Reproduce in new output/receipt paths, retaining the stated caps:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout 180 /usr/bin/prlimit --as=8589934592 --cpu=150 --fsize=268435456 -- \
  /home/bizon/anaconda3/bin/python scripts/check_independent_joint_ancestral_frames.py \
  --output NEW_AUDIT_DIRECTORY --receipt NEW_SOFTWARE_RECEIPT.json
```
