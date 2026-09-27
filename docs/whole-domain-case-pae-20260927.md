# PAE for the whole-protein/domain inspection cases

All 39 exact model/version combinations used by the 13 inspection cases now have
retrieved PAE matrices. The preparation script checked each complete local
coordinate sequence against the triad sequence hash and verified the PDB hash;
only matching versions were requested. One cache receipt was initially present.
Retrieval used two download workers, one CPU and a 2 GiB memory cap, with no GPU
prediction or paid resources.

The retrieval and independent readback both completed successfully. All 39
matrices were verified, with zero failed models; 11,389,053 entries were checked
for dimensions, type, finiteness, nonnegative values and the declared maximum.
Compressed and uncompressed hashes, versioned URLs, sequence provenance and
lengths were also checked. The matrices occupy approximately 3.8 MB compressed
and 27.7 MB as JSON.

This supplies confidence information for the next regional analysis; it does
not establish confidence in the observed domain arrangement or provide
independent validation. Domain-core versus outside-region PAE summaries still
need to use the actual retained residue sets and both directional matrix blocks.
The case coordinates and mappings remain unchanged.

Preparation: `scripts/prepare_case_pae_snapshot.py`.
Retrieval: `scripts/retrieve_marker_pae.py --snapshot
results/structural_comparisons/whole-domain-case-pae-inputs-20260927-v1 --output
<new-output-path>`.
Independent audit for this fixed run: `scripts/readback_case_pae.py`.
The audit waits for the recorded retrieval process and requires its terminal
success before checking the complete manifest and all matrix entries.

The input snapshot is
`results/structural_comparisons/whole-domain-case-pae-inputs-20260927-v1`;
the retrieval manifest is in
`results/structural_comparisons/whole-domain-case-pae-20260927-v1`.
Versioned evidence is in
`metadata/whole_domain_case_pae_completed_20260927.json` and
`metadata/whole_domain_case_pae_readback_20260927.json`. Per-model receipts and
large matrices remain in `data/structures/pae` outside Git history.

## Regional confidence results

`scripts/summarize_case_regional_pae.py` completed all 7,488 regional rows:
832 partitions × three protein roles × three matrix blocks. Blocks are internal
domain PAE with diagonal entries excluded, rows inside/columns outside, and
rows outside/columns inside. Raw row/column labels preserve both directions
without relying on plot-axis conventions. Mean, median, 90th percentile,
maximum and fractions at or below 5, 10 and 15 Å are retained for every block.
The latter thresholds are descriptive sensitivity settings, not calibrated
acceptance criteria. Independent scalar indexing, sums, rank interpolation and
threshold counts verified every unique region; complete serialized readback
passed. All 117 case/role/block summaries remain available.

For the three cases passing both anchor and outside n50/c70 coverage, ranges of
the **regional 90th percentile** across both inter-region blocks, roles and
mapping alternatives are:

| Case | Inter-region PAE p90 range (Å) |
| --- | ---: |
| Heliocybe OG0000054 | 5–8 |
| Jaapia OG0000294 | 13–20 |
| Phycomyces OG0002812 | 9–10 |

The Jaapia B model has only about 1.2–3.8% of these directional inter-region
entries at or below 5 Å, depending on direction and alternative. Its apparently
arrangement-sensitive coordinate contrast therefore requires particular caution.
This is an inference about confidence limitations, not proof that the contrast
is an artifact. AlphaFold's documentation explains that high PAE between domains
indicates uncertain relative positions/orientations:
[AlphaFold DB confidence guidance](https://alphafold.ebi.ac.uk/faq).
Here the outside region can include more than one domain or linker, so a
single whole-block statistic should not be treated as a specific domain-pair
orientation assessment.

PAE and observed RMSD contrasts are different quantities. These PAE summaries
are not contrast error bars, p-values or independent evidence; lower values in
Heliocybe and Phycomyces do not validate their small differences. Independent
prediction or experimental evidence remains necessary before assigning a
mechanism. No alternative prediction has been started while GPUs are paused.

Results are in `results/structural_comparisons/case-regional-pae-20260927-v1`;
versioned validation is in `metadata/case_regional_pae_completed_20260927.json`.
