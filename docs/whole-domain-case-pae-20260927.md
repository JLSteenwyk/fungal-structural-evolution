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
