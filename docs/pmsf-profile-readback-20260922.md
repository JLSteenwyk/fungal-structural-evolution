# PMSF readback and guide sensitivity on the profile alignment

Both completed full-taxon PMSF runs on the 49,027-site profile alignment passed
saved-profile, tree/report and empirical bootstrap readback on September 22.
Each contains all 526 taxa, 49,027 ordered site-frequency vectors, and 1,000
bootstrap trees with exact taxon identities and finite nonnegative branches.
All empirical split counts agree with NEXUS weights and saved tree support at
the precision of the printed support values. The NEXUS files round percentages
to integers, so exact replicate frequencies are retained separately.

| Check | Profile guide | MAFFT guide |
|---|---:|---:|
| Distinct observed bootstrap splits, including terminal splits | 1,247 | 1,256 |
| Internal ML branches | 523 | 523 |
| SH-aLRT ≥80 and empirical UFBoot ≥95 | 456 | 461 |
| ML–consensus Robinson–Foulds distance | 0 | 0 |
| Reported consensus likelihood improvement | 0.697786 | 0.292607 |

Within each run, the saved ML and reoptimized consensus trees have identical
topology. Their branch lengths differ. Both variants are retained; support is
attached to explicit splits rather than transferred by node number.
SH-aLRT labels were checked for presence, range and agreement between report
and tree, but the SH-aLRT procedure was not independently rerun. Likelihoods
were read back, not recomputed. Both reports warn that 28 sequences contain
more than 50% gaps or ambiguity.

The two guide-conditioned ML trees share 516 of 523 internal splits (RF=14),
confirmed independently using the native OrthoFinder tree parser. No pair of
incompatible splits meets both stated support thresholds in both runs. This
does not establish identical trees or biological correctness. The mean
per-site L1 distance between the two estimated amino-acid profiles is 0.017408;
the maximum is 1.003811. Profile differences mean that cross-run likelihoods
must not be treated as a direct preference test between the resulting trees.

The [IQ-TREE PMSF documentation](https://iqtree.github.io/doc/Complex-Models#posterior-mean-site-frequency-model)
describes guide-based profile estimation. This comparison addresses guide
sensitivity for one alignment. The remaining crossed alignment/guide runs,
rooting, gene-tree discordance, taxon and marker sensitivities, and model
adequacy remain necessary before selecting a final species framework.

Reproduction:

```bash
python scripts/check_species_pmsf_audit.py
python scripts/audit_species_pmsf.py --run results/phylogeny/pmsf-profile-profile-recovery-v2 --output NEW_PROFILE_AUDIT
python scripts/audit_species_pmsf.py --run results/phylogeny/pmsf-profile-mafft-v1 --output NEW_MAFFT_GUIDE_AUDIT
python scripts/compare_pmsf_guide_sensitivity.py --run-a results/phylogeny/pmsf-profile-profile-recovery-v2 --run-b results/phylogeny/pmsf-profile-mafft-v1 --audit-a NEW_PROFILE_AUDIT --audit-b NEW_MAFFT_GUIDE_AUDIT --output NEW_COMPARISON
```

Completed outputs are in `results/phylogeny/pmsf-profile-profile-readback-v1`,
`pmsf-profile-mafft-readback-v1`, and `pmsf-profile-guide-sensitivity-v1`.
Receipts are archived as `metadata/pmsf_profile_profile_readback.json`,
`metadata/pmsf_profile_mafft_readback.json`, and
`metadata/pmsf_profile_guide_sensitivity.json`. These are small CPU-only
readbacks using existing files; no inference was restarted. The comparison
tables preserve split identities and missing branches explicitly.
