# Matched domain measurements before effect inference

Production and full independent verification completed for all 150,080 shared boundary-domain matches across
both masks and all four target/background alignment-order combinations:
1,200,640 rows. No preferred order is selected or averaged. All
configuration links remain traceable to the audited selected-control inventory;
configurations without shared domains remain in upstream coverage accounting.

Each side retains aligned length, recomputed RMSD, sequence identity, separate
endpoint TM scores, original and retained-input coverage, and joint confidence.
Native left/right quantities are normalized to canonical interval A/B according
to input order; target and background A/B labels do not imply cross-comparison
residue correspondence. Six same-mask and six both-mask qualification flags
are attached using the full verified coverage/geometry tables. Unusable side
metrics are blank and numerical/geometry statuses remain explicit.

Differences in RMSD, sequence identity and aligned length are computed only for
numerically usable target/background comparisons. An RMSD difference compares
separately aligned cores: it is not a four-protein common-residue fit or a
physical displacement attributable to duplication. Different mapped residues,
coverage and prediction quality must remain part of interpretation and modeling.
Original failed strict audits are preserved; no metric is silently corrected.

Plan: `metadata/matched_domain_measurements_plan_20260927.json`.
Output: `results/structural_comparisons/matched-domain-measurements-20260927-v1`.
Producer: `scripts/prepare_matched_domain_contrasts.py`.
Endpoint-orientation fixtures passed unequal lengths, both orders, distinct
original/retained coverage and missing inputs. Full independent verification passed all 1,200,640 rows, 27,604,480 numeric
values and 14,407,680 qualification flags. All eight mask/order combinations
per match were present exactly once. There are 600,320 full-mask and 599,680
pLDDT70 numerically computable rows, with 640 pLDDT70 rows excluded.
Proof: `metadata/matched_domain_measurements_readback_20260927.json`.
Completion: `metadata/matched_domain_measurements_completed_20260927.json`.
Numerical computability alone is not scientific qualification.

Resources are one CPU, 16 GiB, no swap, 4 GiB output and 0.1–2 hours planned.
No GPU prediction or paid resources. Next steps are selected-event projection, sequence/coverage balance checks and
family/phylogeny-aware effect estimation with control reuse accounted for.
