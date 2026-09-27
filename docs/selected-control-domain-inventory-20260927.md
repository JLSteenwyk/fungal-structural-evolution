# Domain comparisons available for selected controls

All 2,786,912 selected metadata-control records have been linked to the
single-copy Pfam-domain comparison inventories. The full independent verifier
checked every original selection field, target/background node identity,
policy, source protein pair, domain accession, both boundary definitions,
configuration key, and aggregate. Producer and verifier terminated successfully.

The 104,918 distinct target-protein-pair/background-protein-pair/policy
configurations comprise 60,828 with shared domain comparisons, 36,584 with no
shared comparison in these inventories, and 7,506 involving identical models
that require separate handling. Configurations are reused across guides,
selection scenarios and biological records; they are not independent events.
No shared comparison does not establish biological domain absence.

Only primary duplicate-pair links (including primary+reference roles) enter
the target side. Reference-only comparisons are excluded from that role.
Every shared Pfam accession retains both alignment and envelope boundaries
and its exact target/background domain-pair keys. Original selection rows
remain intact in `selection_domain_links.tsv.gz`; repeated calculations are
represented once in `domain_configurations.jsonl` without discarding their
many-to-one associations. Identical models receive no invented zero distance.

Output: `results/structural_comparisons/selected-control-domain-inventory-20260927-v1`.
Source plan: `metadata/selected_control_domain_inventory_plan_20260927.json`.
Full proof: `metadata/selected_control_domain_inventory_readback_20260927.json`.
Scripts: `inventory_selected_control_domains.py` and
`readback_selected_control_domains.py`. Each used one CPU, 16 GiB and no swap;
the prelaunch budget was 4 GiB output and 0.1–2 hours. No predictions or charges.

Next, join these configurations to the audited coverage/geometry tables and
retain all scenarios, masks, boundary alternatives and missingness categories.
That will establish the actual matched domain cohort before effect estimation.
Outcome availability, repeated controls, family and phylogenetic dependence,
and prediction uncertainty still require explicit treatment. No effect is
estimated by this inventory.
