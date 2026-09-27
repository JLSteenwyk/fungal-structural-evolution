# Duplication event/domain robustness across references and annotation alternatives

This stage links interval-triad robustness back to its full event context.
It retains all 34,909 provisionally referenced guide/event records, including
28,058 without a comparable common-domain candidate. The 6,851 events with
candidates yield 8,463 event–Pfam combinations and 50,778 rows across six
coverage screens. The two guides are sensitivity alternatives, not independent
replication; multiple Pfams within an event are also not independent events.

For each guide/event/Pfam, the required universe includes every retained tied
reference, all four annotation policies and both domain-boundary choices.
An event/domain direction label requires a domain match in every such unit
and successful coverage/geometry screening in **all 32** fitting alternatives
within every unit. Missing units and incomplete units remain explicit.
The summary reports observed/expected/missing/complete unit counts and the
range of eligible signed contrasts. Ranges in incomplete rows describe only
the available subset and cannot support a consistent-direction label.

Pfam candidates remain separate rather than averaging across domains. The
original whole-protein sequence-tree covariates and chosen-reference identity
are copied exactly with lexical gene A/B orientation. No sequence rate is
re-estimated for a domain, and no tied reference receives another reference's
distance. The whole-protein covariates are context for subsequent analysis,
not an adjustment already performed here.

## Independently verified results

At the 30-residue/70%-original-interval-coverage screen and descriptive 0.1 Å
margin:

| Guide | A farther throughout | B farther throughout | Extends beyond both directions | Touches/within margin band | Incomplete |
|---|---:|---:|---:|---:|---:|
| MAFFT | 487 | 481 | 97 | 2,573 | 592 |
| Profile | 487 | 481 | 97 | 2,574 | 594 |

Thus 968 event/domain combinations per guide retain one direction beyond
0.1 Å across every retained alternative. Matching aggregate counts do not
establish identity agreement between guides. The labels describe similarity
to provisional references, not post-duplication direction of evolution,
selection or a statistically established duplication effect. The 0/0.01/0.1 Å
margins are descriptive sensitivities and not prediction-error calibrations.
All six screens and all margins remain available in the full output.

The independent dataframe checker reconstructed all source joins, the complete
reference/policy/boundary universe, every event availability row, original
sequence covariate, extremum, completeness flag and direction label. See the
[completion record](../metadata/duplication_domain_event_robustness_completed_20260927.json)
and [pinned plan](../metadata/duplication_domain_event_robustness_plan_20260927.json).

```bash
python scripts/summarize_domain_event_robustness.py --plan metadata/duplication_domain_event_robustness_plan_20260927.json
python scripts/check_domain_event_robustness.py --plan metadata/duplication_domain_event_robustness_plan_20260927.json --output metadata/duplication_domain_event_robustness_completed_20260927.json
```

Use new output paths when rerunning; artifacts refuse overwrite. This stage
only joins existing tables and uses no GPUs or external paid resources.
Cross-guide identity checks, phylogenetic/family dependence, sequence-divergence
adjustment, matched controls and prediction-source sensitivity remain required
before biological inference.
