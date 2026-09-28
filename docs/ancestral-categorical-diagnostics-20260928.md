# Categorical ancestral-state diagnostic

`scripts/ancestral_categorical_diagnostics.py` screens one source-node/extant-
residue anchor across four independently seeded chains. It receives categorical
state codes and constructs a binary indicator for each declared amino-acid,
X or gap state. It never applies a scalar diagnostic directly to arbitrary
integer amino-acid labels. The caller must verify matching model/input,
source-node and anchor identities, independent seeds and burn-in handling.

For nonconstant observed indicators it reuses the existing project's
rank-normalized R-hat, bulk/tail effective sample size and mean Monte Carlo
standard error screen. Its thresholds remain R-hat < 1.01 and both ESS values
at least 400. These diagnostics follow the methods described by
[Vehtari and colleagues](https://arxiv.org/abs/1903.08008); the
[Stan diagnostic reference](https://mc-stan.org/rstan/reference/Rhat.html)
describes rank/folded split diagnostics and approximately 100 effective draws
per chain. Applying them to individual state indicators is this project's
marginal screening design, not a theorem establishing convergence of the
joint ancestral-sequence/alignment posterior.

Constant agreement is reported as requiring review. A state not observed in
any chain has unresolved probability rather than inferred probability zero.
A state that is constant in one chain inherits the conservative review flag
from the existing scalar screen. This is an explicit project diagnostic
policy, not a claim that every constant chain must be scientifically wrong.
We must distinguish strong concentration from failure to explore states.

The report retains each state's four counts and frequencies, all six empirical
pairwise total-variation distances, and every indicator's diagnostic status.
Frequency distances are descriptive: there is no iid categorical test or
confidence interval that ignores chain autocorrelation. Passing observed
indicators is labeled `observed_state_indicator_screens_pass_only`; unseen
states remain explicit and the report never labels a posterior qualified.

The initial saved-sample horizon gives only 75 or 50 retained observations per
chain after the two burn-in cutoffs. The ESS target is not reduced to make
these short traces pass. Adequate mixing may require longer sampling. The
current inference horizon and its pinned scripts are unchanged.

Eight tests passed in the locked ancestral-diagnostic environment: independent
categorical draws with an unseen state; constant agreement; distinct stuck
modes; chain frequency shifts; slow autocorrelation; a rare state absent in
three chains; amino-acid label reordering; and invalid/insufficient draws.
The rare and constant cases cannot return an observed-indicator pass.

```bash
SOFTWARE/ancestral-diagnostics-20260927/bin/python -m unittest discover \
  -s scripts -p test_ancestral_categorical_diagnostics.py
```

This completes the tested diagnostic primitive. Full production integration
still requires provenance-checked four-chain arrays, both cutoffs, preserved
node/anchor denominators, computational planning, and all unresolved states.
No scientific chain quartet was replaced by synthetic replicates, and no
production posterior was assessed or qualified by these fixture results.

## Exact trace reuse

`scripts/ancestral_state_patterns.py` groups only byte-identical complete
chain-by-draw state traces, retaining a pattern identifier for every original
node/anchor coordinate. Chain order, temporal order, state labels and all
coordinates are preserved. Equal frequencies with different draw order remain
different patterns. Identical traces at different anchors do not establish
homology or biological independence: this is computation reuse only.

Five tests passed: full array reconstruction with duplicates; equal frequencies
but different temporal sequences; preserved chain order and state labels;
identity of reused and separately calculated categorical reports at every
fixture coordinate; and rejected malformed inputs. Coordinate multiplicity
must be retained in downstream reporting, even when computation is shared.

The seven completed production chains were each inventoried independently
at both cutoffs. Across those 14 chain/cutoff records, 943,656 node/anchor
records map to 22,129 distinct single-chain temporal patterns. Every input
state reconstructed exactly; no anchor, state, or draw was dropped. The
inventory completed in 1.50 seconds (about 70 MiB measured peak RSS).
Cutoffs overlap, so these totals are not independent biological counts.
Four-chain pattern cardinality can be larger; single-chain compression cannot
be used as a direct estimate of full-grid diagnostic cost. No artificial
quartet was assembled from the seven different input models.

Inventory evidence: `metadata/ancestral_state_pattern_inventory_20260928.json`.
Reproduce with a new receipt path:

```bash
python scripts/inventory_ancestral_state_patterns.py \
  --completion metadata/baliphy_first_anchored_states_completed_20260928.json \
  --output metadata/ancestral_state_pattern_inventory_NEW.json
```
