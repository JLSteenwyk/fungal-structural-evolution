# Full matched-model workload and exact input reuse

The working linear model grid contains 82,944 record-level settings across
192 measurement/qualification settings and 432 guide/policy/scenario strata.
Using all five species trees gives 414,720 fits before exact reuse. At 22
boundary/start candidates per fit, that is 9,123,840 optimizer attempts. This
is a workload count, not a wall-time forecast or a statement that all attempts
converge. Verified record summaries range from 148 to 10,963 records per fit;
there are 266,729,088 record occurrences across the 82,944 settings.

The inventory retains a row for every setting and seeks only exact duplicate
likelihood inputs. Each input is ordered by target node ID and represented by:

- Target and background node identities, family-component label and species
  pattern identity for every record.
- Domain-averaged RMSD response and the four candidate covariate differences.
- Record count, numeric-column order and the explicit working model specification.

Ordered identities and little-endian 64-bit numeric arrays are SHA-256 hashed.
Signed zero is canonicalized; values are otherwise never rounded or matched
within a tolerance. Node identities are retained even if observations happen
to share numerical measurements. Settings with different rows or identifiers
therefore remain distinct. Each unique input has a source partition/stratum
recipe, and the complete mapping preserves every original scenario, policy,
boundary, mask, cohort, threshold and input order. The five trees remain
separate fit alternatives.

This inventory applies only to the specified equal-record linear Gaussian
working model. Changing covariates, response construction, dependence terms,
weighting or optimization specification requires reassessing reuse. It is not
a mechanism for merging biological replicates or declaring settings equivalent
under every future model.

Producer: `scripts/inventory_matched_model_fits.py`.
Independent checker: `scripts/readback_matched_fit_inventory.py`.
Plan: `metadata/full_matched_fit_inventory_plan_20260927.json`.
Launches: `metadata/full_matched_fit_inventory_launch_20260927.json` and
`metadata/full_matched_fit_inventory_readback_launch_20260927.json`.
Output: `results/model_validation/full-matched-fit-inventory-20260927-v1`.
Expected proof: `metadata/full_matched_fit_inventory_readback_20260927.json`.

The checker waits for the exact producer process and successful authoritative
terminal state, then reconstructs every fingerprint using independent SQL
joins and identity construction. Every setting label, source recipe and record
denominator must agree. Both stages use one CPU, 12 GiB memory and no swap, with
0.1–2 hours planned per stage and at most 1 GiB producer output. At launch the
full inventory and readback are pending. Full-grid fitting, runtime estimates,
model adequacy and inferential calibration remain outstanding. GPU prediction
remains paused.
