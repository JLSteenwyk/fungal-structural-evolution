# Whole-protein/domain opposition: inspection cases

All 13 event/domain cases opposing at the descriptive 0.1 Å margin in any of the
six coverage screens are retained. Selection did not use function or apparent
novelty. The dataset includes all 234 case/screen/margin outcomes, 13 whole-model
triads (416 fit alternatives) and 26 domain-interval triads (832 alternatives).
Every fit passes n30/c50 and the unique-geometry check; stricter-screen failures
remain in the sensitivity table. Exact reference choices are retained.

Pfam 38.2 labels below are computational family annotations, not verified
activities of the particular proteins. Coverage is not a calibrated uncertainty
measure. Signed ranges are RMSD(A, reference) minus RMSD(B, reference), in Å,
over every included alternative; they are not confidence intervals.

| Family | Taxon | Pfam label | Highest opposing coverage | Whole contrast range | Domain contrast range |
| --- | --- | --- | ---: | ---: | ---: |
| OG0000054 | Heliocybe sulcata | GST_C | 70% | -0.356 to -0.113 | 0.162 to 0.165 |
| OG0000095 | Rhodonia placenta | FAD_binding_8 | 50% | 0.208 to 1.073 | -0.809 to -0.196 |
| OG0000107 | Cryoendolithus antarcticus | WHD_MCM6 | 50% | -0.638 to -0.395 | 0.140 to 0.621 |
| OG0000152 | Scedosporium apiospermum | ADH_zinc_N | 50% | -0.498 to -0.211 | 0.267 to 0.609 |
| OG0000230 | Furculomyces boomerangus | Acyl-CoA_ox_N | 70% | -1.214 to -0.221 | 0.141 to 0.190 |
| OG0000294 | Jaapia argillacea | GTP_EFTU | 70% | -0.879 to -0.695 | 0.162 to 0.293 |
| OG0000336 | Synchytrium microbalum | FA_hydroxylase | 50% | 0.169 to 0.312 | -0.105 to -0.105 |
| OG0000972 | Neolecta irregularis | FA_desaturase | 70% | 0.300 to 0.466 | -0.151 to -0.140 |
| OG0001082 | Leucosporidium creatinivorum | HHH_8 | 50% | 0.258 to 0.505 | -0.140 to -0.140 |
| OG0001200 | Smittium simulii | LeuA_dimer | 50% | -0.243 to -0.169 | 0.107 to 0.541 |
| OG0001203 | Leucosporidium creatinivorum | GST_N_3 | 50% | -0.891 to -0.700 | 0.279 to 0.360 |
| OG0002650 | Piloderma croceum | rva_4 | 70% | 2.488 to 2.587 | -0.160 to -0.154 |
| OG0002812 | Phycomyces blakesleeanus | Meth_synt_1 | 90% | 0.136 to 0.151 | -0.130 to -0.129 |

## Next discriminating analyses

The *Phycomyces* case (OG0002812) is the only case retained at 90% coverage.
Its approximately +0.14 Å whole-protein and −0.13 Å domain contrasts are small;
independent prediction settings or experimental structures are needed before
interpreting that persistence biologically. It is a stringent-coverage control,
not automatically the strongest functional candidate.

For all 13 cases, the next coordinate analysis should fit the matched domain and
measure displacement in the remaining matched residues under that same transform.
Compare this with independently fitted domain and whole-protein residuals, and
retain confidence masks, mappings, boundaries and reference choices. This can
distinguish localization outside the domain from changes within it; direct
orientation claims additionally need confidently mapped second domains and
interdomain uncertainty checks.

The Piloderma case (OG0002650) has a large whole-protein contrast but a much
smaller reversed domain contrast. Its integrase-core Pfam label makes gene-model,
repeat/mobile-element and duplication-context review necessary before choosing
an experimental mechanism. Annotation alone does not establish a host cellular
function. The other cases also require sequence/model provenance and architecture
review before selecting biochemical assays.

If coordinate and independent-prediction checks support changes in relative
domain arrangement while domain cores remain similar, experimentally testable
hypotheses can focus on interdomain contacts or linker substitutions. If the
changes instead localize within a conserved domain, tests should target those
residues, with activity and substrate choices informed by independently supported
function. These are conditional experimental proposals, not results.

## Reproduction

Run `scripts/prepare_whole_domain_case_dossiers.py` using one CPU. Output creation
is exclusive. Inputs are hash checked, all required triads have exactly 32 fit
alternatives, all selected fits pass the declared baseline, and every output
cell passes serialized readback. The full input fits were independently checked
in their production stages.

Outputs are in
`results/structural_comparisons/whole-domain-case-dossiers-20260927-v1`:
case annotations and extrema, all sensitivity settings, all linked references,
and full whole/domain fit alternatives. Hashes and provenance are recorded in
`metadata/whole_domain_case_dossiers_completed_20260927.json`. GPU inference
remains paused; independent prediction experiments have not been run for these
cases. See [full integration](whole-domain-contrast-integration-20260927.md).


## Integrated ancestral uncertainty for every structural case

The [updated 13-case table](tables/case_ancestral_uncertainty_20260927.tsv)
preserves all 60 original structural and annotation fields and appends 15
ancestral uncertainty fields. These cover optimization sensitivity, whole/domain
state disagreement with its denominator, opposing high-confidence calls with
and without local domain coverage, and the number of exact coordinate and
source-node groups underlying repeated calls. Every original field was checked
against the serialized joined table; no case was removed or reordered.

The table supports case review, not a biological ranking. All cases retain an
explicit statement that their current amino-acid diagnostics are conditional
and joint indel/convergence requirements are unresolved. Zero high-confidence
conflict does not qualify a case for ancestral structure prediction. The
reproducible join is `scripts/integrate_case_ancestral_uncertainty.py`; input and
output hashes are recorded in
`metadata/case_ancestral_uncertainty_integration_completed_20260927.json`.
# Combined experimental and ancestral evidence

The [complete case-evidence table](tables/case_evidence_all_settings_20260927.tsv)
joins each original case dossier and its ancestral-uncertainty fields to all
experimental-reference sensitivity settings. It retains 13 cases × 72 settings
(four metrics, six coverage screens, three margins), with 90 fields per row.
All 936 rows passed a separate relational-merge check of every field; 648
settings have no qualifying experimental reference units and remain explicit.

Case and ancestral fields repeat across settings and must not be summed.
Reference units share structures, sequences and publication/dependence groups;
their counts are not independent evolutionary replicates. This table supports
case review without selecting a favorable screen. It supplies no ranking,
significance test, ancestral polarity or qualified ancestral structural ensemble.
Reproduce with `scripts/integrate_case_experimental_ancestral_evidence.py`;
source checksums and validation scope are recorded in
`metadata/case_experimental_ancestral_integration_completed_20260927.json`.
