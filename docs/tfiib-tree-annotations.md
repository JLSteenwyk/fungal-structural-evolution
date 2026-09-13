# TFIIB family-tree input annotations

The exact 1,040 combined-domain tree tips (512 taxa) are now joined to source
protein/gene identities, domain-hit counts, BUSCO marker selection, taxonomy and
curated hybrid flags. The same tip identities underlie both repeat-specific
alignments. Receipts pin all inputs, and all 1,040 triplets of TFIIB, BRF1 and
Zn_Ribbon_TF hit counts were independently checked against source domain-hit rows.

There are 473 BRF1-detected and 567 BRF1-undetected input entries; 446 taxa have
both classes in these tree inputs. Among tree-eligible entries selected for
BUSCO marker 4986044at2759, 431 are BRF1-detected and 33 are BRF1-undetected.
Thirteen of those 33 taxa also have an aligned BRF1-detected candidate. These
counts extend the initial Aspergillus copy-selection concern to the full family
input universe. They are not a genome-wide paralogy estimate: domain-pair
eligibility excluded other candidates, and profile non-detection is not proof
of domain loss or of protein function.

The combined and repeat trees are still computing. Completed trees must retain
expected tip identities and valid support records before testing how these
evidence classes map onto clades. Tree topology, species-tree reconciliation,
alignment sensitivity and gene/subgenome assignment remain necessary to resolve
copy selection. Hybrid flags do not identify parental origin of individual tips.
No number or timing of ancestral duplications is inferred from these annotations.

```bash
python scripts/annotate_tfiib_tree_inputs.py \
  --output results/phylogeny/tfiib-tip-annotations-v1
```

Full annotation tables and source receipts are versioned under `metadata/`.
The classes deliberately use detected/undetected terminology and do not assign
experimentally validated TFIIB/BRF1 function from domain counts alone.


## Completed combined-domain tree and exploratory class split

The combined 184-column tree completed for all 1,040 tips with 1,000 bootstrap
trees. Independent split analysis checked all bootstrap tip identities and
finite/nonnegative branch lengths and parsed 1,037 internal splits. The 48
unlabeled internal groups each contain a single distinct aligned sequence;
no support is imputed for their identical-tip expansions.

The split minimizing disagreement with the domain-detection classes has
SH-aLRT/UFBoot support 100/100 and occurs in all 1,000 saved bootstrap trees:

| Side of unrooted split | BRF1 detected | BRF1 undetected | Selected focal BUSCO copies |
| --- | ---: | ---: | ---: |
| BRF1-enriched | 473 | 12 | 437 |
| Other side | 0 | 555 | 27 |

This is an exploratory split selected using annotations, not a prespecified
association test. Bootstrap support concerns the split, not protein function
or independence of duplication events. Nevertheless, the original marker
selection draws copies from both deeply separated parts of this family tree.
`config/marker_orthology_review.json` records a confirmatory-analysis caveat for
4986044at2759: gene-copy reconciliation and repeat-tree sensitivity are required
before single-ortholog rate interpretation. The unrooted branch does not date
or locate the ancestral duplication, and BRF1 non-detection is not domain loss.

```bash
python scripts/assess_tfiib_domain_class_split.py \
  --output results/phylogeny/tfiib-class-split-v1
```

Two tests check split invariance to root placement and rejection of duplicate
tips/negative branches. Complete split counts, exact tip-side annotations and
receipts are versioned. The separate repeat trees are still computing.


## Effect on current paired structural analyses

The caveat is now attached to all 250 marker/source review rows (125 markers in
each source), with no implied validation of unflagged markers. In the expanded
AlphaFold input, marker 4986044at2759 is eligible with 116 taxa: exact
 taxon/protein/sequence joins place 108 on the BRF1-enriched family-tree side,
four on the other side and four outside the eligible family-tree tip set.
The local snapshot has one eligible taxon on the other side and fails the
four-taxon marker gate. The earlier 52-marker completed fits did not include
this marker. Thus the active expanded AlphaFold input itself exhibits mixed
copy selection; this is not merely a concern inferred from unrelated genomes.

`metadata/paired_marker_review.tsv` records
`withhold_pending_copy_reconciliation` for this marker. Existing numerical fits
remain diagnostic and their immutable inputs are preserved. No automatic
pruning or reclassification of the four unresolved candidates is performed.
Downstream confirmatory summaries must consume the explicit review status.

```bash
python scripts/apply_marker_review_caveats.py \
  --output results/phylogeny/marker-review-overlay-v1
```
