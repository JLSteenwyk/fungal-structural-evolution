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
