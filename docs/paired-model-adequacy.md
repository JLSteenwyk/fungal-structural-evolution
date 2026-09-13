# Paired-marker model diagnostics

Completed unfiltered IQ-TREE symmetry screens for all 72 local-source markers
in both alphabets (144 alignments). These diagnostics complement the existing
composition/missingness warnings; they do not establish overall model adequacy.
The [IQ-TREE documentation](https://iqtree.github.io/doc/Assessing-Phylogenetic-Assumptions)
and [Naser-Khdour et al. 2019](https://doi.org/10.1093/gbe/evz193) describe the
maximum matched-pairs tests. Original matched observation masks and identical
taxa were retained. No marker was removed based on these results.

| Alphabet | Symmetry maximum test available | Marginal maximum test available | Internal maximum test available |
| --- | ---: | ---: | ---: |
| AA | 72/72 | 44/72 | 44/72 |
| 3Di | 72/72 | 7/72 | 7/72 |

In total, 186 of 432 maximum-test results are nonfinite/unavailable. Undefined
results are not failures to reject a model and cannot be counted as evidence
of adequacy. The very limited marginal/internal availability in 3Di is a
constraint on this screen. Sparse contingency tables and context-dependent
characters require attention before inferential use. Even finite results are
nominal: no multiple-testing correction or simulation calibration has been
performed, and no biological rejection decisions are made here.

## Validation

Every raw result, input/output hash, alignment dimension and summary field was
checked. A separate implementation reconstructed each pair's Bowker statistic
using observed off-diagonal state counts and a chi-square tail probability.
Across all 533,186 taxon pairs, the resulting nominal significant/nonsignificant
counts matched IQ-TREE. Zero-degree pairs are counted separately rather than
assigned a p-value. Maximum-test p-values and the marginal/internal statistics
were not independently recalculated. Matching software arithmetic does not
validate asymptotic calibration for these data.

The raw test counts and maximum-test p-values remain in
`metadata/esmfold_paired_symmetry_raw.tsv`, with unavailable/tested-pair accounting
in `metadata/esmfold_paired_symmetry_availability.tsv`. Receipts and configurations
use the same metadata prefix. Production outputs remain outside Git.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/screen_paired_symmetry.py --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --fit-config results/phylogeny/paired-marker-fits-esmfold-partial-v1/config.json --output results/phylogeny/paired-symmetry-esmfold-v1
OPENBLAS_NUM_THREADS=1 python scripts/audit_paired_symmetry.py --screen results/phylogeny/paired-symmetry-esmfold-v1 --inputs results/phylogeny/paired-inputs-esmfold-partial-v1 --output results/phylogeny/paired-symmetry-audit-esmfold-v1
```

The prelaunch resource plan used four single-thread jobs, 8 GB RAM and 2 GB
disk headroom on the existing host. Further work includes rate-heterogeneity
sensitivity, simulation-based adequacy checks that respect missingness and
character dependence, and larger-source comparisons. Nonstationary modeling
and predictor controls remain separate requirements.
