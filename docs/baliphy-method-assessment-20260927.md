# BAli-Phy assessment for ancestral uncertainty

The [author guide](https://www.bali-phy.org/README.html) describes joint samples
of ancestral residues and gaps, matched to sampled trees and alignments.
It supports fixed topology or fixed tree analyses. Internal node identifiers
must be interpreted through each matching tree. These capabilities address
uncertainty missing from the current Historian outputs; suitability for our
largest families and MCMC convergence remain untested.

The official [4.3 release](https://github.com/bredelings/BAli-Phy/releases/tag/4.3)
was installed locally using its Ubuntu 24.04 archive. The archive matches the
release SHA-256. Version/help checks pass; no project MCMC was launched.
`install_baliphy_local.py` and its versioned plan reproduce installation without
system changes. Files and help captures are under
`data/software_audits/baliphy-4.3-20260927`.

Before inference, specify priors and indel-rate assumptions, verify all-family
input preservation and capacity, and establish multiple-chain convergence and
sampling diagnostics. Existing Historian checks continue as method sensitivity.

## Initialization checks across the complete selected input grid

`run_baliphy_initialization_grid.py` runs `--test` for all 324 derived-tree
configurations, covering all 78 alignments and retaining all protein copies.
The two floors and all 15 resolutions remain computational sensitivities.
This command initializes a model and exits; it does not sample a posterior.
Resources: two CPU workers, 16 GiB aggregate RAM, no swap, 6 GiB sampled RSS
and 10 minutes per process. The full planning allowance is 0.1–28 hours and
5 GiB output. No GPUs or paid resources are used.

The explicit model arguments are `LG +> ASRV.Gamma(4,alpha=1)`,
`RS07(rate=0.01,meanLength=3)`, constant indel rates, a fixed tree and seed
20260927. The generated Haskell code reveals an additional default: LG
exchangeabilities receive a sampled frequency vector under a symmetric
Dirichlet(1) prior. These frequencies are not the published LG frequencies
or the IQ-TREE empirical estimates. Original internal labels are dropped.
Source data are converted to unaligned character data; source column numbers
must not be used as posterior coordinates. These behaviors are retained and
reported, not silently equated with the preceding model fits.

`audit_baliphy_initializations.py` checks generated model code, input paths,
finite initial scores, prior-plus-likelihood arithmetic, normalization of all
20 frequencies and agreement of total tree length. Its first snapshot verifies
18 jobs and retains 306 pending dispositions. This does not independently
replay likelihoods, verify runtime ancestral identities, establish sample
preservation, measure MCMC mixing or show convergence. Before posterior use,
samples must be linked to the software's internal tree identities and validated
against exact extant sequence identities and candidate descendant sets.
