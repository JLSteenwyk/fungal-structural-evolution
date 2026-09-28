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

## Same-process sample mapping

The installed Haskell library constructs ancestral names with
`Graph.addAncestralLabel`; its node IDs arise from a process-local ID supply.
A separately loaded tree therefore cannot establish sample identity merely
by matching numerical names. The generated fixed-tree model also omits a tree
logger. For the diagnostic chains, we insert one tree export immediately
before `makeMCMCState`, using `addInternalLabels` and `writeNewick_rooted` on
the actual model tree. No model expressions are changed. The exact original
and modified programs and their hashes are retained for every job.

`run_baliphy_sample_mapping_grid.py` runs 20 iterations for all 324 configurations,
waiting for each verified initialization receipt. Unsuccessful initializations
remain explicit dispositions. This uses two CPU workers, 32 GiB aggregate RAM,
zero swap, 12 GiB sampled RSS and 30 minutes per process, with 20 GiB storage
and a 1–82 hour planning allowance. These are short output-integrity checks,
not final posterior chains, and do not establish MCMC capacity at convergence.

`audit_baliphy_sample_mapping.py` checks the sole tree-export modification,
rooted clades and branch lengths against each specified input tree, all leaf
identities and known residues, and all sample/node associations. The 21 logged
states (iterations 0–20) and three saved alignments (0,10,20) are retained.
Candidate ancestors are identified by exact descendant sets, not numerical
node names. The first four complete chains pass all checks, providing 48
verified candidate-node/sample mappings, with 320 jobs pending in the snapshot.
The degree-two root is retained only conditional on the specified root position.
No early sample is a qualified posterior draw for biological conclusions.
