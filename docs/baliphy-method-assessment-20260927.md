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
