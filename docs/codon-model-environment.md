# Codon-model environment

HyPhy 2.5.101 is installed from commit
`646e16eb4243e5e5588ac340483fcb987552751d` in the shared SOFTWARE directory.
CMake 3.31.6 resides in an isolated build environment; GCC 13.3 builds the
executable with four workers. System BLAS and ZLIB are available; CURL and MPI
were absent. Analysis commands will use local input files and explicit CPU limits.

The exact configure/build/install commands and paths are recorded in
`metadata/hyphy_build_config.json`; source retrieval is
`git clone --depth 1 --branch 2.5.101 https://github.com/veg/hyphy.git`, followed
by verification of the full commit above. The executable and 457 installed-file
checksums are recorded in the installation receipt and its referenced manifest.

The pinned upstream mitochondrial codon test passed its likelihood, dN/dS and
tree-length assertions using two CPU threads. See
`metadata/hyphy_numerical_test_receipt.json`. The installation receipt records
the earlier installation checkpoint; this separate receipt establishes the
subsequent test completion. This is one numerical test, not exhaustive validation.

For method configuration, use the [pinned BUSTED source](https://github.com/veg/hyphy/blob/646e16eb4243e5e5588ac340483fcb987552751d/res/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf):
`--multiple-hits`, `--mss` and `--syn-rates` control distinct options. This version
rejects MSS combined with multiple hits or branch-site synonymous-rate variation.
Do not substitute generic web examples for the selected version's CLI. The
[official build instructions](https://github.com/veg/hyphy/blob/646e16eb4243e5e5588ac340483fcb987552751d/README.md)
describe local installation.

Fungal selection tests have not started. Existing genus/code coverage screening
is insufficient by itself: gene-copy identity, topology, alignment and synonymous
divergence remain analysis-specific requirements. Preserve mixed-code splits,
recorded annotation/taxon exceptions and alternative alignments, and assess
recombination and model adequacy before interpreting selection tests.
