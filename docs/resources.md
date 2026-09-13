# Resources and cost model

Observed at initialization: 192 logical CPUs, approximately 1 TiB RAM, 12 TiB available on workspace volume, two NVIDIA RTX 6000 Ada GPUs (~48 GiB each). One GPU had ~21 GiB in use. These are capacity observations, not reserved allocations.

Planning scenarios only: 525 proteomes × 8,000/12,000/20,000 proteins gives 4.2/6.3/10.5 million proteins. Actual counts and sequence lengths must replace these assumptions after inventory. Outgroups may differ substantially.

Prediction GPU-hours = missing protein count × measured seconds/protein / 3600. At 6.3 million missing proteins, illustrative 1/10/60 seconds per protein implies 1,750/17,500/105,000 GPU-hours. These are sensitivity scenarios, not benchmark claims. Wall time also depends on GPU allocation, length distribution, MSA generation, retries and utilization.

At 6.3 million proteins, illustrative 50/200 kB per stored structure implies 315 GB/1.26 TB for coordinates alone (decimal units). Budget separately for sequence databases, MSAs, intermediates, search indexes, ensemble predictions and backups. Full all-pairs search is impractical; use scalable clustering and targeted within-family comparisons.

No paid resources provisioned and no prediction jobs launched. Before expensive stages: measure full-dataset counts, inventory existing predictions, assess tool/database requirements and available allocation, then publish updated estimates. Bounded throughput measurements during production are QC, not a separate biological pilot.

Full-scale QC acquisition authorized on existing disk: 500 selected compressed proteomes, conservatively budgeted at up to 10 GB total before observed counts replace the scenario. Two concurrent download workers; no paid services or GPU prediction jobs. Publisher checksums and FASTA parsing occur as each download completes. `metadata/proteome_qc_snapshot.json` tracks measured bytes, residues and protein counts for completed files; not a prediction-throughput benchmark.

External supplements: seven Figshare articles provide two fungal annotations/genomes plus five candidate outgroup genome bundles; total retrieval budget below 1 GB. BUSCO 6.1.0 installation is isolated under ignored .cache/envs/busco. Environment specification follows https://busco.ezlab.org/busco_userguide.html; installation success and package lock must be verified before biological runs.

Initial completeness batch: BUSCO protein mode, eukaryota_odb12.2 (125 markers), 4 concurrent taxa × 4 threads = maximum 16 configured BUSCO threads. Scenario of 1–20 minutes per taxon implies roughly 2–42 hours wall time for 502 taxa at four-way concurrency; this is a planning range, not measured throughput. Existing CPU/disk only. Subsequent batches will reuse verified success receipts, and more specific lineage assessments remain required.
