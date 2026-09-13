# Resources and cost model

Observed at initialization: 192 logical CPUs, approximately 1 TiB RAM, 12 TiB available on workspace volume, two NVIDIA RTX 6000 Ada GPUs (~48 GiB each). One GPU had ~21 GiB in use. These are capacity observations, not reserved allocations.

Planning scenarios only: 525 proteomes × 8,000/12,000/20,000 proteins gives 4.2/6.3/10.5 million proteins. Actual counts and sequence lengths must replace these assumptions after inventory. Outgroups may differ substantially.

Prediction GPU-hours = missing protein count × measured seconds/protein / 3600. At 6.3 million missing proteins, illustrative 1/10/60 seconds per protein implies 1,750/17,500/105,000 GPU-hours. These are sensitivity scenarios, not benchmark claims. Wall time also depends on GPU allocation, length distribution, MSA generation, retries and utilization.

At 6.3 million proteins, illustrative 50/200 kB per stored structure implies 315 GB/1.26 TB for coordinates alone (decimal units). Budget separately for sequence databases, MSAs, intermediates, search indexes, ensemble predictions and backups. Full all-pairs search is impractical; use scalable clustering and targeted within-family comparisons.

No paid resources provisioned and no prediction jobs launched. Before expensive stages: measure full-dataset counts, inventory existing predictions, assess tool/database requirements and available allocation, then publish updated estimates. Bounded throughput measurements during production are QC, not a separate biological pilot.
