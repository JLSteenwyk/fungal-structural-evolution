# Coding-sequence acquisition and validation

The selection analyses require coding sequences tied to the same assembly, annotation and protein identifiers as the structural and orthology inputs. Acquisition has started for all **519 NCBI-backed taxa** in the 526-taxon working manifest. Every CDS URL must have the same assembly directory and prefix as its recorded protein FASTA. The prelaunch plan permits two download/validation workers, 50 GB additional disk and approximately 0.5–8 hours on the existing host; these are planning allowances, not measured completion costs. No paid resources are used.

```bash
python scripts/download_coding_sequences.py
python scripts/inventory_external_cds.py
python -m unittest discover -s tests -p test_cds_fasta_validation.py
```

For each NCBI assembly, the downloader archives the publisher checksum listing, requires exactly one checksum for the CDS file, checks MD5 before accepting the download and records local SHA256. It checks gzip integrity, unique FASTA record identifiers, nonempty sequences and the IUPAC DNA alphabet. Counts of sequences without a protein identifier, non-triplet lengths and ambiguous bases remain explicit. Those properties are flags for later interpretation, not automatic deletion criteria. Three tests confirm that partial/unlinked records remain flagged, repeated IDs are rejected and non-DNA input fails validation.

The downloader holds a single-writer lock and resumes validated sources only after checking their sequence and checksum-listing hashes against the cached receipt and exact taxon/assembly/URL identity. SIGTERM or SIGINT stops new submissions and drains the two active downloads. Download errors and incomplete taxa remain in the event log and final status. Changed cached sources or publisher checksum listings require review rather than overwriting previously validated evidence. A completed download does not establish translation correctness.

Files remain under `data/cds`; each taxon's result is appended to `data/raw/cds_downloads.jsonl`. The eventual batch summary is `metadata/cds_download_receipts.json`, with explicit pending/error taxa and the seven external-source taxa. That completion file is written when the batch ends. Execution log: `logs/cds_downloads_v1.log`. The launch plan is `metadata/cds_acquisition_resource_plan.json`.

## External-source taxa

The external-source inventory checks actual CDS file hashes, unique sequence IDs and the DNA alphabet against the existing published/extraction receipts. It identifies six available CDS sources and one extraction gap:

| Taxon | Available CDS records | Current status |
| --- | ---: | --- |
| Amoeboradix gromovi | 7,220 | Previously extracted with exact protein translation |
| Sanchytrium tribonematis | 9,367 | Previously extracted with exact protein translation; one source protein remains out of assembly bounds |
| Corallochytrium limacisporum | 7,535 | Published CDS FASTA; translation/correspondence validation pending |
| Chromosphaera perkinsii | 12,463 | Published CDS FASTA; translation/correspondence validation pending |
| Pirum gemmata | 21,835 | Published CDS FASTA; translation/correspondence validation pending |
| Abeoforma whisleri | 17,283 | Published CDS FASTA; translation/correspondence validation pending |
| Creolimax fragrantissima | — | Annotation-guided CDS extraction pending |

Sources and artifact hashes are in `metadata/external_cds_source_inventory.json`. The two sanchytrid CDS sets inherit the detailed coordinate, strand and translation audit described in the gene/isoform workflow; their provisional ORF status and the Sanchytrium exception remain explicit. Raw external source files are preserved. The available published CDS counts do not prove one-to-one protein mappings or complete gene annotations.

## Requirements before selection tests

Next, link CDS records to exact versioned protein identifiers and gene-representative decisions. Validate translations under the documented genetic code, with terminal-stop handling, partial CDSs, ambiguous residues, frameshifts and annotation-specific exceptions explicit. Creolimax requires ordered strand-aware extraction from its deposited genome and CDS annotation, with phase and translation checks. Never force translation agreement by trimming unexplained residues or picking a different code solely to improve the match.

Codon alignments must preserve validated amino-acid correspondence and gene-copy identity. Test suitability at the family/clade scale, including alignment quality and synonymous divergence, before choosing selection models. Full-genome downloads or raw sequence availability do not establish suitability for deep fungal comparisons. Selection tests, multiple-testing correction and mapping supported sites to structures remain pending; structural acceleration alone is not evidence of positive selection.
