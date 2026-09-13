# Domain annotation workflow

The first production layer annotates all 59,840 complete marker-protein records spanning the 526-taxon extraction. Exact sequence deduplication yields 58,883 search targets (33,047,630 residues), with all marker, taxon, protein and source checksums retained in a link table. Deduplication reduces repeated search work; it does not collapse species-specific records for evolutionary analysis. Full representative-proteome annotation remains required and is a subsequent stage.

## Pinned resources and execution

Pfam 38.2 contains 30,134 families and is based on UniProtKB 2025_03 according to its version file. The [release-specific EBI archive](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/) supplies profile HMMs, entry types, clan assignments and active-site metadata. `prepare_pfam.py` pins publisher MD5 values, verifies decompression and records local SHA256 hashes. Active-site metadata are acquired for later use; no site transfer is implemented or claimed here. All HMMs must have unique accessions, match-state lengths and gathering thresholds. Their total count must agree with the release. Profiles are assigned deterministically to 64 chunks balanced by summed match-state lengths.

The environment is specified in `environments/domains.yml`; the current installed HMMER binary is additionally pinned by hash at execution. HMMER 3.4 `hmmsearch --cut_ga --cpu 2 --noali` searches each HMM chunk against the entire deduplicated marker database. Four jobs run concurrently. E-values refer to this unique-sequence target database, not all fungal proteins; Pfam-specific sequence and domain gathering scores control acceptance. See the [Pfam glossary](https://pfam-docs.readthedocs.io/en/latest/glossary.html) for the two curated gathering thresholds. Search commands, binary hash/version, source hashes, runtimes and output hashes are recorded. Resume accepts only matching completed receipts; unreceipted outputs require inspection after interruption. A lock prevents concurrent writers.

```bash
python scripts/prepare_pfam.py
python scripts/prepare_marker_domain_inputs.py
python scripts/run_marker_domains.py --output results/domains/marker-search-v1
python scripts/summarize_marker_domains.py --search results/domains/marker-search-v1 --output results/domains/marker-annotations-v1
```

Marker inputs and interpreted annotations are immutable snapshots; prepare them in a new versioned location when inputs change. The search can resume its existing output directory under an unchanged configuration. Raw files remain under ignored `data/pfam/`, `data/domains/` and `results/domains/`. Published metadata receipts contain retrieval and verification information.

## Interpretation and downstream gates

Only a complete search across all pinned profiles can enter the annotation summary. The parser verifies target lengths, model identities and lengths, inclusive alignment/envelope coordinates, finite scores and gathering-threshold consistency (allowing 0.051 bits for printed-score rounding). It preserves Pfam versioned accession, entry type, clan, HMM coverage, scores, E-values and coordinates. Types such as Family, Repeat and Domain remain distinct; a Pfam match is not automatically a discrete structural domain or an orthogroup.

All gathering-threshold matches are retained. Overlaps in alignment coordinates, including a shared boundary residue, are explicitly recorded with same-clan and same-family flags. This is not a resolved domain architecture. Overlap resolution, fragmented matches, nested domains, repeated domains, alignment-boundary uncertainty and independent structural support must be assessed before counting gains, losses, fusions or duplications. Unknown or undetected domains are not treated as confirmed absences. Preserving the raw hits supports later comparisons of alternative overlap policies.

The parser and overlap-coordinate tests are implementation checks, not biological validation. Domain-restricted structural comparisons, domain architecture reconstruction, functional-site transfer and full-proteome annotation remain pending until their execution is separately documented.
