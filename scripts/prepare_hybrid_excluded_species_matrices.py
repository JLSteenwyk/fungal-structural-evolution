#!/usr/bin/env python3
"""Remove only the two declared hybrid entries from both complete matrices."""
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    omitted={'F27292':'Saccharomyces pastorianus','F332112':'Saccharomyces cerevisiae x Saccharomyces kudriavzevii'}
    with Path('metadata/analysis_manifest.tsv').open() as f:manifest={r['taxon_id']:r for r in csv.DictReader(f,delimiter='\t')}
    assert all(manifest[t]['species_name']==name for t,name in omitted.items())
    root=Path('results/phylogeny/hybrid-excluded-matrices-20260927-v1');root.mkdir(parents=True,exist_ok=False)
    result=[];expected=None
    for label in ['profile','mafft']:
        source=Path('results/phylogeny')/(label+'-matrix-50-v1')
        data=source/'matrix.faa';records=list(SeqIO.parse(data,'fasta'))
        assert len(records)==526 and len({r.id for r in records})==526
        assert set(omitted)<=set(r.id for r in records)
        keep=[r for r in records if r.id not in omitted]
        assert len(keep)==524
        if expected is None:expected={r.id for r in keep}
        assert {r.id for r in keep}==expected
        columns=len(records[0].seq);assert all(len(r.seq)==columns for r in records)
        out=root/label;out.mkdir()
        with (out/'matrix.faa').open('w') as f:
            for r in keep:f.write('>'+r.id+'\n'+str(r.seq)+'\n')
        # Column identity and partition/site coordinates stay unchanged.
        for name in ['partitions.nex','site_mapping.tsv']:(out/name).write_bytes((source/name).read_bytes())
        result.append(dict(label=label,taxa=524,columns=columns,source_matrix=str(data),source_matrix_sha256=sha(data),source_receipt_sha256=sha(source/'receipt.json'),output_matrix_sha256=sha(out/'matrix.faa'),artifacts={p.name:sha(p) for p in out.iterdir()}))
    receipt=dict(status='complete_hybrid_excluded_full_matrices_pending_readback',omitted_taxa=omitted,retained_taxa=524,runs=result,manifest_sha256=sha('metadata/analysis_manifest.tsv'),script_sha256=sha(__file__),scope='Only two hybrid entries removed; all retained sequence strings, column order, partition/site mapping and outgroups unchanged. No realignment or new tree inference yet. Other taxon-identity/contamination sensitivities remain separate.')
    (root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
