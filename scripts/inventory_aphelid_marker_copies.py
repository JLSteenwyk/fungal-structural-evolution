"""Retain all BUSCO hit proteins and join existing gene-family partitions."""
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
from Bio.SeqIO.FastaIO import SimpleFastaParser
from readback_whole_proteome_catalog import sha
from ancestral_chain_attempt import write_json


def table(path):
    with open(path) as h:return list(csv.DictReader(h,delimiter='\t'))


def main():
    taxon='F1243177';base=Path('results/qc/busco-gene-copies-v1')
    receipt=json.loads((base/'receipt.json').read_text())
    marker_path=base/'marker_gene_copies.tsv';assert sha(marker_path)==receipt['artifacts'][marker_path.name]
    selected=[r for r in table(marker_path) if r['taxon_id']==taxon];assert len(selected)==125
    expected={p for r in selected for p in json.loads(r['protein_ids_json'])}
    raw=Path('results/busco')/taxon/'run_eukaryota_odb12.2/full_table.tsv'
    binding=next(r for r in receipt['sources'] if r['taxon_id']==taxon);assert sha(raw)==binding['busco_table_sha256']
    raw_ids={line.split('\t')[2] for line in raw.read_text().splitlines() if line and not line.startswith('#') and line.split('\t')[1]!='Missing'}
    assert raw_ids==expected
    decisions_path=Path(f'results/gene_representatives/full-v2/{taxon}.decisions.tsv')
    assert sha(decisions_path)==binding['representative_decisions_sha256']
    decisions={r['protein_id']:r for r in table(decisions_path)}
    fasta=Path('data/qc_proteomes')/(taxon+'.faa');sequences={}
    with fasta.open() as h:
        for title,seq in SimpleFastaParser(h):
            protein=title.split()[0]
            if protein in expected:
                assert protein not in sequences;sequences[protein]=seq
    assert set(sequences)==expected
    dbpath=Path('results/domains/family-domain-bridge-v1/family_domain_bridge.sqlite')
    br=json.loads((dbpath.parent/'receipt.json').read_text())
    assert sha(dbpath)==br['bridge_sha256']
    audit_path=Path('metadata/family_domain_bridge_completed_readback.json')
    audit=json.loads(audit_path.read_text())
    assert audit['status']=='passed_complete_independent_family_domain_bridge_readback'
    assert audit['producer_receipt_sha256']==sha(dbpath.parent/'receipt.json')
    db=sqlite3.connect('file:'+str(dbpath.resolve())+'?mode=ro',uri=True)
    guides={r[0] for r in db.execute('select distinct guide from families')};assert len(guides)==2
    records=[];families=[]
    for marker in selected:
        for protein in json.loads(marker['protein_ids_json']):
            decision=decisions[protein];genes=json.loads(decision['gene_ids_json']);seq=sequences[protein]
            assert len(seq)==int(decision['protein_length'])
            native=db.execute('select native_gene_id from native_proteins where taxon_id=? and protein_id=?',(taxon,protein)).fetchall()
            assert len(native)<=1
            memberships=[item for guide in sorted(guides) for item in db.execute('select guide,family from assignments where guide=? and native_gene_id=?',(guide,native[0][0])).fetchall()] if native else []
            assert len({x[0] for x in memberships})==len(memberships)
            records.append(dict(taxon_id=taxon,marker=marker['marker'],busco_status=marker['busco_status'],protein_id=protein,
                gene_ids_json=json.dumps(genes),representative_decision=decision['decision'],length=len(seq),
                sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),native_gene_id=native[0][0] if native else '',
                family_guides_found=len(memberships)))
            by=dict(memberships)
            for guide in sorted(guides):families.append(dict(marker=marker['marker'],taxon_id=taxon,protein_id=protein,guide=guide,family=by.get(guide,''),status='assigned' if guide in by else 'no_family_assignment'))
    assert len(records)==sum(int(r['raw_hit_proteins']) for r in selected)
    assert len({(r['marker'],r['protein_id']) for r in records})==len(records)
    out=Path('results/ecology/aphelid-marker-copy-inventory-20260928-v1');out.mkdir(parents=True,exist_ok=False)
    for name,rows in [('marker_protein_copies.tsv',records),('family_assignments.tsv',families)]:
        with (out/name).open('w') as h:
            w=csv.DictWriter(h,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    with (out/'all_hit_proteins.faa').open('w') as h:
        for protein,seq in sorted(sequences.items()):h.write('>'+taxon+'_'+protein+'\n'+seq+'\n')
    with (out/'all_hit_proteins.faa').open() as h:
        back={title:seq for title,seq in SimpleFastaParser(h)}
    assert back=={taxon+'_'+p:s for p,s in sequences.items()}
    result=dict(status='complete_all_aphelid_marker_copy_family_inventory',taxon_id=taxon,marker_universe=125,
        markers_with_hits=sum(bool(json.loads(r['protein_ids_json'])) for r in selected),marker_protein_rows=len(records),unique_hit_proteins=len(sequences),
        sequence_unique_count=len({r['sequence_sha256'] for r in records}),guides=sorted(guides),family_assignment_rows=len(families),
        missing_assignment_rows=sum(r['status']!='assigned' for r in families),
        families_by_guide={g:len({r['family'] for r in families if r['guide']==g and r['status']=='assigned'}) for g in guides},
        sources={str(p):sha(p) for p in [marker_path,raw,decisions_path,fasta,dbpath,dbpath.parent/'receipt.json',audit_path]},script_sha256=sha(__file__),
        artifacts={p.name:sha(p) for p in out.iterdir()},scope='Every BUSCO hit retained with annotation gene identities and both existing family assignments; no best-copy selection, new structure prediction, biological duplication classification or orthology claim.')
    write_json(out/'receipt.json',result);write_json(Path('metadata/aphelid_marker_copy_inventory_20260928.json'),result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['sources','artifacts']},indent=2))

if __name__=='__main__':main()
