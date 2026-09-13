#!/usr/bin/env python3
"""Audit every exported CA correspondence grid and selected raw mmCIF atom readbacks."""
import argparse,csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from audit_busco_gene_copies import sha,read_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['mapping','screen','coordinates','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new audit output')
    r=json.loads((a.mapping/'receipt.json').read_text());c=json.loads((a.mapping/'config.json').read_text())
    if r['config_sha256']!=sha(a.mapping/'config.json') or c['screen_receipt_sha256']!=sha(a.screen/'receipt.json') or c['coordinate_receipt_sha256']!=sha(a.coordinates/'receipt.json'):raise ValueError('Source lineage differs')
    sr=json.loads((a.screen/'receipt.json').read_text())
    if sha(a.screen/'sequence_correspondence.tsv')!=sr['artifacts']['sequence_correspondence.tsv']:raise ValueError('Screen changed')
    targets={x['entity_id']:x for x in read_table(a.screen/'sequence_correspondence.tsv') if x['sequence_class']=='exact_full_sequence'}
    chosen=sorted(r['entry_receipt_sha256'],key=lambda x:hashlib.sha256(x.encode()).hexdigest())[:5]
    totals=Counter();raw_rows=0;all_atoms=0
    for entry,h in r['entry_receipt_sha256'].items():
        rp=a.mapping/(entry+'.receipt.json');assert sha(rp)==h;er=json.loads(rp.read_text());table=a.mapping/(entry+'.residues.tsv.gz');assert sha(table)==er['table_sha256']
        grids=defaultdict(dict);counts=Counter();raw={}
        if entry in chosen:
            path=a.coordinates/(entry+'.cif.gz');cr=json.loads((a.coordinates/(entry+'.receipt.json')).read_text());assert sha(path)==cr['gzip_sha256']
            with gzip.open(path,'rt') as f:d=MMCIF2Dict(f)
            cols={k.removeprefix('_atom_site.'):v for k,v in d.items() if k.startswith('_atom_site.')}
            for i,name in enumerate(cols['label_atom_id']):
                if name!='CA':continue
                key=(cols['pdbx_PDB_model_num'][i],cols['label_asym_id'][i],cols['label_seq_id'][i])
                raw.setdefault(key,[]).append({k:v[i] for k,v in cols.items()})
        with gzip.open(table,'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                assert row['entry_id']==entry;key=(row['model_number'],row['label_asym_id']);pos=int(row['label_seq_id']);assert pos not in grids[key]
                grids[key][pos]=(row['target_aa'],row['entity_id']);atoms=json.loads(row['atom_records_json']);assert len(atoms)==int(row['CA_record_count']);all_atoms+=len(atoms);counts[row['CA_status']]+=1
                for atom in atoms:
                    assert atom['label_atom_id']=='CA' and atom['label_entity_id']==row['entity_id'] and atom['label_asym_id']==row['label_asym_id'] and atom['pdbx_PDB_model_num']==row['model_number'] and int(atom['label_seq_id'])==pos
                if row['CA_status']=='no_CA_observation':assert not atoms
                if row['CA_status']=='multiple_CA_records':assert len(atoms)>1
                if entry in chosen:
                    source=raw.get((row['model_number'],row['label_asym_id'],row['label_seq_id']),[]);assert len(source)==len(atoms)
                    for x,y in zip(atoms,source):assert all(y[k]==v for k,v in x.items())
                    raw_rows+=1
        for grid in grids.values():
            seq=''.join(grid[i][0] for i in range(1,len(grid)+1));entities={v[1] for v in grid.values()};assert len(entities)==1
            target=targets[entry+'_'+next(iter(entities))];assert len(seq)==int(target['target_length']) and hashlib.sha256(seq.encode()).hexdigest()==target['model_sequence_sha256']
        assert dict(counts)==er['status_counts'] and sum(counts.values())==er['position_rows'];totals.update(counts)
    assert dict(totals)==r['status_counts'] and sum(totals.values())==r['position_rows']
    a.output.mkdir(parents=True);result={'status':'complete_CA_grid_and_sampled_raw_atom_readback','mapping_receipt_sha256':sha(a.mapping/'receipt.json'),'script_sha256':sha(Path(__file__)),'entries_checked':len(r['entry_receipt_sha256']),'position_rows_checked':sum(totals.values()),'exported_CA_atoms_checked':all_atoms,'raw_mmcif_entries_checked':chosen,'raw_mmcif_position_rows_checked':raw_rows,'status_counts':dict(totals),'interpretation':'All exported grids, sequence hashes, embedded atom identities and count partitions verified. Raw mmCIF atom fields independently read back for five identity-hash-selected entries only. Geometry and experimental quality not assessed.'};(a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
