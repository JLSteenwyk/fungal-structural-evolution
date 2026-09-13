#!/usr/bin/env python3
"""Normalize paired-site ASA and quantify reference-scale threshold sensitivity."""
import argparse,csv,gzip,json,math
from collections import Counter
from pathlib import Path
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha
from prepare_paired_phylogenetic_inputs import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['projection','snapshot','config','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable output')
    source=checked_receipt(a.projection);checked_receipt(a.snapshot);c=json.loads(a.config.read_text())
    if source['source_receipts']['snapshot']!=sha(a.snapshot/'receipt.json'):raise ValueError('Snapshot differs')
    if sha(Path(c['source_path']))!=c['source_sha256']:raise ValueError('Reference source changed')
    scales=c['scales'];names=list(scales)
    if len(names)!=2 or any(set(v)!=set('ACDEFGHIKLMNPQRSTVWY') or any(not math.isfinite(x) or x<=0 for x in v.values()) for v in scales.values()):raise ValueError('Invalid normalization scales')
    lengths={m['model_id']:m['length'] for m in json.loads((a.snapshot/'model_provenance.json').read_text())}
    thresholds=c['diagnostic_thresholds'];counts=Counter();cross=Counter();n=0;terminals=0
    a.output.mkdir(parents=True)
    with gzip.open(a.projection/'paired_site_accessibility.tsv.gz','rt') as f,gzip.open(a.output/'normalized_paired_sites.tsv.gz','wt') as out:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames+['normalization_status']+['rsa_'+x for x in names]
        writer=csv.DictWriter(out,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in reader:
            n+=1;pos=int(row['protein_residue_1based']);length=lengths[row['model_id']];aa=row['amino_acid'];area=float(row['sasa_angstrom_squared'])
            if not 1<=pos<=length or not math.isfinite(area) or area<0:raise ValueError('Invalid source values')
            if pos in (1,length):
                terminals+=1;row['normalization_status']='terminal_not_normalized';row.update({'rsa_'+x:'' for x in names})
            else:
                values={name:area/scales[name][aa] for name in names};row['normalization_status']='internal_residue_normalized'
                for name,value in values.items():
                    row['rsa_'+name]=value;counts[name,'rows']+=1;counts[name,'above_one']+=int(value>1)
                    for t in thresholds:counts[name,str(t)]+=int(value<t)
                for t in thresholds:cross[str(t)]+=int((values[names[0]]<t)!=(values[names[1]]<t))
            writer.writerow(row)
    if n!=source['observed_sites_linked']:raise ValueError('Source row count differs')
    summaries=[{'scale':name,'normalized_rows':counts[name,'rows'],'above_one':counts[name,'above_one'],**{'below_'+str(t):counts[name,str(t)] for t in thresholds}} for name in names]
    write_table(a.output/'scale_summary.tsv',summaries)
    write_table(a.output/'threshold_sensitivity.tsv',[{'threshold':t,'eligible_rows':n-terminals,'scale_disagreements':cross[str(t)],'fraction_disagreement':cross[str(t)]/(n-terminals) if n>terminals else ''} for t in thresholds])
    r={'status':'complete_reference_normalization','source_receipt_sha256':sha(a.projection/'receipt.json'),'snapshot_receipt_sha256':sha(a.snapshot/'receipt.json'),'normalization_config_sha256':sha(a.config),'script_sha256':sha(Path(__file__)),'rows':n,'terminal_rows_not_normalized':terminals,'interpretation':'Unclipped reference-normalized isolated-chain ASA. Threshold counts measure sensitivity to normalization convention, not biological core/surface truth or independent evolutionary events. Amino-acid identity, prediction confidence, domain orientation and missing partners require controls. DSSP and ShrakeRupley implementation equivalence is unestablished.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    (a.output/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
