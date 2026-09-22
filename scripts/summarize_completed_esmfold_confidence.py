#!/usr/bin/env python3
"""Describe confidence-filter retention; not prediction accuracy or evolution."""
import argparse
import csv
import json
from pathlib import Path
from readback_whole_proteome_family_coverage import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--encodings',type=Path,required=True)
    ap.add_argument('--output-prefix',type=Path,required=True)
    args=ap.parse_args()
    table=Path(str(args.output_prefix)+'.tsv');receipt=Path(str(args.output_prefix)+'_receipt.json')
    if table.exists() or receipt.exists():raise FileExistsError(args.output_prefix)
    source=args.encodings/'receipt.json';r=json.loads(source.read_text())
    assert r['status']=='complete_union_of_previously_qualified_encodings'
    summary=args.encodings/'model_summary.tsv'
    assert sha(summary)==r['artifacts']['model_summary.tsv']
    keys=list(r['totals']);cohorts={};paths={};seen=set()
    for c in r['source_cohorts']:
        p=Path(c['encodings'])/'receipt.json'
        assert sha(p)==c['encoding_receipt_sha256']
        cohorts[c['name']]={'cohort':c['name'],'models':0,**{k:0 for k in keys},
                           'models_without_joint_qualified_states':0,'models_below_half_joint_retention':0}
        paths[str(Path(c['encodings']))]=c['name']
    with summary.open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            assert row['model_name'] not in seen;seen.add(row['model_name'])
            c=cohorts[paths[str(Path(row['encoding_path']).parent)]]
            c['models']+=1
            v={key:int(row[key]) for key in keys}
            assert v['length']==v['valid_states']+v['invalid_states']
            assert 0<=v['valid_feature_plddt70_pae10']<=v['valid_feature_plddt70']<=v['valid_focal_plddt70']<=v['valid_states']<=v['length']
            for key,value in v.items():c[key]+=value
            c['models_without_joint_qualified_states']+=v['valid_feature_plddt70_pae10']==0
            c['models_below_half_joint_retention']+=2*v['valid_feature_plddt70_pae10']<v['length']
    assert len(seen)==r['models']
    for key in keys:assert sum(c[key] for c in cohorts.values())==r['totals'][key]
    for c in r['source_cohorts']:
        values=cohorts[c['name']]
        assert values['models']==c['models']
        original=json.loads((Path(c['encodings'])/'receipt.json').read_text())
        for key in keys:assert values[key]==original['totals'][key]
    for values in cohorts.values():
        values['joint_qualified_fraction_of_residues']=values['valid_feature_plddt70_pae10']/values['length']
    with table.open('x') as h:
        w=csv.DictWriter(h,fieldnames=list(next(iter(cohorts.values()))),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(cohorts.values())
    result={'status':'complete_descriptive_cohort_confidence_retention','models':len(seen),
            'source_receipt':str(source),'source_receipt_sha256':sha(source),
            'source_summary_sha256':sha(summary),'script_sha256':sha(__file__),
            'artifacts':{str(table):sha(table)},
            'scope':'Full protein residue counts, before marker alignment eligibility. Joint filter requires native validity, all six context residues pLDDT>=70 and maximum directional context PAE<=10. Cohorts differ in length and sampling; differences are not predictor accuracy, biological divergence or independent evolutionary observations.'}
    with receipt.open('x') as h:json.dump(result,h,indent=2);h.write('\n')
    print(table.read_text())


if __name__=='__main__':main()
