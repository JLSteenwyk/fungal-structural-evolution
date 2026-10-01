#!/usr/bin/env python3
"""Exercise variable-tip complete readers with rehashed false native exports.

Synthetic likelihoods, profiles and bootstrap trees are software fixtures only.
No model inference, biological validation, provenance closure or pilot claim.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
from readback_native_taxon_pmsf_collection import inspect_run, SUPPORT_FIELDS
from run_ortholog_pair_guide_comparison import sha


def write(path, record):
    Path(path).write_text(json.dumps(record, indent=2)+'\n')


def make_case(root, n):
    root.mkdir(parents=True)
    run, audit = root/'run', root/'run'/'audit'
    audit.mkdir(parents=True)
    taxa = ['T_'+str(i) for i in range(n)]
    matrix = root/'matrix.faa'
    matrix.write_text(''.join('>'+t+'\n'+'A'*17+'\n' for t in taxa))
    sides = [tuple(taxa[i:]) for i in range(2, n-1)]
    def tree(label):
        part = f'({taxa[-2]}:0.1,{taxa[-1]}:0.1){label}:0.1'
        for i in range(n-3, 1, -1):
            part = f'({taxa[i]}:0.1,{part}){label}:0.1'
        return f'({taxa[0]}:0.1,{taxa[1]}:0.1,{part});'
    ml, consensus, guide = tree("'80/100'"), tree('100'), tree('')
    for name, value in [('pmsf.treefile', ml), ('pmsf.contree', consensus), ('input_guide.treefile', guide)]:
        (run/name).write_text(value+'\n')
    (run/'pmsf.ufboot').write_text((guide+'\n')*1000)
    (run/'pmsf.sitefreq').write_text(''.join(str(i)+' '+' '.join(['0.05']*20)+'\n' for i in range(1,18)))
    (run/'pmsf.log').write_text('Computing posterior mean site frequencies\n')
    (run/'pmsf.iqtree').write_text(f'Input data: {n} sequences with 17 amino-acid sites\n'
        'Model of substitution: LG+SSF+F+G4\nSH-aLRT support (%) / ultrafast bootstrap support (%)\n'
        f'Log-likelihood of the tree: -100.00\nLog-likelihood of consensus tree: -100.00\n'
        f'Total tree length (sum of branch lengths): {(2*n-3)/10:.4f}\n'
        f'Robinson-Foulds distance between ML tree and consensus tree: 0\n{ml}\n{consensus}\n')
    def canon(side):
        opposite = tuple(t for t in taxa if t not in side)
        return min(tuple(sorted(side)), tuple(sorted(opposite)), key=lambda x:(len(x),x))
    splits = [canon([t]) for t in taxa]+[canon(side) for side in sides]
    (run/'pmsf.splits.nex').write_text('#NEXUS\n'+''.join(f"[{i}] '{t}'\n" for i,t in enumerate(taxa,1))+
        'MATRIX\n'+''.join('100.0000 '+' '.join(str(taxa.index(t)+1) for t in split)+',\n' for split in splits)+';\n')
    rows = [dict(tree=kind, split_taxa_json=json.dumps(canon(side)), branch_length=0.1,
        sh_alrt_percent=80.0 if kind=='ml' else None, reported_ufboot_percent=100.0,
        empirical_ufboot_percent=100.0) for kind in ['ml','consensus'] for side in sides]
    with (audit/'branch_support.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,SUPPORT_FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    command=['synthetic-IQTREE-NOT-EXECUTED','-s',str(matrix.resolve()),'-m','LG+C20+F+G4',
        '--tree-freq',str((run/'input_guide.treefile').resolve()),'--alrt','1000','-B','1000','--bnni','--boot-trees']
    write(run/'config.json',dict(command=command,pinned_files={str(matrix):sha(matrix)}))
    write(run/'receipt.json',dict(status='complete_pmsf_execution_pending_full_audit',returncode=0,taxa=n,columns=17,artifacts={}))
    write(audit/'receipt.json',dict(status='synthetic_first_audit_not_production',taxa=n,sites=17,
        bootstrap_trees=1000,empirical_splits=2*n-3,internal_support_rows=2*(n-3),ml_consensus_rf=0,
        reported_ml_log_likelihood=-100.0,reported_consensus_log_likelihood=-100.0,artifacts={}))
    rehash(root)


def rehash(root):
    run=root/'run';audit=run/'audit';matrix=root/'matrix.faa'
    c=json.loads((run/'config.json').read_text());c['command'][c['command'].index('-s')+1]=str(matrix.resolve())
    c['command'][c['command'].index('--tree-freq')+1]=str((run/'input_guide.treefile').resolve())
    c['pinned_files']={str(matrix):sha(matrix)};write(run/'config.json',c)
    r=json.loads((run/'receipt.json').read_text());r['config_sha256']=sha(run/'config.json')
    r['artifacts']={p.name:sha(p) for p in run.iterdir() if p.is_file() and p.name!='receipt.json'};write(run/'receipt.json',r)
    a=json.loads((audit/'receipt.json').read_text());a['source_receipt_sha256']=sha(run/'receipt.json')
    a['artifacts']={'branch_support.tsv':sha(audit/'branch_support.tsv')};write(audit/'receipt.json',a)


def check(root,n):
    return inspect_run(root/'run',root/'matrix.faa',n,17,root/'run'/'audit',{},'synthetic_first_audit_not_production')


def change(path, old, new):
    text=path.read_text();assert old in text;path.write_text(text.replace(old,new,1))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    summaries=[]
    for n in [5,6,8,11]:
        root=args.output/('valid-'+str(n));make_case(root,n);summary,_,_=check(root,n);summaries.append(summary)
    base=args.output/'valid-6'
    cases={
        'missing_bootstrap':lambda r:(r/'run'/'pmsf.ufboot').write_text('\n'.join((r/'run'/'pmsf.ufboot').read_text().splitlines()[:-1])+'\n'),
        'wrong_bootstrap_tip':lambda r:change(r/'run'/'pmsf.ufboot','T_5','T_wrong'),
        'negative_bootstrap_edge':lambda r:change(r/'run'/'pmsf.ufboot',':0.1',':-0.1'),
        'wrong_report_tree':lambda r:change(r/'run'/'pmsf.iqtree',':0.1',':0.2'),
        'profile_normalization':lambda r:change(r/'run'/'pmsf.sitefreq','0.05','0.5'),
        'profile_negative':lambda r:change(r/'run'/'pmsf.sitefreq','0.05','-0.05'),
        'profile_nan':lambda r:change(r/'run'/'pmsf.sitefreq','0.05','NaN'),
        'profile_index':lambda r:change(r/'run'/'pmsf.sitefreq','2 ','1 '),
        'profile_missing_site':lambda r:(r/'run'/'pmsf.sitefreq').write_text('\n'.join((r/'run'/'pmsf.sitefreq').read_text().splitlines()[:-1])+'\n'),
        'nexus_weight':lambda r:change(r/'run'/'pmsf.splits.nex','100.0000','95.0000'),
        'nexus_duplicate_label':lambda r:change(r/'run'/'pmsf.splits.nex',"[6] 'T_5'","[6] 'T_4'"),
        'nexus_missing_split':lambda r:change(r/'run'/'pmsf.splits.nex','100.0000 1,\n',''),
        'nexus_duplicate_split':lambda r:change(r/'run'/'pmsf.splits.nex','100.0000 1,\n','100.0000 1,\n100.0000 1,\n'),
        'wrong_saved_support':lambda r:change(r/'run'/'audit'/'branch_support.tsv','80.0','79.0'),
        'missing_saved_support':lambda r:(r/'run'/'audit'/'branch_support.tsv').write_text('\n'.join((r/'run'/'audit'/'branch_support.tsv').read_text().splitlines()[:-1])+'\n'),
        'duplicate_saved_support':lambda r:(r/'run'/'audit'/'branch_support.tsv').write_text((r/'run'/'audit'/'branch_support.tsv').read_text()+(r/'run'/'audit'/'branch_support.tsv').read_text().splitlines()[1]+'\n'),
        'wrong_report_rf':lambda r:change(r/'run'/'pmsf.iqtree','consensus tree: 0','consensus tree: 2'),
        'wrong_report_total_length':lambda r:change(r/'run'/'pmsf.iqtree','0.9000','1.0000'),
        'wrong_matrix_length':lambda r:change(r/'matrix.faa','A'*17,'A'*16),
        'wrong_model':lambda r:change(r/'run'/'config.json','LG+C20+F+G4','LG+F+G4'),
    }
    rejected=[]
    for name,mutate in cases.items():
        root=args.output/name;shutil.copytree(base,root);mutate(root);rehash(root)
        try:check(root,6)
        except (ValueError,AssertionError,KeyError,IndexError,TypeError) as e:rejected.append(dict(case=name,reason=str(e)))
        else:raise AssertionError('False export accepted: '+name)
    proof=dict(status='passed_variable_tip_native_reader_software_cases',valid_full_cases=summaries,
        rehashed_false_exports_rejected=rejected,script_sha256=sha('scripts/readback_native_taxon_pmsf_collection.py'),
        fixture_script_sha256=sha(Path(__file__)),scope=__doc__)
    write(args.output/'receipt.json',proof)
    print(json.dumps(dict(status=proof['status'],valid_cases=len(summaries),false_exports_rejected=len(rejected))),flush=True)


if __name__=='__main__':main()
