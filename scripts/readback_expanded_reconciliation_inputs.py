#!/usr/bin/env python3
"""Independently check staged expanded partitions, copied trees and native paths."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
from orthofinder.tools import mcl, tree
from orthofinder.utils import files


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    receipt=json.loads((a.inputs/'receipt.json').read_text())
    summaries=[]
    for guide in receipt['guides']:
        base=(a.inputs/guide['guide']).resolve()
        wd=base/'Source/WorkingDirectory'
        manifest=base/'copied_files.tsv'
        if sha(manifest)!=guide['manifest_sha256']:
            raise ValueError('Manifest changed')
        with manifest.open() as handle:
            copied=list(csv.DictReader(handle,delimiter='\t'))
        expected_files={r['relative_path'] for r in copied}
        if len(expected_files)!=len(copied):
            raise ValueError('Duplicated copied path')
        actual_files={str(p.relative_to(base)) for p in base.rglob('*') if p.is_file()}
        if actual_files != expected_files|{'copied_files.tsv','Source/Log.pending.txt'}:
            raise ValueError('Unexpected file universe')
        for row in copied:
            dest=base/row['relative_path']
            source=Path(row['source'])
            if (dest.is_symlink() or dest.stat().st_ino==source.stat().st_ino
                    or dest.stat().st_size!=int(row['bytes'])
                    or sha(dest)!=row['sha256'] or sha(source)!=row['sha256']):
                raise ValueError('Copy integrity mismatch')
        cluster=wd/'clusters_OrthoFinder.txt_id_pairs.txt'
        groups=mcl.GetPredictedOGs(str(cluster))
        if len(groups)!=guide['families'] or sum(map(len,groups))!=5815847:
            raise ValueError('Wrong partition dimensions')
        species={line.split(':',1)[0] for line in (wd/'SpeciesIDs.txt').read_text().splitlines() if line.strip()}
        if len(species)!=526:
            raise ValueError('Wrong species universe')
        pending={r['family'] for r in guide['pending_families']}
        expected_trees={f'OG{i:07d}' for i,g in enumerate(groups) if len(g)>=3}
        observed={p.stem for p in (wd/'Trees_ids').glob('*.txt')}
        if observed!=expected_trees-pending or observed&pending:
            raise ValueError('Wrong tree universe')
        checked=tips_total=0
        for i,genes in enumerate(groups):
            if not {g.split('_',1)[0] for g in genes}<=species:
                raise ValueError('Unknown taxon')
            name=f'OG{i:07d}'
            if len(genes)<3 or name in pending:
                continue
            text=(wd/'Trees_ids'/f'{name}.txt').read_text().strip()
            if text.count(';')!=1 or not text.endswith(';'):
                raise ValueError('Incomplete tree')
            parsed=tree.Tree(text,format=0)
            tips=parsed.get_leaf_names()
            if len(tips)!=len(set(tips)) or set(tips)!=set(genes):
                raise ValueError('Native tree differs from actual expanded family: '+name)
            edges=0
            for node in parsed.traverse():
                if node.is_root():
                    continue
                edges+=1
                if not math.isfinite(node.dist) or node.dist<0:
                    raise ValueError('Invalid branch length')
            if text.count(':') not in [edges,edges+1]:
                raise ValueError('Missing explicit lengths')
            checked+=1
            tips_total+=len(tips)
        descriptor=base/'Source/Log.pending.txt'
        if sha(descriptor)!=guide['pending_descriptor_sha256'] or (base/'Source/Log.txt').exists():
            raise ValueError('Unexpected runnable or changed descriptor')
        with tempfile.TemporaryDirectory(prefix='expanded-native-locator-') as tmp:
            (Path(tmp)/'Log.txt').write_text(descriptor.read_text())
            locator=files.PreviousFilesLocator_new(SimpleNamespace(qStartFromFasta=False,qStartFromBlast=False),tmp)
            if ([Path(p) for p in locator.wd_base_prev]!=[wd]
                    or Path(locator.wd_trees)!=wd or Path(locator.clustersFilename_pairs)!=cluster):
                raise ValueError('Native paths escape isolated expanded inputs')
        summary=dict(guide=guide['guide'],families=len(groups),proteins=sum(map(len,groups)),
                     trees_checked=checked,tree_tips_checked=tips_total,pending=sorted(pending),
                     physical_files_checked=len(copied),native_locator_isolated=True)
        summaries.append(summary)
        print(json.dumps(summary),flush=True)
        del groups
    complete=all(not g['pending'] for g in summaries)
    result=dict(status='passed_complete_expanded_input_readback' if complete else 'passed_complete_available_expanded_input_readback',guides=summaries,
                input_receipt_sha256=sha(a.inputs/'receipt.json'),script_sha256=sha(Path(__file__)),
                launch_ready=False,scope='All physical copies and native tree tips checked against full expanded partitions; pending families, if any, remain explicit. No reconciliation or final species-tree validation.')
    a.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
