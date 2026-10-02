#!/usr/bin/env python3
"""Compute native gCF for all eight supported baseline views and both marker sets.

All125 marker trees per alignment are frozen audited inputs. Reference internal
labels are removed so native rerooting cannot move existing bootstrap/SH labels;
original support remains available through canonical split joins. This does not
reinfer trees, assign roots or explain biological causes of discordance.
"""
import argparse
import io
import json
import os
from pathlib import Path
import subprocess
import time
from Bio import Phylo
import psutil
from species_gcf_sources import load
from audit_species_guide import edges
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def write(path,record):
    with Path(path).open('x') as handle:handle.write(json.dumps(record,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    plan=json.loads(args.plan.read_text());sources,universe,markers,bindings=load(plan,args.plan)
    out=Path(plan['output']);out.mkdir(parents=True,exist_ok=False);write(out/'config.json',dict(plan_sha256=sha(args.plan),source_hashes=bindings))
    prepared=out/'inputs';prepared.mkdir();reference_paths={};gene_paths={}
    for alignment,items in markers.items():
        path=prepared/(alignment+'.genes.tree')
        with path.open('x') as handle:
            for item in items:
                text=Path(item['path']).read_text().strip();assert text.endswith(';') and text.count(';')==1
                handle.write(text+'\n')
        gene_paths[alignment]=path;bind(bindings,path)
        write(prepared/(alignment+'.gene_order.json'),items);bind(bindings,prepared/(alignment+'.gene_order.json'))
    for label,source in sources.items():
        for kind,file in [('ml','pmsf.treefile'),('consensus','pmsf.contree')]:
            original=Path(source['spec']['run'])/file;tree=Phylo.read(original,'newick');before=edges(tree,universe)
            for node in tree.get_nonterminals():node.name=None;node.confidence=None
            stream=io.StringIO();Phylo.write(tree,stream,'newick',format_branch_length='%.10f')
            text=stream.getvalue();assert edges(Phylo.read(io.StringIO(text),'newick'),universe)==before
            path=prepared/(label+'.'+kind+'.reference.tree');path.write_text(text);bind(bindings,path);reference_paths[label,kind]=path
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');completed=[]
    for label,source in sources.items():
        for kind in ['ml','consensus']:
            for alignment in ['profile','mafft']:
                name=label+'-'+kind+'-'+alignment;folder=out/name;folder.mkdir()
                command=[plan['executable'],'-t',str(reference_paths[label,kind]),'--gcf',str(gene_paths[alignment]),
                    '--cf-verbose','-T',str(plan['resources']['cpu']),'--mem','8G','--seed',str(plan['seed']),'--prefix',str(folder/'gcf')]
                with (folder/'stdout.log').open('x') as handle:
                    start=time.monotonic();process=subprocess.Popen(command,stdout=handle,stderr=subprocess.STDOUT,env=env);p=psutil.Process(process.pid)
                    write(folder/'native_launch.json',dict(pid=p.pid,created=p.create_time(),cmdline=p.cmdline(),command=command,plan_sha256=sha(args.plan)))
                    print('native_gene_concordance_started',name,p.pid,command,flush=True);code=process.wait()
                assert code==0,(name,code)
                names=['gcf.cf.stat','gcf.cf.stat_tree','gcf.cf.branch','gcf.cf.tree','gcf.log'];assert all((folder/f).is_file() for f in names)
                record=dict(label=name,species_run=label,tree_type=kind,marker_alignment=alignment,reference=str(reference_paths[label,kind]),
                    genes=str(gene_paths[alignment]),gene_order=str(prepared/(alignment+'.gene_order.json')),native_command=command,
                    marker_count=125,taxa=526,expected_internal_branches=523,elapsed_seconds=time.monotonic()-start,returncode=0,
                    artifacts={p.name:sha(p) for p in folder.iterdir() if p.is_file()})
                write(folder/'receipt.json',dict(status='complete_native_gene_concordance_pending_full_independent_readback',**record))
                bind(bindings,folder/'receipt.json')
                for file,digest in record['artifacts'].items():bind(bindings,folder/file,digest)
                completed.append(record);print('native_gene_concordance_completed',name,len(completed),'/16',flush=True)
    assert len(completed)==16;verify(bindings)
    write(out/'receipt.json',dict(status='complete_full_native_species_gene_concordance_batch_pending_independent_readback',plan_sha256=sha(args.plan),
        references=8,marker_alignments=2,markers_per_alignment=125,native_runs=16,branch_marker_cells=1046000,runs=completed,
        source_hashes=bindings,scientific_eligibility=False,scope=plan['scope']))


if __name__=='__main__':main()
