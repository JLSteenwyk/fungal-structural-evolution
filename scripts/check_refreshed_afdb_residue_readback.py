#!/usr/bin/env python3
"""Exercise the full AFDB residue auditor on hand-specified synthetic mappings."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

AUDITOR=Path(__file__).with_name('readback_refreshed_afdb_residues.py').resolve()


def main():
    outcomes=[]
    for corruption in [None,'position','confidence','missing','extra','duplicate_ca']:
        with tempfile.TemporaryDirectory(prefix='afdb-residue-readback-') as temp:
            root=Path(temp)
            def write(path,text):
                p=root/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text);return p
            def js(path,value):return write(path,json.dumps(value))
            def sha(path):return hashlib.sha256((root/path).read_bytes()).hexdigest()
            def table(path,rows):
                p=root/path;p.parent.mkdir(parents=True,exist_ok=True)
                with p.open('w') as f:
                    w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
            mapping='results/structural_markers/mapped'
            catalog='results/structural_markers/catalog'
            markers=[];models=[];links=[];selected=[];residues=[]
            seqs=['ACDE','ACDF']
            hand_positions=[[(1,1),(2,2),(3,3),(4,4)],[(1,1),(3,3),(4,4)]]
            for i,seq in enumerate(seqs):
                taxon='F'+str(i);ident='AF-X'+str(i);path='coordinates/'+ident+'.cif'
                confidence=[80,60,90,70]
                atoms=''.join(f'CA {j+1} {v}\n' for j,v in enumerate(confidence))
                if corruption=='duplicate_ca' and i==0:atoms+='CA 1 80\n'
                write(path,'data_test\n_entity_poly.pdbx_seq_one_letter_code_can '+seq+'\nloop_\n_atom_site.label_atom_id\n_atom_site.label_seq_id\n_atom_site.B_iso_or_equiv\n'+atoms+'#\n')
                sequence_sha=hashlib.sha256(seq.encode()).hexdigest()
                markers.append(dict(marker='M',taxon_id=taxon,protein_id='p'+str(i),sequence_sha256=sequence_sha))
                models.append(dict(path=path,sha256=sha(path),model_id=ident,version=1))
                selected.append(dict(markers[-1],model_id=ident,version=1,model_path=path))
                values=[confidence[r-1] for _,r in hand_positions[i]]
                links.append(dict(markers[-1],model_id=ident,model_version=1,model_path=path,
                    retained_marker_residues=len(values),mean_retained_ca_plddt=sum(values)/len(values),
                    fraction_retained_ca_plddt_ge70=sum(x>=70 for x in values)/len(values)))
                for column,pos in hand_positions[i]:
                    residues.append(dict(marker='M',taxon_id=taxon,protein_id='p'+str(i),model_id=ident,model_version=1,
                        matrix_column_1based=column,protein_residue_1based=pos,ca_plddt=confidence[pos-1]))
            if corruption=='position':residues[0]['protein_residue_1based']=2
            if corruption=='confidence':residues[0]['ca_plddt']=79
            if corruption=='missing':residues.pop()
            if corruption=='extra':residues.append(dict(residues[-1]))
            table('results/phylogeny/markers-full-v1/protein_mapping.tsv',markers)
            write('results/phylogeny/markers-full-v1/unaligned/M.faa','>F0\nACDE\n>F1\nACDF\n')
            write('results/phylogeny/profile-alignments-full-v1/M.sto','# STOCKHOLM 1.0\nF0 AC-DE\nF1 A-CDF\n//\n')
            js('results/phylogeny/profile-alignments-full-v1/M.receipt.json',dict(retained_stockholm_columns_1based=[1,2,4,5]))
            table('results/phylogeny/profile-matrix-50-v1/site_mapping.tsv',[dict(marker='M',alignment_column_1based=i,matrix_column_1based=i) for i in range(1,5)])
            table(catalog+'/marker_model_links.tsv',selected)
            table(mapping+'/marker_structure_links.tsv',links)
            js(mapping+'/model_provenance.json',models)
            path=root/mapping/'matrix_to_structure_residues.tsv.gz'
            with gzip.open(path,'wt') as f:
                w=csv.DictWriter(f,list(residues[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(residues)
            artifact_names=['marker_structure_links.tsv','model_provenance.json','matrix_to_structure_residues.tsv.gz']
            js(mapping+'/receipt.json',dict(artifacts={x:sha(mapping+'/'+x) for x in artifact_names},matrix_residue_links=7,marker_proteins_linked=2))
            pinned={str(p.relative_to(root)):sha(str(p.relative_to(root))) for p in root.rglob('*') if p.is_file()}
            js('preparation.json',dict(pins=pinned,controller_output='controller',mapping_output=mapping,catalog=catalog))
            js('controller/mapping/receipt.json',dict(status='complete_refreshed_mapping_and_catalog_agreement_pending_residue_readback',result_receipt_sha256=sha(mapping+'/receipt.json'),plan_sha256=sha('preparation.json')))
            js('plan.json',dict(preparation_plan='preparation.json',output='audit',pins={'preparation.json':sha('preparation.json')},predecessor_pid=0,predecessor_start_ticks='0',resources=dict(minimum_free_disk_gib=0)))
            process=subprocess.run([sys.executable,str(AUDITOR),'--plan','plan.json'],cwd=root,capture_output=True,text=True,timeout=30)
            if corruption is None:
                assert process.returncode==0,process.stderr
                result=json.loads((root/'audit/receipt.json').read_text())
                assert result['matrix_residue_links']==7 and result['models']==2 and result['marker_links']==2
            else:
                assert process.returncode!=0,corruption+' accepted'
                assert not (root/'audit/receipt.json').exists()
                assert json.loads((root/'audit/state.json').read_text())['status']=='failed_requires_review'
            outcomes.append(dict(case=corruption or 'known_valid_projection',passed=True))
    print(json.dumps(dict(status='passed_full_residue_auditor_synthetic_checks',cases=outcomes),indent=2))


if __name__=='__main__':main()
