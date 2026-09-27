"""Render audited domain intervals while retaining full-protein residue numbering."""
from duplication_alignment_inputs import render_ca


def render_interval(record,start,end,threshold=None):
    if record['status']!='validated':raise ValueError('Validated source required')
    n=record['length']
    if not 1<=start<=end<=n or any(len(record[k])!=n for k in ['sequence','ca_xyz','ca_plddt']):raise ValueError('Invalid source dimensions or interval bounds')
    subset=dict(status='validated',sequence=record['sequence'][start-1:end],ca_xyz=record['ca_xyz'][start-1:end],ca_plddt=record['ca_plddt'][start-1:end])
    blob,sequence,positions=render_ca(subset,threshold);positions=[p+start-1 for p in positions]
    lines=[];index=0
    for line in blob.decode().splitlines():
        if line.startswith('ATOM  '):
            field=f'{positions[index]:4d}'
            if len(field)!=4:raise ValueError('Original residue position overflow')
            line=line[:22]+field+line[26:];index+=1
        lines.append(line+'\n')
    if index!=len(positions):raise ValueError('Serialized residue count differs')
    return ''.join(lines).encode(),sequence,positions
