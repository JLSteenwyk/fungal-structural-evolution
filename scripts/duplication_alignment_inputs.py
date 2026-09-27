"""C-alpha alignment inputs with explicit original-residue mappings."""
import math
from Bio.Data.PDBData import protein_letters_1to3,protein_letters_3to1


def render_ca(record,threshold=None):
    if record['status']!='validated':raise ValueError('Validated coordinate record required')
    sequence=record['sequence'];xyz=record['ca_xyz'];confidence=record['ca_plddt']
    if len(sequence)!=len(xyz) or len(sequence)!=len(confidence):raise ValueError('Coordinate dimensions differ')
    if threshold not in (None,70):raise ValueError('Unplanned confidence mask')
    selected=[];lines=[];letters=[]
    for pos,(aa,point,b) in enumerate(zip(sequence,xyz,confidence),1):
        if len(point)!=3 or not all(math.isfinite(float(v)) for v in [*point,b]) or not 0<=b<=100:raise ValueError('Invalid coordinate/confidence')
        if aa not in protein_letters_1to3:raise ValueError('Noncanonical sequence')
        if threshold is not None and b<threshold:continue
        serial=len(selected)+1;comp=protein_letters_1to3[aa].upper();x,y,z=point
        line=f'ATOM  {serial:5d}  CA  {comp:3s} A{pos:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{b:6.2f}           C  '
        if len(line)!=80:raise ValueError('PDB coordinate or residue field overflow')
        if int(line[22:26])!=pos or protein_letters_3to1[line[17:20]]!=aa:raise ValueError('Serialized identity differs')
        if any(abs(float(line[start:stop])-value)>.000501 for start,stop,value in [(30,38,x),(38,46,y),(46,54,z)]):raise ValueError('Coordinate rounding error')
        if abs(float(line[60:66])-b)>.005001:raise ValueError('Confidence rounding error')
        selected.append(pos);letters.append(aa);lines.append(line+'\n')
    return (''.join(lines)+'TER\nEND\n').encode(),''.join(letters),selected
