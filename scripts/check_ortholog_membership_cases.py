#!/usr/bin/env python3
"""Compare native stream membership to a Python set; reject corrupt inputs."""
import argparse,json,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);a=p.parse_args()
def key(x,y):return f'{x:06x}{y:06x}'
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);stream=root/'stream';queries=root/'queries';output=root/'output'
    pairs=[(1,3),(2,9),(10,11)];targets=[(0,1),(1,3),(1,9),(2,9),(5,7),(10,11),(12,13)]
    native=''.join(key(x,y)+'0\n'+key(x,y)+'1\n' for x,y in pairs)
    query=''.join(key(x,y)+'\n' for x,y in targets)
    def run(s,q):
        stream.write_text(s);queries.write_text(q)
        return subprocess.run([str(a.binary.resolve()),str(stream),str(queries),str(output)],capture_output=True,text=True)
    r=run(native,query);assert r.returncode==0,r.stderr
    assert output.read_text()==''.join(key(x,y)+'\t'+str(int((x,y) in set(pairs)))+'\n' for x,y in targets)
    assert json.loads(r.stdout)==dict(stream_pairs=3,queries=7,present=3,absent=4)
    assert run('',query).returncode==0
    assert run(native,'').returncode==0
    for s,q in [(native[:-1],query),(native.replace('1\n','0\n',1),query),(native+native[:28],query),(native,query+query[:13]),(native,query[:-1]),(native,key(9,2)+'\n')]:
        assert run(s,q).returncode!=0
print('Exact membership, leading/interior/trailing absence, empty streams/queries and six corruption cases passed.')
