#!/usr/bin/env python3
"""Check binary-search membership against all ordered pairs in a small universe."""
from io import BytesIO
from readback_background_orthology import membership
pairs={(1,3),(2,9),(10,11)}
def key(x,y):return f'{x:06x}{y:06x}'.encode()
stream=BytesIO(b''.join(key(x,y)+b'0\n'+key(x,y)+b'1\n' for x,y in sorted(pairs)))
for x in range(14):
    for y in range(x+1,14):assert membership(stream,len(pairs),key(x,y))==int((x,y) in pairs)
assert membership(BytesIO(),0,key(0,1))==0
try:membership(BytesIO(b'x'*28),1,key(0,1))
except AssertionError:pass
else:raise AssertionError('Corrupt reciprocal record accepted')
print('All 91 ordered queries, empty stream and corrupt-record rejection passed.')
