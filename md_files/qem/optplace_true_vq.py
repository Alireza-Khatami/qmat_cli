"""Verify vertex quadrics still bit-match with OptimalPlacement=true (they should — placement mode doesn't affect quadric build)."""
import json, os, re, math

BASE = r'C:\Users\alirz\Projects\Graphics\QMAT_old working version  exe file\qmat_x64\qmat\output\01_00040057_f8f78dbd17414efda75bc437_trimesh_000'
P = '01_00040057_f8f78dbd17414efda75bc437_trimesh_000'

def load(path):
    raw = open(path, 'rb').read().decode('utf-8', errors='replace')
    raw = re.sub(r'(?<![A-Za-z_0-9])(-?)inf(?![A-Za-z_0-9])', r'\1Infinity', raw)
    raw = re.sub(r'(?<![A-Za-z_0-9])nan(?![A-Za-z_0-9])', 'NaN', raw)
    return json.loads(raw)

ours = load(os.path.join(BASE, P + '_vertex_quadrics_.json'))
ml   = load(os.path.join(BASE, P + '_mat_initial.off_vertex_quadrics_.json'))

def normalize(blob):
    if isinstance(blob, list): return blob
    if isinstance(blob, dict):
        for k in ('entries', 'quadrics', 'vertex_quadrics', 'data', 'verts', 'vertices'):
            if k in blob and isinstance(blob[k], list): return blob[k]
    return []

oe = normalize(ours)
me = normalize(ml)
print(f'ours: {len(oe)}  ml: {len(me)}')
if oe: print(f'ours[0]: {oe[0]}')
if me: print(f'ml[0]: {me[0]}')

o_by = {e['id']: e for e in oe}
m_by = {e['id']: e for e in me}
common = set(o_by) & set(m_by)
print(f'common vertices: {len(common)}')

max_a = 0.0; max_b = 0.0; max_c = 0.0
arg_a = arg_b = arg_c = None
for v in common:
    o = o_by[v]; m = m_by[v]
    for i in range(6):
        d = abs(o['a'][i] - m['a'][i])
        if d > max_a: max_a, arg_a = d, (v, i, o['a'][i], m['a'][i])
    for i in range(3):
        d = abs(o['b'][i] - m['b'][i])
        if d > max_b: max_b, arg_b = d, (v, i, o['b'][i], m['b'][i])
    d = abs(o['c'] - m['c'])
    if d > max_c: max_c, arg_c = d, (v, o['c'], m['c'])

print(f'\nmax abs diff a[i]: {max_a:.3e}  at v={arg_a[0]} i={arg_a[1]} ours={arg_a[2]:.6e} ml={arg_a[3]:.6e}')
print(f'max abs diff b[i]: {max_b:.3e}  at v={arg_b[0]} i={arg_b[1]} ours={arg_b[2]:.6e} ml={arg_b[3]:.6e}')
print(f'max abs diff c   : {max_c:.3e}  at v={arg_c[0]} ours={arg_c[1]:.6e} ml={arg_c[2]:.6e}')
