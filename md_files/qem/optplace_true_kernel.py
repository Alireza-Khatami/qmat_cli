"""For each cost-divergent edge:
  - sum the bit-identical vertex quadrics
  - evaluate Apply at ours_opt and at ml_opt
  - compare both to applyMid (= the midpoint Apply)
Decides whether our solver lands on a different kernel point (Apply(our_opt) ≈ Apply(ml_opt))
or a genuinely better point (Apply(our_opt) << Apply(ml_opt)).
"""
import json, os, re, math

BASE = r'C:\Users\alirz\Projects\Graphics\QMAT_old working version  exe file\qmat_x64\qmat\output\01_00040057_f8f78dbd17414efda75bc437_trimesh_000'
P = '01_00040057_f8f78dbd17414efda75bc437_trimesh_000'

def load(path):
    raw = open(path, 'rb').read().decode('utf-8', errors='replace')
    raw = re.sub(r'(?<![A-Za-z_0-9])(-?)inf(?![A-Za-z_0-9])', r'\1Infinity', raw)
    raw = re.sub(r'(?<![A-Za-z_0-9])nan(?![A-Za-z_0-9])', 'NaN', raw)
    return json.loads(raw)

ours_h = load(os.path.join(BASE, P + '_initial_heap_state_.json'))
ml_h   = load(os.path.join(BASE, P + '_mat_initial.off_initial_heap_state_.json'))
vq = load(os.path.join(BASE, P + '_vertex_quadrics_.json'))['verts']

qbyid = {q['id']: q for q in vq}

def add_q(qa, qb):
    return {
        'a': [qa['a'][i] + qb['a'][i] for i in range(6)],
        'b': [qa['b'][i] + qb['b'][i] for i in range(3)],
        'c':  qa['c'] + qb['c'],
    }

def apply_q(q, p):
    a = q['a']; b = q['b']; c = q['c']
    # Quadric::Apply: p[0]^2*a[0] + 2*p[0]*p[1]*a[1] + 2*p[0]*p[2]*a[2] + p[0]*b[0]
    #              +  p[1]^2*a[3] + 2*p[1]*p[2]*a[4] + p[1]*b[1]
    #              +  p[2]^2*a[5] + p[2]*b[2] + c
    return (p[0]*p[0]*a[0] + 2*p[0]*p[1]*a[1] + 2*p[0]*p[2]*a[2] + p[0]*b[0]
          + p[1]*p[1]*a[3] + 2*p[1]*p[2]*a[4] + p[1]*b[1]
          + p[2]*p[2]*a[5] + p[2]*b[2] + c)

ours = {(e['v0_id'], e['v1_id']): e for e in ours_h['entries']}
ml   = {(e['v0_id'], e['v1_id']): e for e in ml_h['entries']}

def vdist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

divergent = []
for k in ours:
    if k not in ml: continue
    rel = abs(ours[k]['cost'] - ml[k]['cost']) / max(abs(ours[k]['cost']), abs(ml[k]['cost']), 1e-300)
    if rel >= 0.1:
        divergent.append((rel, k))
divergent.sort(reverse=True)

print(f'{len(divergent)} edges; comparing solver outputs.')
print()
print('relerr edge          dist(our,mid)  dist(ml,mid)  dist(our,ml)  Apply(our)   Apply(ml)    ratio')

# Also: is ml_opt exactly equal to v0 or v1?
def vdist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

for rel, k in divergent[:40]:
    v0, v1 = k
    q0 = qbyid[v0]; q1 = qbyid[v1]
    p0 = q0['pos']; p1 = q1['pos']
    mid = [(p0[i] + p1[i]) * 0.5 for i in range(3)]
    qsum = add_q(q0, q1)
    a_our = apply_q(qsum, ours[k]['opt'])
    a_ml  = apply_q(qsum, ml[k]['opt'])
    do_mid = vdist(ours[k]['opt'], mid)
    dm_mid = vdist(ml[k]['opt'], mid)
    do_ml  = vdist(ours[k]['opt'], ml[k]['opt'])
    ratio = (a_our / a_ml) if a_ml != 0 else float('inf')
    print(f'{rel:.3f}  ({v0:5d},{v1:5d})  {do_mid:.3e}  {dm_mid:.3e}  {do_ml:.3e}  {a_our:.3e}  {a_ml:.3e}  {ratio:.3e}')

# Summary
mids = [(vdist(ml[k]['opt'], [(qbyid[k[0]]['pos'][i]+qbyid[k[1]]['pos'][i])*0.5 for i in range(3)])) for _, k in divergent]
ours_d = [(vdist(ours[k]['opt'], [(qbyid[k[0]]['pos'][i]+qbyid[k[1]]['pos'][i])*0.5 for i in range(3)])) for _, k in divergent]
print(f'\nML opt -> midpoint distance:  median={sorted(mids)[len(mids)//2]:.3e}  max={max(mids):.3e}  min={min(mids):.3e}')
print(f'OUR opt -> midpoint distance: median={sorted(ours_d)[len(ours_d)//2]:.3e}  max={max(ours_d):.3e}  min={min(ours_d):.3e}')
