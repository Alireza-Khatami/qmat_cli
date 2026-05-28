"""For the cost-divergent edges, see whether opt[] also diverges, and what
Apply(opt) looks like on both sides — this isolates whether the divergence is
the position solver or something downstream."""
import json, os, re, math

BASE = r'C:\Users\alirz\Projects\Graphics\QMAT_old working version  exe file\qmat_x64\qmat\output\01_00040057_f8f78dbd17414efda75bc437_trimesh_000'
P = '01_00040057_f8f78dbd17414efda75bc437_trimesh_000'

def load(path):
    raw = open(path, 'rb').read().decode('utf-8', errors='replace')
    raw = re.sub(r'(?<![A-Za-z_0-9])(-?)inf(?![A-Za-z_0-9])', r'\1Infinity', raw)
    raw = re.sub(r'(?<![A-Za-z_0-9])nan(?![A-Za-z_0-9])', 'NaN', raw)
    return json.loads(raw)

ours = load(os.path.join(BASE, P + '_initial_heap_state_.json'))
ml   = load(os.path.join(BASE, P + '_mat_initial.off_initial_heap_state_.json'))

o = {(e['v0_id'], e['v1_id']): e for e in ours['entries']}
m = {(e['v0_id'], e['v1_id']): e for e in ml['entries']}

# Tabulate cost-divergent edges (relerr >= 0.1)
divergent = []
for k in o:
    if k not in m: continue
    oe, me = o[k], m[k]
    oc, mc = oe['cost'], me['cost']
    rel = abs(oc - mc) / max(abs(oc), abs(mc), 1e-300)
    if rel >= 0.1:
        divergent.append((rel, k, oe, me))
divergent.sort(reverse=True)

def vdist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

print(f'\nEdges with cost relerr >= 0.1: {len(divergent)}')
print('-- ours heap has only {rank,v0,v1,cost,opt}; ML heap has breakdown.')
print('relerr  edge          ours_cost   ml_cost     opt_dist   ours_opt(x,y,z)               ml_opt(x,y,z)                 ml_applyOpt  ml_applyMid  ml_gate     ml_quadErr  ml_newQual')
for rel, k, oe, me in divergent[:30]:
    od = vdist(oe['opt'], me['opt'])
    print(f'{rel:.3f}  ({k[0]:5d},{k[1]:5d})  {oe["cost"]:.3e}  {me["cost"]:.3e}  {od:.3e}  ({oe["opt"][0]:+.5f},{oe["opt"][1]:+.5f},{oe["opt"][2]:+.5f})  ({me["opt"][0]:+.5f},{me["opt"][1]:+.5f},{me["opt"][2]:+.5f})  {me.get("applyOpt", 0):.3e}  {me.get("applyMid", 0):.3e}  {me.get("gate", 0):.3e}  {me.get("quadErr", 0):.3e}  {me.get("newQual", 0):.3e}')

# For each divergent edge, what does applyOpt show?
# If our applyOpt << ml applyOpt and opt_dist > 0, that means our solver found a
# "better" minimum (closer to true min of quadric) than vcg's — i.e., we're not
# wrong, we're just more aggressive.
print('\nSummary:')
print(f'  edges where opt[] differs > 1e-4 : {sum(1 for _,_,oe,me in divergent if vdist(oe["opt"],me["opt"]) > 1e-4)}')
print(f'  edges where opt[] differs < 1e-4 : {sum(1 for _,_,oe,me in divergent if vdist(oe["opt"],me["opt"]) <= 1e-4)}')
# How often is ours applyOpt < ml applyOpt?
# Skip: we don't have ours applyOpt
# To get it, we'd need to extend WriteInitialHeapState.

# Compute the relative size of opt vs midpoint distance
print('\nFor worst-10, how does opt compare to midpoint?')
def vmid(oe):  # use ours' applyMid to back out — we don't have midpoint coords
    return None
print('(we have applyOpt, applyMid — but not midpoint coords directly)')

# Per-vertex check: do these involve specific 'hot' vertices?
from collections import Counter
verts = Counter()
for _, (a,b), _, _ in divergent:
    verts[a] += 1
    verts[b] += 1
print(f'\nTop-10 vertices appearing in divergent edges:')
for v, c in verts.most_common(10):
    print(f'  v={v:5d}  count={c}')
