"""Compare initial heaps with OptimalPlacement=true on both sides.

Both JSONs schema:
  { mesh, heap_size, order, entries: [ {rank,v0_id,v1_id,cost,opt[3],quadErr,applyOpt,applyMid,gate,newQual,minCos}, ... ] }
ML emits `inf` literal which we patch to `Infinity` before json.loads.
"""
import json, os, re, math

BASE = r'C:\Users\alirz\Projects\Graphics\QMAT_old working version  exe file\qmat_x64\qmat\output\01_00040057_f8f78dbd17414efda75bc437_trimesh_000'
P = '01_00040057_f8f78dbd17414efda75bc437_trimesh_000'

OURS = os.path.join(BASE, P + '_initial_heap_state_.json')
ML   = os.path.join(BASE, P + '_mat_initial.off_initial_heap_state_.json')

def load(path):
    raw = open(path, 'rb').read().decode('utf-8', errors='replace')
    # patch bare inf / nan tokens
    raw = re.sub(r'(?<![A-Za-z_0-9])(-?)inf(?![A-Za-z_0-9])', r'\1Infinity', raw)
    raw = re.sub(r'(?<![A-Za-z_0-9])nan(?![A-Za-z_0-9])', 'NaN', raw)
    return json.loads(raw)

ours = load(OURS)
ml   = load(ML)

print(f'ours: heap_size={ours.get("heap_size")} entries={len(ours["entries"])}')
print(f'ml  : heap_size={ml.get("heap_size")} entries={len(ml["entries"])}')

def keyed(entries):
    """Map (min,max) -> list of entries (both directions are separate entries when asymmetric)."""
    by_und = {}
    by_dir = {}
    for e in entries:
        a, b = e['v0_id'], e['v1_id']
        ukey = (min(a,b), max(a,b))
        dkey = (a, b)
        by_und.setdefault(ukey, []).append(e)
        by_dir[dkey] = e
    return by_und, by_dir

o_und, o_dir = keyed(ours['entries'])
m_und, m_dir = keyed(ml['entries'])

print(f'\nUndirected unique edges: ours={len(o_und)}, ml={len(m_und)}')
print(f'Directed entries       : ours={len(o_dir)}, ml={len(m_dir)}')

common_und = set(o_und) & set(m_und)
only_o     = set(o_und) - set(m_und)
only_m     = set(m_und) - set(o_und)
print(f'common undirected: {len(common_und)}  only-ours: {len(only_o)}  only-ml: {len(only_m)}')

# Directional comparison
common_dir = set(o_dir) & set(m_dir)
only_o_dir = set(o_dir) - set(m_dir)
only_m_dir = set(m_dir) - set(o_dir)
print(f'common directed  : {len(common_dir)}  only-ours: {len(only_o_dir)}  only-ml: {len(only_m_dir)}')

if only_o_dir:
    print(f'  sample only-ours: {list(only_o_dir)[:5]}')
if only_m_dir:
    print(f'  sample only-ml  : {list(only_m_dir)[:5]}')

# Cost agreement on common directed edges
def relerr(a, b):
    d = abs(a - b)
    s = max(abs(a), abs(b), 1e-300)
    return d / s

bins = [0, 1e-12, 1e-9, 1e-6, 1e-3, 1e-2, 0.1, 1.0, 10.0, math.inf]
labels = ['1e-12','1e-9','1e-6','1e-3','1e-2','0.1','1.0','10.0','>10']
counts = [0] * (len(bins) - 1)
big = []
for d in common_dir:
    ec = o_dir[d]['cost']; mc = m_dir[d]['cost']
    r = relerr(ec, mc)
    for i in range(len(bins) - 1):
        if r < bins[i + 1]:
            counts[i] += 1
            break
    if r > 1e-3:
        big.append((r, d, ec, mc))

print('\nRelative cost-diff histogram (common directed edges):')
for lab, c in zip(labels, counts):
    print(f'  < {lab:<6} : {c:6d}  ({100.0*c/max(1,len(common_dir)):.2f}%)')

big.sort(reverse=True)
print(f'\nWorst {min(10,len(big))} cost differences (rel >1e-3):')
for r, d, ec, mc in big[:10]:
    print(f'  ({d[0]:5d},{d[1]:5d})  ours={ec:.6e}  ml={mc:.6e}  relerr={r:.3g}')

# Compare opt[] positions on common directed edges
def vdist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

opt_diffs = []
for d in common_dir:
    if 'opt' in o_dir[d] and 'opt' in m_dir[d]:
        opt_diffs.append((vdist(o_dir[d]['opt'], m_dir[d]['opt']), d))
opt_diffs.sort(reverse=True)
print(f'\nopt[] position diffs ({len(opt_diffs)} common edges):')
if opt_diffs:
    print(f'  max  : {opt_diffs[0][0]:.6e}  edge={opt_diffs[0][1]}')
    print(f'  p99  : {opt_diffs[len(opt_diffs)//100][0]:.6e}')
    print(f'  med  : {opt_diffs[len(opt_diffs)//2][0]:.6e}')
    print(f'  min  : {opt_diffs[-1][0]:.6e}')

# Heap-order disagreement: how far apart are common edges in the sorted list?
o_rank = {(e['v0_id'], e['v1_id']): e['rank'] for e in ours['entries']}
m_rank = {(e['v0_id'], e['v1_id']): e['rank'] for e in ml['entries']}
rank_diffs = []
for d in common_dir:
    rank_diffs.append((abs(o_rank[d] - m_rank[d]), d))
rank_diffs.sort(reverse=True)
print(f'\nHeap-rank disagreement (common directed edges):')
print(f'  max abs-rank-diff: {rank_diffs[0][0]}  edge={rank_diffs[0][1]}')
print(f'  p99       : {rank_diffs[len(rank_diffs)//100][0]}')
print(f'  med       : {rank_diffs[len(rank_diffs)//2][0]}')

# Top of heap agreement
def top_n(entries, n):
    return [(e['v0_id'], e['v1_id']) for e in entries[:n]]
for N in (5, 10, 20, 50, 100):
    o_top = top_n(ours['entries'], N)
    m_top = top_n(ml['entries'], N)
    overlap = len(set(o_top) & set(m_top))
    print(f'top-{N:<3} directed overlap: {overlap}/{N}')

# Where the FIRST disagreement happens
first_disagree = None
for i in range(min(len(ours['entries']), len(ml['entries']))):
    o_e = ours['entries'][i]
    m_e = ml['entries'][i]
    if (o_e['v0_id'], o_e['v1_id']) != (m_e['v0_id'], m_e['v1_id']):
        first_disagree = i
        break
print(f'\nFirst rank where ordered edge differs: {first_disagree}')
if first_disagree is not None:
    for i in range(max(0, first_disagree-1), min(len(ours['entries']), first_disagree + 4)):
        o_e = ours['entries'][i]; m_e = ml['entries'][i]
        print(f'  rank {i:5d}  ours ({o_e["v0_id"]:5d},{o_e["v1_id"]:5d})={o_e["cost"]:.4e}  ml ({m_e["v0_id"]:5d},{m_e["v1_id"]:5d})={m_e["cost"]:.4e}')
