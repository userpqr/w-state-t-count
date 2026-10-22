import sys, pickle, numpy as np
sys.path.insert(0, sys.argv[1])
from gadget_search import search
w = np.exp(1j*np.pi/4); T1 = np.array([1, w])/np.sqrt(2)
W = np.zeros(8, dtype=complex); W[[1,2,4]] = 1
cnt, best = search(3, 3, W, T1)
ps = {}
for p, r, v in best: ps[round(p, 6)] = ps.get(round(p, 6), 0) + 1
print('states', cnt, 'solutions', len(best))
print('p distribution', sorted(ps.items(), reverse=True)[:15])
top = [b for b in best if b[2] is not None][:20]
pickle.dump(top, open(sys.argv[1] + '/t3_top.pkl', 'wb'))
