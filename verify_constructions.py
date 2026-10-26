"""Statevector checks of the constructions in the N-qubit W-state paper."""
import itertools
import numpy as np

def w(n):
    v = np.zeros(2 ** n); v[[1 << j for j in range(n)]] = 1
    return v / np.linalg.norm(v)

# --- Proposition (step W_n -> W_{n+1}) ------------------------------------
# qubit j is bit j of the index; qubit n is the new qubit.
for n in range(2, 11):
    N = n + 1
    s = np.kron(np.array([1, 1]) / np.sqrt(2), w(n))       # new qubit is bit n
    t = np.zeros_like(s)
    for i, a in enumerate(s):
        t[i ^ (((i >> n) & 1) << 0)] += a                   # CNOT new -> qubit 0
    flag = np.array([((i >> 0) & 1) & ((i >> n) & 1) for i in range(2 ** N)])
    ok = t * (flag == 0); p = ok @ ok
    rest = t * (flag == 1)
    assert abs(p - (n + 1) / (2 * n)) < 1e-12
    assert np.allclose(ok / np.sqrt(p), w(N))
    # failure branch: qubits 0 and n are 1, qubits 1..n-1 hold W_{n-1}
    exp = np.zeros(2 ** N)
    for j in range(1, n):
        exp[(1 << 0) | (1 << n) | (1 << j)] = 1
    assert np.allclose(rest / np.linalg.norm(rest), exp / np.linalg.norm(exp))
print('step proposition: OK for n = 2..10')

# --- New 1:8 split: |T>^4, parity into output, CZ on (2,3),(2,4),(3,4), X-measure
om = np.exp(1j * np.pi / 4)
q = lambda x: (x[1] & x[2]) ^ (x[1] & x[3]) ^ (x[2] & x[3])
ptot = 0
for s in itertools.product((0, 1), repeat=4):
    o = np.zeros(2, complex)
    for x in itertools.product((0, 1), repeat=4):
        o[sum(x) % 2] += om ** sum(x) * (-1) ** (q(x) + np.dot(s, x)) / 16
    pn = (np.abs(o) ** 2).sum()
    if s in [(0,) * 4, (1,) * 4]:
        o2 = o * np.array([1, 1 if s == (0,) * 4 else -1])      # Z correction
        v = o2 / np.sqrt(pn); v = v / (v[0] / abs(v[0]))         # remove global phase
        assert np.allclose(v, np.array([1, 2 * np.sqrt(2)]) / 3)
        ptot += pn
assert abs(ptot - 9 / 32) < 1e-12
print('1:8 split circuit: OK, p =', ptot)

# --- Two-AND W5 circuit (qubits 1..5 data, 4 = AND target, 6 = herald) -----
n = 7; st = np.zeros(2 ** n, complex); st[0] = 1
def bit(i, k): return (i >> k) & 1
def perm(f):
    global st
    t = np.zeros_like(st)
    for i in range(2 ** n): t[f(i)] += st[i]
    st = t
def h(k):
    global st
    t = np.zeros_like(st)
    for i in range(2 ** n):
        t[i] += st[i] * (-1) ** bit(i, k) / np.sqrt(2) + st[i ^ (1 << k)] / np.sqrt(2)
    st = t
for k in (2, 3, 5): h(k)
perm(lambda i: i ^ (1 << 1))
for c in (2, 3, 5): perm(lambda i, c=c: i ^ (bit(i, c) << 1))
perm(lambda i: i ^ ((bit(i, 2) & bit(i, 3)) << 4))
for tq in (3, 2, 1): perm(lambda i, tq=tq: i ^ (bit(i, 4) << tq))
perm(lambda i: i ^ ((bit(i, 1) & bit(i, 5)) << 6))
good = np.array([st[i] if bit(i, 6) == 0 else 0 for i in range(2 ** n)])
p = (np.abs(good) ** 2).sum()
exp = np.zeros(2 ** n); exp[[1 << k for k in range(1, 6)]] = 1
assert abs(p - 5 / 8) < 1e-12 and np.allclose(good / np.sqrt(p), exp / np.sqrt(5))
print('two-AND W5 circuit: OK, p =', p)
