"""Exhaustive search for post-selected Clifford gadgets that consume t |T> states
and output a target state exactly.

Any post-selected stabilizer operation (Clifford gates, stabilizer ancillas,
Pauli measurements with post-selection) mapping t input qubits to d output
qubits is proportional to a stabilizer tensor Psi on t+d qubits. Feeding
|T>^t gives out = sum_x tau_x Psi[x, :]. Stabilizer tensors have flat singular
values, so the best success probability of a gadget built from Psi is
p = rank(Psi as 2^t x 2^d matrix) * ||out||^2  (Psi normalised).
"""
import itertools, sys
import numpy as np


def rref_subspaces(n, k):
    """Yield k x n binary RREF basis matrices, one per k-dim subspace of F2^n."""
    for pivots in itertools.combinations(range(n), k):
        free = [(r, c) for r in range(k) for c in range(n)
                if c > pivots[r] and c not in pivots]
        for bits in itertools.product((0, 1), repeat=len(free)):
            B = np.zeros((k, n), dtype=np.int64)
            for r, p in enumerate(pivots):
                B[r, p] = 1
            for (r, c), b in zip(free, bits):
                B[r, c] = b
            yield B, pivots


def search(t, d, target, tau1, tol=1e-9):
    n = t + d
    target = target / np.linalg.norm(target)
    tau = tau1
    for _ in range(t - 1):
        tau = np.kron(tau, tau1)
    if t == 0:
        tau = np.ones(1)
    weights = 1 << np.arange(n - 1, -1, -1)
    best = []
    count = 0
    for k in range(n + 1):
        ys = np.array(list(itertools.product((0, 1), repeat=k)), dtype=np.int64).reshape(2 ** k, k)
        pairs = list(itertools.combinations(range(k), 2))
        # all phase exponents (mod 4) over y: a.y + 2 * sum Q_jl y_j y_l
        A = np.array(list(itertools.product(range(4), repeat=k)), dtype=np.int64).reshape(4 ** k, k)
        lin = (A @ ys.T) % 4                                   # (4^k, 2^k)
        if pairs:
            Qs = np.array(list(itertools.product((0, 1), repeat=len(pairs))), dtype=np.int64)
            yy = np.array([[y[a] * y[b] for (a, b) in pairs] for y in ys], dtype=np.int64)
            quad = (2 * (Qs @ yy.T)) % 4                       # (2^|pairs|, 2^k)
        else:
            quad = np.zeros((1, 2 ** k), dtype=np.int64)
        for B, pivots in rref_subspaces(n, k):
            span = (ys @ B) % 2                                 # (2^k, n)
            nonpiv = [c for c in range(n) if c not in pivots]
            for cbits in itertools.product((0, 1), repeat=len(nonpiv)):
                c = np.zeros(n, dtype=np.int64)
                c[nonpiv] = cbits
                idx = ((span + c) % 2) @ weights                # basis indices
                for qrow in quad:
                    ph = (lin + qrow) % 4                       # (4^k, 2^k)
                    amps = (1j) ** ph / np.sqrt(2 ** k)
                    V = np.zeros((amps.shape[0], 2 ** n), dtype=complex)
                    V[:, idx] = amps
                    count += V.shape[0]
                    M = V.reshape(-1, 2 ** t, 2 ** d)
                    out = np.einsum('x,mxd->md', tau, M)
                    nrm = np.einsum('md,md->m', out.conj(), out).real
                    ov = np.abs(out @ target.conj()) ** 2
                    ok = np.where((nrm > tol) & (ov > (1 - 1e-9) * nrm))[0]
                    for m in ok:
                        r = np.linalg.matrix_rank(M[m], tol=1e-9)
                        best.append((r * nrm[m], r, V[m].copy() if r * nrm[m] > 0.6 else None))
    best.sort(key=lambda z: -z[0])
    return count, best


if __name__ == '__main__':
    w = np.exp(1j * np.pi / 4)
    T1 = np.array([1, w]) / np.sqrt(2)
    W = np.zeros(8, dtype=complex); W[[1, 2, 4]] = 1
    # Validation: 1-qubit target (|0> + i sqrt2 |1>)/sqrt3, the output of the
    # Paetznick-Svore 2-T RUS gate (I + i sqrt2 X)/sqrt3 on |0> (p = 3/4).
    R = np.array([1, 1j * np.sqrt(2)])
    for t in (0, 1, 2):
        cnt, best = search(t, 1, R, T1)
        print(f'validation t={t}: states={cnt} solutions={len(best)} best_p={best[0][0] if best else 0:.4f}')
    sys.stdout.flush()
    for t in (0, 1, 2):
        cnt, best = search(t, 3, W, T1)
        ps = sorted({round(b[0], 6) for b in best}, reverse=True)
        print(f'W t={t}: states={cnt} solutions={len(best)} distinct_p={ps[:8]}')
        sys.stdout.flush()
