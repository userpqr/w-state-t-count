"""Expected T-counts of exact W_N preparation routes, N = 3..NMAX.
Reproduces Table 2 and Fig. 2 of the N-qubit W-state paper.

Routes (all exact, all costs in T gates):
  step      : W2 (free) followed by heralded steps W_n -> W_{n+1}
              (4 T, success (n+1)/(2n), failure leaves W_{n-1}).
  hybrid    : optimal mix of CH-tree rebuilds and steps (value iteration).
  tree      : CH tree, deterministic, N = 2^k only, 2(N-2) T.
  truncate  : best exact W_M with M > N (any route), measure M - N qubits and
              keep the all-zero outcome (probability N/M); restart on failure.
  rus       : CH tree whose root split is an exact RUS split
              (1:2 -> 2 T, p=3/4; 1:4 -> 4 T, 5.26 expected [Paetznick-Svore];
               1:8 -> 4 T, p=9/32 [this work]).
"""
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
MMAX = 2 * NMAX          # largest truncation source considered
POW2 = [1 << k for k in range(12)]

def tree(n):
    return 2 * (n - 2) if n >= 2 else 0

def p_step(n):
    return (n + 1) / (2 * n)

def step_only(N):
    C = {1: 0.0}; E = 0.0
    for n in range(2, N):
        C[n] = (4 + (1 - p_step(n)) * C[n - 1]) / p_step(n)
        E += C[n]
    return E

def hybrid(N):
    V = {m: 1e6 for m in range(N + 1)}; V[N] = 0.0
    for it in range(200000):
        R = min(tree(k) + V[k] for k in POW2 if k <= N)
        new = {N: 0.0, 0: R}
        for m in range(1, N):
            st = V[2] if m == 1 else 4 + p_step(m) * V[m + 1] + (1 - p_step(m)) * V[m - 1]
            new[m] = min(st, R)
        done = it > 1000 and all(abs(new[m] - V[m]) < 1e-10 for m in V)
        V = new
        if done:
            break
    return V[0]

# root splits a:b with a,b powers of two (ratio 1:r), exact RUS cost
ROOT = {2: (2, 8 / 3), 4: (4, 5.26), 8: (4, 4 / (9 / 32))}   # r: (T per run, expected T)

def rus(N):
    m = N
    while m % 2 == 0:
        m //= 2
    r = m - 1                                   # N = (1 + r) 2^j
    if r in ROOT:
        t, e = ROOT[r]
        return t + 2 * (N - 2), e + 2 * (N - 2)
    return None

def direct(N):
    """Best expected cost of W_N without truncation, and the route used."""
    cands = {}
    if N in POW2:
        cands['tree'] = tree(N)
    else:
        ru = rus(N)
        if ru:
            cands['rus'] = ru[1]
    cands['hybrid'] = hybrid(N)
    cands['step'] = step_only(N)
    best = min(cands, key=cands.get)       # ties go to the first entry
    return cands[best], best

def best_table(nmax, mmax):
    """best[N] = (expected T, route); trunc[N] = (source M, p, expected T).
    Truncating W_M to W_K and then W_K to W_N costs the same as truncating W_M
    to W_N directly, so sources are taken from the non-truncation routes only."""
    dir_ = {M: direct(M) for M in range(3, mmax + 1)}
    best, trunc = {}, {}
    for N in range(3, nmax + 1):
        M = min(range(N + 1, mmax + 1), key=lambda M: dir_[M][0] * M / N)
        trunc[N] = (M, N / M, dir_[M][0] * M / N)
        e, route = dir_[N]
        if trunc[N][2] < e - 1e-9:
            e, route = trunc[N][2], f'truncate W{M}'
        best[N] = (e, route)
    return best, trunc

if __name__ == '__main__':
    best, trunc = best_table(NMAX, MMAX)
    print('N lb step_run step_p step_E hybrid_E tree trunc_src trunc_p trunc_E rus_run rus_E best')
    for N in range(3, NMAX + 1):
        so, hy, ru, tr = step_only(N), hybrid(N), rus(N), trunc[N]
        det = N in POW2
        print(N, max(4, N - 1), 4 * (N - 2), f'{N}/{2**(N-1)}', f'{so:.1f}', f'{hy:.1f}',
              tree(N) if det else '-', f'W{tr[0]}', f'{N}/{tr[0]}', f'{tr[2]:.1f}',
              ru[0] if ru and not det else '-', f'{ru[1]:.1f}' if ru and not det else '-',
              f'{best[N][0]:.1f}({best[N][1]})')
