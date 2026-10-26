"""Multi-target version of gadget_search.search: best success probability of a
post-selected Clifford gadget consuming t |T> states whose d-qubit output is
proportional to each target."""
import itertools, sys, time, json
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from gadget_search import rref_subspaces

def search_multi(t, d, targets, tol=1e-9):
    n=t+d
    w8=np.exp(1j*np.pi/4); tau1=np.array([1,w8])/np.sqrt(2)
    tau=np.ones(1)
    for _ in range(t): tau=np.kron(tau,tau1)
    T=np.array([v/np.linalg.norm(v) for v in targets.values()]).conj()   # (K, 2^d)
    names=list(targets); best={k:0.0 for k in names}; nsol={k:0 for k in names}
    weights=1<<np.arange(n-1,-1,-1)
    for k in range(n+1):
        ys=np.array(list(itertools.product((0,1),repeat=k)),dtype=np.int64).reshape(2**k,k)
        pairs=list(itertools.combinations(range(k),2))
        A=np.array(list(itertools.product(range(4),repeat=k)),dtype=np.int64).reshape(4**k,k)
        lin=(A@ys.T)%4
        if pairs:
            Qs=np.array(list(itertools.product((0,1),repeat=len(pairs))),dtype=np.int64)
            yy=np.array([[y[a]*y[b] for (a,b) in pairs] for y in ys],dtype=np.int64)
            quad=(2*(Qs@yy.T))%4
        else: quad=np.zeros((1,2**k),dtype=np.int64)
        for B,piv in rref_subspaces(n,k):
            span=(ys@B)%2; nonpiv=[c for c in range(n) if c not in piv]
            for cb in itertools.product((0,1),repeat=len(nonpiv)):
                c=np.zeros(n,dtype=np.int64); c[nonpiv]=cb
                idx=((span+c)%2)@weights
                for q in quad:
                    amps=(1j)**((lin+q)%4)/np.sqrt(2**k)
                    V=np.zeros((amps.shape[0],2**n),complex); V[:,idx]=amps
                    M=V.reshape(-1,2**t,2**d)
                    out=np.einsum('x,mxd->md',tau,M)
                    nrm=np.einsum('md,md->m',out.conj(),out).real
                    ov=np.abs(out@T.T)**2          # (m,K)
                    good=(nrm[:,None]>tol)&(ov>(1-1e-9)*nrm[:,None])
                    for mi,ki in zip(*np.where(good)):
                        r=np.linalg.matrix_rank(M[mi],tol=1e-9); p=r*nrm[mi]
                        nm=names[ki]; nsol[nm]+=1
                        if p>best[nm]+1e-12: best[nm]=p
    return best,nsol

if __name__=='__main__':
    mode=sys.argv[1]; tmax=int(sys.argv[2])
    if mode=='split':
        d=1
        targets={f'1:{m}':np.array([1,np.sqrt(m)]) for m in (2,3,4,8,16)}
        targets.update({'3:4':np.array([np.sqrt(3),2]),'2:3':np.array([np.sqrt(2),np.sqrt(3)]),'3:5':np.array([np.sqrt(3),np.sqrt(5)])})
    else:
        d=3
        def ww(a,b,c):
            v=np.zeros(8); v[4]=np.sqrt(a); v[2]=np.sqrt(b); v[1]=np.sqrt(c); return v
        targets={'1:1:1':ww(1,1,1),'1:2:4':ww(1,2,4),'1:2:8':ww(1,2,8),'1:4:8':ww(1,4,8),'1:2:16':ww(1,2,16),'1:4:16':ww(1,4,16),'1:1:2':ww(1,1,2),'1:1:4':ww(1,1,4)}
    for t in range(0,tmax+1):
        t0=time.time(); best,ns=search_multi(t,d,targets)
        print(json.dumps({'t':t,'best':{k:round(v,6) for k,v in best.items()},'nsol':ns,'sec':round(time.time()-t0)}),flush=True)
