import itertools, sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from gadget_search import rref_subspaces
t=int(sys.argv[1]); a,b=map(int,sys.argv[2].split(':')); d=1; n=t+d
w8=np.exp(1j*np.pi/4); tau1=np.array([1,w8])/np.sqrt(2)
tau=np.ones(1)
for _ in range(t): tau=np.kron(tau,tau1)
target=np.array([np.sqrt(a),np.sqrt(b)]); target/=np.linalg.norm(target)
# single-qubit Cliffords
H=np.array([[1,1],[1,-1]])/np.sqrt(2); S=np.diag([1,1j])
gens=[H,S]; cl=[np.eye(2)]
frontier=[np.eye(2)]
def key(U):
    i=np.argmax(np.abs(U.ravel())>1e-9); U=U/ (U.ravel()[i]/abs(U.ravel()[i])); return tuple(np.round(U.ravel(),6))
seen={key(np.eye(2))}
while frontier:
    nf=[]
    for U in frontier:
        for g in gens:
            V=g@U; k=key(V)
            if k not in seen: seen.add(k); cl.append(V); nf.append(V)
    frontier=nf
cltargets=np.array([C@target for C in cl])   # 24
P1={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
paulis=[]
for lab in itertools.product('IXYZ',repeat=t):
    M=np.ones((1,1))
    for l in lab: M=np.kron(M,P1[l])
    paulis.append(M)
sols=[]
weights=1<<np.arange(n-1,-1,-1)
for k in range(n+1):
    ys=np.array(list(itertools.product((0,1),repeat=k)),dtype=np.int64).reshape(2**k,k)
    pairs=list(itertools.combinations(range(k),2))
    A=np.array(list(itertools.product(range(4),repeat=k)),dtype=np.int64).reshape(4**k,k)
    lin=(A@ys.T)%4
    if pairs:
        Qs=np.array(list(itertools.product((0,1),repeat=len(pairs))),dtype=np.int64)
        yy=np.array([[y[a_]*y[b_] for (a_,b_) in pairs] for y in ys],dtype=np.int64)
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
                ov=np.abs(out@target.conj())**2
                for m in np.where((nrm>1e-9)&(ov>(1-1e-9)*nrm))[0]: sols.append(M[m].copy())
print('solutions',len(sols))
best=0; bestinfo=None
for M in sols:
    G=M.T   # (2^d, 2^t) operator input->output
    s=np.linalg.svd(G,compute_uv=False); G=G/s.max()
    Pi=G.conj().T@G
    classes={}
    for P in paulis:
        Q=P@Pi@P
        kk=tuple(np.round(Q.ravel(),6))
        if kk not in classes: classes[kk]=P
    reps=list(classes.values())
    tot=sum(P@Pi@P for P in reps)
    assert np.allclose(tot,np.eye(2**t)), 'not complete'
    p=0; pall=0
    for P in reps:
        o=G@P@tau; pn=np.vdot(o,o).real; pall+=pn
        if pn<1e-12: continue
        if np.max(np.abs(cltargets.conj()@o)**2)>(1-1e-9)*pn: p+=pn
    if p>best+1e-12: best=p; bestinfo=(len(reps),pall); bestM=M
print('target',a,':',b,'t',t,'best adaptive p',round(best,6),'outcome classes',bestinfo)

import json
nz=[(i,complex(v)) for i,v in enumerate(bestM.reshape(-1)) if abs(v)>1e-12]
json.dump({'t':t,'split':f'{a}:{b}','p':best,'qubit_order':'first t qubits = |T> inputs, last qubit = output',
           'amplitudes':{format(i,'0%db'%n):[round(v.real,10),round(v.imag,10)] for i,v in nz}},
          open(f'gadget_{a}to{b}_t{t}.json','w'),indent=1)
print('saved gadget tensor')
