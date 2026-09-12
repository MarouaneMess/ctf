from sage.all import *
import json, logging, sys, traceback
sys.path.insert(0, "/tmp/cuso/src")
import cuso

logging.basicConfig(level=logging.INFO)

with open("tmp_onehit/leaks.json") as f:
    d=json.load(f)
p=ZZ(d["p"]); n=int(d["bits"]); k=int(d["keep"])
low=n-k
X=ZZ(1)<<(low-1)
hs=[ZZ(q["y"]) for q in d["leaks"]]
A=[(h<<low)+X for h in hs]

R=PolynomialRing(ZZ, names=[f"e{i}" for i in range(8)])
e=list(R.gens())
rels=[]
for i in range(1,8):
    for j in range(i+1,8):
        f=(i-j)*(A[i]+e[i])*(A[j]+e[j]) + j*(A[0]+e[0])*(A[j]+e[j]) - i*(A[0]+e[0])*(A[i]+e[i])
        rels.append(R(f))

bounds={x:(-X, X) for x in e}
print(f"[+] pbits={p.nbits()} low={low} Xbits={X.nbits()} rels={len(rels)}", flush=True)

modes=[
    dict(use_graph_optimization=True, use_intermediate_sizes=True),
    dict(use_graph_optimization=False, use_intermediate_sizes=True),
]
roots=[]
for kw in modes:
    try:
        print("[+] trying",kw,flush=True)
        roots=cuso.find_small_roots(relations=rels,bounds=bounds,modulus=p,**kw)
        print("[+] roots",roots,flush=True)
        if roots: break
    except Exception:
        traceback.print_exc()

if not roots:
    raise SystemExit("no roots")

for sol in roots:
    try:
        es=[ZZ(sol[e[i]]) for i in range(8)]
        z=[(A[i]+es[i])%p for i in range(8)]
        if not all((z[i]>>(n-k)) == hs[i] for i in range(8)): continue
        u=(z[1]*inverse_mod(z[0]-z[1],p))%p
        v=(u*z[0])%p
        if all((v*inverse_mod((u+i)%p,p))%p == z[i] for i in range(8)):
            with open("tmp_onehit/answer.json","w") as f:
                json.dump({"p":int(p),"u":int(u),"v":int(v)},f)
            print("[+] SOLVED",u,v,flush=True)
            break
    except Exception:
        traceback.print_exc()
else:
    raise SystemExit("roots did not validate")
