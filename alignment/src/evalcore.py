import numpy as np, pandas as pd
def load():
    idx=pd.read_csv("results/e3_index.csv")
    pairs=pd.read_excel("data/E3_substrate_pairs_and_e3_roster.xlsx",sheet_name="Known pairs")
    pairs.columns=["e3","sub","in_roster","databases","pmids"]
    pairs=pairs[pairs.in_roster=="Yes"]
    # collapse duplicate accessions (TRIM49D1/D2) into one entity
    idx["entity"]=idx["UniProt accession"].map({a:i for i,a in reversed(list(enumerate(idx["UniProt accession"])))})
    keep=idx.index[idx.entity==idx.index]
    return idx,pairs,keep
def coarse_class(row):
    if row["E3-ome role"]=="receptor": return str(row["Receptor family"])
    c=str(row["E3-ome class"]).split(",")[0].strip()
    for k,v in [("RING-between-RING","RBR"),("HECT","HECT"),("Degenerate RING","Degenerate RING"),("Atypical","Atypical")]:
        if c.startswith(k): return v
    return "RING" if c.startswith("RING") else c
def topk_neighbors(sim,kmax):
    s=sim.copy(); np.fill_diagonal(s,-np.inf)
    order=np.argsort(-s,axis=1,kind="stable")[:,:kmax]
    return order
def share_matrix(sets):
    n=len(sets); S=np.zeros((n,n),dtype=bool)
    inv={}
    for i,s in enumerate(sets):
        for x in s: inv.setdefault(x,[]).append(i)
    for lst in inv.values():
        for a in lst: S[a,lst]=True
    np.fill_diagonal(S,False); return S
def metrics(nb,S,ks):
    out={}
    for k in ks:
        sh=np.take_along_axis(S,nb[:,:k],axis=1)
        out[f"prec@{k}"]=sh.mean(); out[f"hit@{k}"]=sh.any(axis=1).mean()
    return out
