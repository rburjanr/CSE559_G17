"""Build deduplicated human E3-substrate pair table from raw UbiBrowser + UbiNet and
cross-check against the mentor's pre-cleaned workbook."""
import pandas as pd, sys
D="data/"
ub=pd.read_csv(D+"UbiBrowser_known_E3_substrate_records_raw.txt",sep="\t",dtype=str,quoting=3)
un=pd.read_csv(D+"UbiNet_E3_substrate_archive_2024-09-20.txt",sep="\t",dtype=str,quoting=3)
print("UbiBrowser raw rows:",len(ub),"| species:",ub["species"].value_counts().head(6).to_dict())
print("UbiNet raw rows:",len(un),"| E3 org:",un["E3_ORGANISM"].value_counts().head(4).to_dict())
print("UbiBrowser type values:",ub["type"].value_counts().to_dict())
print("UbiBrowser SOURCE values:",ub["SOURCE"].value_counts().to_dict())
# human filters: E3 and substrate both human
ubh=ub[ub["species"].eq("H.sapiens")].copy()
unh=un[un["E3_ORGANISM"].str.startswith("Homo sapiens",na=False) & un["SUB_ORGANISM"].str.startswith("Homo sapiens",na=False)].copy()
print("human UbiBrowser rows:",len(ubh),"| human UbiNet rows (E3 & sub human):",len(unh))
a=ubh.rename(columns={"SwissProt AC (E3)":"e3","SwissProt AC (Substrate)":"sub"})[["e3","sub","SOURCEID"]].assign(db="UbiBrowser")
a=a.rename(columns={"SOURCEID":"pmid"})
b=unh.rename(columns={"E3_AC":"e3","SUB_AC":"sub","PMID":"pmid"})[["e3","sub","pmid"]].assign(db="UbiNet")
allr=pd.concat([a,b]).dropna(subset=["e3","sub"])
allr["e3"]=allr["e3"].str.strip(); allr["sub"]=allr["sub"].str.strip()
def agg(g):
    pm=sorted({p for x in g["pmid"].dropna() for p in str(x).replace(";",",").split(",") if p.strip().isdigit()} )
    return pd.Series({"databases":"; ".join(sorted(set(g["db"]))),"pmids":"; ".join(pm),"n_records":len(g)})
pairs=allr.groupby(["e3","sub"]).apply(agg,include_groups=False).reset_index()
print("deduplicated human pairs:",len(pairs))
pairs.to_csv("results/pairs_dedup_human.csv",index=False)
m=pd.read_excel(D+"E3_substrate_pairs_and_e3_roster.xlsx",sheet_name="Known pairs")
m.columns=["e3","sub","in_roster","databases","pmids"]
mine=set(zip(pairs.e3,pairs["sub"])); theirs=set(zip(m.e3,m["sub"]))
print("mentor pairs:",len(theirs),"| mine:",len(mine),"| both:",len(mine&theirs),"| mine only:",len(mine-theirs),"| mentor only:",len(theirs-mine))
pd.DataFrame(sorted(mine-theirs),columns=["e3","sub"]).to_csv("results/discrepancy_mine_only.csv",index=False)
pd.DataFrame(sorted(theirs-mine),columns=["e3","sub"]).to_csv("results/discrepancy_mentor_only.csv",index=False)
