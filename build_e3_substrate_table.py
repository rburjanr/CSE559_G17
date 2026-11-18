import pandas as pd
from pathlib import Path
from openpyxl import load_workbook

BASE=Path('/mnt/data')
RAW=BASE/'UbiBrowser_known_E3_substrate_records_raw.txt(1).crdownload'
UBINET=BASE/'UbiNet_E3_substrate_archive_2024-09-20.txt(1).crdownload'
CUR=BASE/'unknown.xlsx'
OUT=BASE/'e3_substrate_output'
OUT.mkdir(exist_ok=True)
# Read sources as strings to preserve accessions/PMIDs.
raw=pd.read_csv(RAW, sep='\t', dtype=str, keep_default_na=False)
ub=pd.read_csv(UBINET, sep='\t', dtype=str, keep_default_na=False)

# Human-only raw UbiBrowser records. Proposal specifies human E3s only; substrate records are protein pairs.
raw_h=raw[raw['species'].str.strip().str.lower().isin(['h.sapiens','homo sapiens','human'])].copy()
# UbiNet: require both E3 and substrate organism to be human for a clean human-human pair table.
ub_h=ub[(ub['E3_ORGANISM'].str.contains('Homo sapiens', case=False, na=False)) &
        (ub['SUB_ORGANISM'].str.contains('Homo sapiens', case=False, na=False))].copy()

# Extract mentor curated E3 genes, retaining sheet/class provenance. Exclude explicitly non-protein-coding sheet.
wb=load_workbook(CUR, read_only=True, data_only=True)
cur_rows=[]
for ws in wb.worksheets:
    if ws.title == 'Non-protein coding':
        continue
    rows=ws.iter_rows(values_only=True)
    try: header=next(rows)
    except StopIteration: continue
    # First column is gene symbol/name in these sheets.
    for r in rows:
        if not r: continue
        gene=r[0]
        if gene is None: continue
        gene=str(gene).strip()
        if not gene or gene.lower() in {'gene symbol','gene name','gene'}: continue
        # capture category if a header named Category exists
        category=''
        for i,h in enumerate(header):
            if h is not None and str(h).strip().lower()=='category' and i < len(r) and r[i] is not None:
                category=str(r[i]).strip(); break
        cur_rows.append((gene, ws.title, category))
cur=pd.DataFrame(cur_rows, columns=['e3_gene','curation_sheet','curation_category']).drop_duplicates()
cur_genes=set(cur.e3_gene.str.upper())

# Normalize two sources into common evidence-level schema.
raw_n=pd.DataFrame({
 'e3_accession':raw_h['SwissProt AC (E3)'].str.strip(),
 'e3_id':raw_h['SwissProt ID (E3)'].str.strip(),
 'e3_gene':raw_h['Gene Symbol (E3)'].str.strip(),
 'substrate_accession':raw_h['SwissProt AC (Substrate)'].str.strip(),
 'substrate_id':raw_h['SwissProt ID (Substrate)'].str.strip(),
 'substrate_gene':raw_h['Gene Symbol (Substrate)'].str.strip(),
 'source_database':'UbiBrowser',
 'source_type':raw_h['SOURCE'].str.strip(),
 'pmid_or_source_id':raw_h['SOURCEID'].str.strip(),
 'evidence_type':raw_h['type'].str.strip(),
 'e3_type':raw_h['E3TYPE'].str.strip(),
 'original_record_id':raw_h['NUMBER'].str.strip(),
 'species':raw_h['species'].str.strip()
})
ub_n=pd.DataFrame({
 'e3_accession':ub_h['E3_AC'].str.strip(),
 'e3_id':ub_h['E3_ID'].str.strip(),
 'e3_gene':ub_h['E3_GN'].str.strip(),
 'substrate_accession':ub_h['SUB_AC'].str.strip(),
 'substrate_id':ub_h['SUB_ID'].str.strip(),
 'substrate_gene':ub_h['SUB_GN'].str.strip(),
 'source_database':'UbiNet',
 'source_type':'PMID',
 'pmid_or_source_id':ub_h['PMID'].str.strip(),
 'evidence_type':'',
 'e3_type':'',
 'original_record_id':ub_h['No'].str.strip(),
 'species':ub_h['E3_ORGANISM'].str.strip()
})
evidence=pd.concat([raw_n,ub_n], ignore_index=True)
for c in ['e3_accession','substrate_accession','e3_gene','substrate_gene']:
    evidence[c]=evidence[c].fillna('').astype(str).str.strip()
# Remove rows lacking pair accessions (dedup definition in proposal is by E3 and substrate accession).
evidence=evidence[(evidence.e3_accession!='') & (evidence.substrate_accession!='')].copy()
evidence['in_mentor_curation']=evidence.e3_gene.str.upper().isin(cur_genes)
# Unique pair table: one row per E3/substrate accession, aggregate provenance rather than losing evidence.
def join_unique(s):
    vals=[]
    for x in s:
        x=str(x).strip()
        if x and x not in vals: vals.append(x)
    return ';'.join(vals)
agg={c:join_unique for c in ['e3_id','e3_gene','substrate_id','substrate_gene','source_database','source_type','pmid_or_source_id','evidence_type','e3_type','original_record_id','species']}
pairs=evidence.groupby(['e3_accession','substrate_accession'], as_index=False).agg(agg)
pairs['in_mentor_curation']=pairs.e3_gene.str.split(';').apply(lambda xs:any(str(x).upper() in cur_genes for x in xs))
# Primary analysis table follows project scope: human E3s represented in mentor curation workbook.
pairs_cur=pairs[pairs.in_mentor_curation].copy()

# Summary statistics.
summary=[]
def add(metric,value): summary.append((metric,value))
add('UbiBrowser raw records',len(raw))
add('UbiBrowser human records',len(raw_h))
add('UbiBrowser human unique accession pairs',raw_h[['SwissProt AC (E3)','SwissProt AC (Substrate)']].drop_duplicates().shape[0])
add('UbiNet raw records',len(ub))
add('UbiNet human-human records',len(ub_h))
add('UbiNet human-human unique accession pairs',ub_h[['E3_AC','SUB_AC']].drop_duplicates().shape[0])
add('Combined evidence records',len(evidence))
add('Combined unique accession pairs',len(pairs))
add('Combined unique E3 accessions',pairs.e3_accession.nunique())
add('Combined unique substrate accessions',pairs.substrate_accession.nunique())
add('Pairs whose E3 gene is in mentor curation workbook',len(pairs_cur))
add('Unique E3 accessions after mentor curation match',pairs_cur.e3_accession.nunique())
add('Unique substrate accessions after mentor curation match',pairs_cur.substrate_accession.nunique())
add('Mentor curation unique gene symbols',len(cur_genes))
# source overlap
srcsets=pairs.source_database.str.split(';').apply(set)
add('Unique pairs supported by both UbiBrowser and UbiNet',sum({'UbiBrowser','UbiNet'}.issubset(x) for x in srcsets))
add('Unique pairs UbiBrowser-only',sum(x=={'UbiBrowser'} for x in srcsets))
add('Unique pairs UbiNet-only',sum(x=={'UbiNet'} for x in srcsets))
summary=pd.DataFrame(summary,columns=['metric','value'])

# Validation issues
issues=[]
# Same accession with multiple gene symbols across evidence
for acc,g in evidence.groupby('e3_accession'):
    genes=sorted({x for x in g.e3_gene if x})
    if len(genes)>1: issues.append(('E3 accession maps to multiple gene symbols',acc,';'.join(genes)))
for acc,g in evidence.groupby('substrate_accession'):
    genes=sorted({x for x in g.substrate_gene if x})
    if len(genes)>1: issues.append(('Substrate accession maps to multiple gene symbols',acc,';'.join(genes)))
issues=pd.DataFrame(issues,columns=['issue','accession','values'])
# Save CSVs.
evidence.to_csv(OUT/'e3_substrate_evidence_records.csv',index=False)
pairs.to_csv(OUT/'e3_substrate_pairs_all_human.csv',index=False)
pairs_cur.to_csv(OUT/'e3_substrate_pairs_curated_e3.csv',index=False)
cur.to_csv(OUT/'mentor_e3_curation_genes.csv',index=False)
summary.to_csv(OUT/'validation_summary.csv',index=False)
issues.to_csv(OUT/'validation_issues.csv',index=False)
print(summary.to_string(index=False))
print('\nIssue rows:',len(issues))
print('\nCurated pair preview:')
print(pairs_cur.head(10).to_string(index=False))
