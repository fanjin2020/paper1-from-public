"""Public-input supplementary scientific kernels."""
from pathlib import Path
import csv, gzip, io, json, os, re, shutil, subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import tarfile
from pipeline.io import require as need, sha, dump as writej
from .common import csvout, readj, reserved_bh as bh, RDIR, need


FOCUS = ['GFAP', 'GSDMD', 'IGFBP5']


STUDIES = {265:'GSE168137',266:'GSE157827',267:'GSE92926',268:'GSE134151',269:'GSE243292'}


PRIMARY_COUNTS = {265:1,266:1,267:2,268:1,269:7}


def geo_root(acc): return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/'


def soft_read(path,with_tables=False):
    """Parse repeated SOFT fields; do not infer study identities from filenames."""
    entities=[];current=None;table=None;table_header=None
    with gzip.open(path,'rt',encoding='utf-8-sig') as f:
        for line in f:
            line=line.rstrip('\r\n')
            if line.startswith('^'):
                kind,value=line[1:].split(' = ',1)
                current={'kind':kind,'id':value,'fields':{},'table':[]};entities.append(current)
                table=None
            elif line.startswith('!') and line.endswith('_table_begin'):
                table=with_tables;table_header=None
            elif line.startswith('!') and line.endswith('_table_end'):table=None
            elif table:
                row=next(csv.reader([line],delimiter='\t',quoting=csv.QUOTE_NONE))
                if table_header is None:table_header=row
                else:current['table'].append(dict(zip(table_header,row)))
            elif line.startswith('!') and ' = ' in line and current is not None:
                k,v=line[1:].split(' = ',1);current['fields'].setdefault(k,[]).append(v)
    return entities


def first(e,key,default=''): return e['fields'].get(key,[default])[0]


def characteristics(e):
    d={}
    for value in e['fields'].get('Sample_characteristics_ch1',[]):
        if ':' in value:
            k,v=value.split(':',1);k=k.strip().lower();v=v.strip()
            need(k not in d or d[k]==v,f'Conflicting characteristic {e["id"]}/{k}')
            d[k]=v
    return d


def source(dl,acc,out,tables=False):
    p=dl.get(geo_root(acc)+'soft/'+acc+'_family.soft.gz')
    es=soft_read(p,tables)
    csvout(out/'source_declarations.csv',[
        dict(entity=e['id'],kind=e['kind'],field=k,value=v)
        for e in es for k,vs in e['fields'].items() for v in vs])
    return es


def unique_matrix(x):
    need(x.index.is_unique and x.columns.is_unique,'Duplicate matrix IDs')
    need(x.shape[0]>100 and x.shape[1]>=1,'Insufficient matrix dimensions')
    a=x.to_numpy(dtype=float)
    need(np.isfinite(a).all(),'Nonfinite matrix: no silent imputation')
    return x


def fit_task(out,name,x,meta,kind,formula,coefficient,fixed,species,rscript,primary=True):
    """Full-background moderated fit; retain every endpoint, even untestable ones."""
    task=out/name;task.mkdir()
    if len(meta):
        need(meta.index.is_unique,'Duplicate analysis units')
        need(set(meta.index)<=set(x.columns),'Metadata/matrix mismatch')
        x=x.loc[:,meta.index];unique_matrix(x)
        need(kind in ('counts','estimated_counts','log2','continuous'),'Unknown matrix scale')
        if kind=='counts':
            values=x.to_numpy(dtype=float)
            need((values>=0).all() and (values==np.floor(values)).all(),
                 'Raw-count model requires nonnegative integer counts; normalized values are not counts')
        elif kind=='estimated_counts':
            need((x.to_numpy(dtype=float)>=0).all(),'Negative estimated counts')
    meta.to_csv(task/'metadata.csv',index_label='unit_id')
    if len(meta):x.to_csv(task/'expression.csv.gz',index_label='gene_symbol',compression='gzip')
    endpoints=fixed.copy()
    endpoints['analysis_gene']=endpoints.human_gene_symbol if species=='human' else endpoints.gene_symbol
    endpoints.to_csv(task/'endpoints.csv',index=False)
    spec=dict(kind=kind,formula=formula,coefficient=coefficient,primary=primary,
              contrast=name,ci='pointwise moderated 95%; not selection adjusted',
              missing_covariates='Never inferred from donor names or expression',min_group=3,min_residual_df=3)
    writej(task/'specification.json',spec)
    engine=RDIR/'external.R'
    cmd=[rscript,'--vanilla',str(engine),str(task),kind,formula,coefficient]
    with (task/'engine.log').open('w',encoding='utf-8') as log:
        run=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    need(run.returncode==0,f'R failed for {name}; see {task / "engine.log"}')
    results=pd.read_csv(task/'fixed_endpoint_results.csv')
    results['task']=name;results['primary_task']=primary
    results['cell_type']=str(meta.cell_type.iloc[0]) if len(meta) and 'cell_type' in meta and meta.cell_type.nunique()==1 else ''
    results['formula']=formula;results['coefficient']=coefficient;results['input_kind']=kind
    results['species']=species;results['n_analysis_units']=len(meta)
    results['analysis_role']='exploratory_hypothesis_generating_not_confirmatory_validation'
    results['primary_task_meaning']='locked_multiplicity_slot_not_validation_eligibility'
    results['age_sex_adjusted']=bool(re.search(r'\bage\w*\b',formula,re.I) and re.search(r'\bsex\w*\b',formula,re.I))
    results['unadjusted_case_only']=formula.replace(' ','')=='~case'
    return results


def collapse_symbols(x,symbols,mode):
    need(len(x)==len(symbols),'Annotation length mismatch')
    x=x.copy();x.index=pd.Index(symbols,name='gene_symbol')
    return x.groupby(level=0,sort=True).sum() if mode=='sum' else x.groupby(level=0,sort=True).mean()


def run265(out,dl,fixed,args):
    es=source(dl,'GSE168137',out);ss=[e for e in es if e['kind']=='SAMPLE']
    raw=pd.read_csv(dl.get(geo_root('GSE168137')+'suppl/GSE168137_countList.txt.gz'),sep='\t',index_col=0)
    processing=sorted({v for e in ss for v in e['fields'].get('Sample_data_processing',[])})
    need(all(any('RSEM' in v for v in e['fields'].get('Sample_data_processing',[])) for e in ss),
         'GSE168137 RSEM source declaration changed; review count scale before fitting')
    values=raw.to_numpy(dtype=float)
    need(np.isfinite(values).all() and (values>=0).all(),'Invalid deposited estimated counts')
    writej(out/'265_input_scale.json',dict(matrix='GSE168137_countList.txt.gz',
        input_kind='estimated_counts',quantifier='RSEM',source_processing_declarations=processing,
        fractional_cells=int((values!=np.floor(values)).sum()),rounded=False,
        boundary='Source-specific deposited gene-count interpretation, not integer UMI counts. '
        'No transcript-length correction reconstructed; upstream production history not re-certified.'))
    def norm(x):return x.replace(';BL6','').lower()
    titles={norm(first(e,'Sample_title')):e for e in ss}
    need(len(titles)==len(ss),'Duplicate GEO titles')
    records=[];column_audit=[]
    for col in raw:
        column_audit.append(dict(column=col,exact_normalized_match=norm(col) in titles,
            used_for_model='_cortex_' in col.lower(),normalization='remove literal ;BL6 then lowercase',
            unmatched_policy='No inferred correction; out-of-scope hippocampus columns retained in audit only'))
        if norm(col) not in titles and '_hippocampus_' in col.lower():continue
        need(norm(col) in titles,f'Unmatched target count column: {col}')
        e=titles[norm(col)];c=characteristics(e)
        m=re.fullmatch(r'(5xFAD|BL6)_(cortex|hippocampus)_(\d+)mon_(female|male)_(\d+)',first(e,'Sample_title'),re.I)
        need(m is not None,f'Unexpected animal title {first(e,"Sample_title")}')
        geno,region,age,sex,animal=m.groups()
        need(c['tissue'].lower()==region.lower() and int(c['time (month)'])==int(age) and c['sex'].lower()==sex.lower(),'Title/metadata disagreement')
        records.append(dict(unit_id=col,gsm=e['id'],animal_id=animal,region=region.lower(),age=int(age),
                            sex=sex.lower(),case=int(geno.lower()=='5xfad')))
    csvout(out/'source_column_identity_audit.csv',column_audit)
    meta=pd.DataFrame(records).set_index('unit_id');meta.to_csv(out/'all_sample_manifest.csv')
    # Inspect repeated regional specimens before selecting the primary cortex analysis.
    discord=meta.groupby('animal_id')[['age','sex','case']].nunique().max(axis=1)>1
    csvout(out/'within_animal_metadata_conflicts.csv',[dict(animal_id=i,conflict=True) for i in discord[discord].index])
    meta=meta[meta.region=='cortex'].copy();need(meta.animal_id.is_unique,'Repeated animal in cortex')
    need(set(meta.age)=={4,8,12,18},'Unexpected age strata')
    prior=pd.read_csv(args.mouse_manifest,dtype=str)
    prior.to_csv(out/'current_mouse_manifest.csv',index=False)
    matches=[]
    for col in prior:
        values=set(prior[col].dropna().astype(str))
        for a in sorted(set(meta.animal_id)&values):matches.append(dict(animal_id=a,prior_field=col,
            interpretation='identifier_collision_possible_overlap_not_proof; no automatic independent-cohort claim'))
    csvout(out/'cross_study_identifier_collisions.csv',matches)
    p=dl.get('https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt')
    # The official MGI report has a trailing tab in its records.
    with p.open(encoding='utf-8') as f:
        rr=csv.reader(f,delimiter='\t');header=next(rr);rows=list(rr)
    si=next(i for i,v in enumerate(header) if 'marker symbol' in v)
    ei=next(i for i,v in enumerate(header) if 'Ensembl gene id' in v)
    lookup={}
    for row in rows:
        if len(row)>ei and re.fullmatch(r'ENSMUSG\d+',row[ei]):lookup.setdefault(row[ei],set()).add(row[si])
    symbols=[];audit=[]
    for gene in raw.index.astype(str):
        stable=gene.split('.')[0];hits=lookup.get(stable,set())
        symbol=next(iter(hits)) if len(hits)==1 else '__UNRESOLVED__'+stable
        symbols.append(symbol);audit.append(dict(source_gene=gene,stable_id=stable,symbol=symbol,mapping_count=len(hits)))
    csvout(out/'gene_mapping.csv',audit)
    x=collapse_symbols(raw.loc[:,meta.index],symbols,'sum')
    results=[fit_task(out,'cortex_overall',x,meta,'estimated_counts','~factor(age)+factor(sex)+case','case',fixed,'mouse',args.rscript)]
    for age in [4,8,12,18]:
        m=meta[meta.age==age]
        results.append(fit_task(out,f'cortex_age_{age}',x,m,'estimated_counts','~factor(sex)+case','case',fixed,'mouse',args.rscript,False))
    return pd.concat(results,ignore_index=True), 'bounded_analysis_completed_independence_unresolved'


def run267(out,dl,fixed,args):
    es=source(dl,'GSE92926',out,True)
    platforms=[e for e in es if e['kind']=='PLATFORM'];need(len(platforms)==1,'Multiple platforms')
    ann=pd.DataFrame(platforms[0]['table']);need({'ID','gene_assignment'}<=set(ann),'Missing platform annotation')
    symbols={};audit=[]
    for row in ann.to_dict('records'):
        names=set()
        for assignment in row['gene_assignment'].split(' /// '):
            tokens=assignment.split(' // ')
            if len(tokens)>1 and tokens[1].strip() not in ['','---']:names.add(tokens[1].strip())
        if len(names)==1:symbols[row['ID']]=next(iter(names))
        audit.append(dict(probe=row['ID'],symbols=';'.join(sorted(names)),retained=len(names)==1))
    csvout(out/'all_probe_mapping.csv',audit)
    columns={};records=[]
    group_map={'APPNL-G-F':'KI','APPWT':'KI_control','3xTg-AD-H':'3xTg','non-Tg':'3xTg_control'}
    for e in es:
        if e['kind']!='SAMPLE':continue
        title=first(e,'Sample_title');label=title.split(',')[0]
        need(label in group_map,f'Unknown group {title}')
        need(any('RMA' in v for v in e['fields'].get('Sample_data_processing',[])),'Missing RMA declaration')
        t=pd.DataFrame(e['table']);need(t.ID_REF.is_unique,'Duplicate probe ID')
        columns[e['id']]=pd.Series(pd.to_numeric(t.VALUE).to_numpy(),index=t.ID_REF)
        records.append(dict(unit_id=e['id'],title=title,group=group_map[label],case=int(label in ['APPNL-G-F','3xTg-AD-H'])))
    raw=pd.DataFrame(columns);ids=[i for i in raw.index if i in symbols]
    x=collapse_symbols(raw.loc[ids],[symbols[i] for i in ids],'mean')
    meta=pd.DataFrame(records).set_index('unit_id');meta.to_csv(out/'all_sample_manifest.csv')
    res=[]
    for group in ['KI','3xTg']:
        m=meta[meta.group.isin([group,group+'_control'])]
        res.append(fit_task(out,group+'_vs_matched_control',x,m,'log2','~case','case',fixed,'mouse',args.rscript))
    return pd.concat(res,ignore_index=True),'bounded_analysis_completed_independence_unresolved'


def run268(out,dl,fixed,args):
    es=source(dl,'GSE134151',out)
    p=dl.get(geo_root('GSE134151')+'suppl/GSE134151_RPKM_Values.xlsx')
    raw=pd.read_excel(p)
    csvout(out/'workbook_column_mapping.csv',[dict(original=str(c),trimmed=str(c).strip()) for c in raw.columns])
    raw.columns=[str(c).strip() for c in raw.columns]
    need(raw.columns.is_unique and {'Gene number','Gene'}<=set(raw),'Workbook columns changed')
    records=[]
    for e in es:
        if e['kind']!='SAMPLE':continue
        c=characteristics(e);title=first(e,'Sample_title')
        need(title in raw.columns,f'Unmatched workbook sample {title}')
        need(c['genotype'].lower() in ['5xfad','wild-type'] and c['treatment'] in ['Vehicle','PLX5622'],'Unknown experimental group')
        records.append(dict(unit_id=title,gsm=e['id'],region=c['tissue region'],
                            case=int(c['genotype'].lower()=='5xfad'),treat=int(c['treatment']=='PLX5622')))
    meta=pd.DataFrame(records).set_index('unit_id');meta.to_csv(out/'all_sample_manifest.csv')
    m=meta[meta.region.str.lower()=='cortex'].copy();need(len(m)==16,'Expected 16 cortex samples')
    need(raw['Gene number'].notna().all() and raw['Gene number'].is_unique,'Missing/duplicate source gene identifiers')
    annotation=raw[['Gene number','Gene']].copy()
    annotation['analysis_symbol']=raw.Gene.fillna('__UNANNOTATED_ENTREZ__'+raw['Gene number'].astype(str))
    annotation.to_csv(out/'source_gene_annotation.csv',index=False)
    raw['Gene']=annotation.analysis_symbol
    x=raw.loc[:,m.index].apply(pd.to_numeric)
    # Duplicated symbols are not summed across possibly overlapping RPKM features.
    valid=~raw.Gene.duplicated(keep=False)
    csvout(out/'ambiguous_gene_rows.csv',[dict(symbol=s,excluded=True) for s in raw.loc[~valid,'Gene']])
    x=x.loc[valid];x.index=raw.loc[valid,'Gene'].astype(str)
    need((x.to_numpy()>=0).all(),'Negative RPKM')
    x=np.log2(x+1)
    results=[fit_task(out,'cortex_genotype_by_depletion',x,m,'continuous','~case*treat','case:treat',fixed,'mouse',args.rscript)]
    for t in [0,1]:
        results.append(fit_task(out,f'cortex_genotype_treatment_{t}',x,m[m.treat==t],'continuous','~case','case',fixed,'mouse',args.rscript,False))
    return pd.concat(results,ignore_index=True),'bounded_analysis_completed_not_gene_specific_causality'


def tar_frame(tar,member):
    f=tar.extractfile(member);need(f is not None,'Missing archive entry')
    with gzip.GzipFile(fileobj=f) as g:return pd.read_csv(g,sep='\t',header=None,dtype=str,keep_default_na=False)


def review_266_diagnoses(es,out):
    ss=[e for e in es if e['kind']=='SAMPLE']
    rows=[]
    for e in ss:
        raw=characteristics(e).get('diagnosis','')
        rows.append(dict(gsm=e['id'],donor=first(e,'Sample_title'),diagnosis_raw=raw,
                         diagnosis_group={'AD':'AD','healthy control':'NC'}.get(raw,'UNRECOGNIZED'),
                         case={'AD':1,'healthy control':0}.get(raw,'')))
    csvout(out/'266_diagnosis_register.csv',rows)
    bad=[r for r in rows if r['diagnosis_group']=='UNRECOGNIZED']
    need(not bad,f'Unexpected source diagnosis; inspect 266_diagnosis_register.csv: {bad}')
    need(len(rows)==21 and len({r['gsm'] for r in rows})==21 and len({r['donor'] for r in rows})==21,
         'Expected 21 unique GEO samples and donor labels; source changed')
    need(sum(r['case']==1 for r in rows)==12 and sum(r['case']==0 for r in rows)==9,
         'Expected source diagnosis counts: 12 AD and 9 healthy control')
    return rows


def run266(out,dl,fixed,args):
    from scipy.io import mmread
    es=source(dl,'GSE157827',out)
    samples=review_266_diagnoses(es,out)
    p=dl.get(geo_root('GSE157827')+'suppl/GSE157827_RAW.tar')
    columns={};records=[];seen_features=None;qc=[]
    with tarfile.open(p,'r:') as tar:
        members=tar.getmembers();need(all(m.isfile() or m.isdir() for m in members),'Archive contains special entries')
        for sample in samples:
            gsm=sample['gsm'];title=sample['donor']
            selected={}
            for suffix in ['features.tsv.gz','barcodes.tsv.gz','matrix.mtx.gz']:
                hits=[m for m in members if Path(m.name).name.startswith(gsm+'_') and m.name.endswith(suffix)]
                need(len(hits)==1,f'{gsm}: expected exactly one {suffix}');selected[suffix]=hits[0]
            genes=tar_frame(tar,selected['features.tsv.gz']);barcodes=tar_frame(tar,selected['barcodes.tsv.gz'])[0]
            need(genes.shape[1]>=2 and genes[0].is_unique and barcodes.is_unique,'Invalid 10x feature/barcode IDs')
            with gzip.GzipFile(fileobj=tar.extractfile(selected['matrix.mtx.gz'])) as g:x=mmread(g,spmatrix=True).tocsr()
            need(x.shape==(len(genes),len(barcodes)),'Matrix/barcode/feature dimension mismatch')
            need(np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all(),'Not nonnegative integer counts')
            if seen_features is None:seen_features=genes.iloc[:,:2].copy()
            else:need(genes.iloc[:,:2].equals(seen_features),'Inconsistent feature dictionaries across donors')
            # Use the deposited filtered nuclei. No outcome-driven cell QC or new annotation.
            good=np.asarray(x.sum(axis=0)).ravel()>0
            need(good.sum()>=50,f'Too few nonzero deposited nuclei: {gsm}')
            groups={'all_recovered_nuclei':np.flatnonzero(good)}
            for cell,idx in groups.items():
                qc.append(dict(gsm=gsm,cell_type=cell,n_cells=len(idx),retained=len(idx)>=50))
                if len(idx)<50:continue
                key=gsm+'__'+cell;columns[key]=np.asarray(x[:,idx].sum(axis=1)).ravel()
                records.append(dict(unit_id=key,gsm=gsm,donor=title,cell_type=cell,case=sample['case'],
                                    diagnosis_raw=sample['diagnosis_raw'],n_cells=len(idx)))
            print(f'266 {gsm}: {int(good.sum())} nuclei',flush=True)
    csvout(out/'cell_count_inclusion.csv',qc)
    raw=pd.DataFrame(columns,index=seen_features[0]);symbols=seen_features[1].astype(str).to_list()
    x=collapse_symbols(raw,symbols,'sum');meta=pd.DataFrame(records).set_index('unit_id')
    meta.to_csv(out/'all_sample_manifest.csv');results=[]
    for i,cell in enumerate(sorted(meta.cell_type.unique())):
        m=meta[meta.cell_type==cell]
        need(m.donor.is_unique,'Repeated donor within cell context')
        results.append(fit_task(out,f'context_{i:02d}_unadjusted',x,m,'counts','~case','case',fixed,'human',args.rscript,
                                primary=(cell=='all_recovered_nuclei')))
    writej(out/'cell_type_scope.json',dict(supplied_annotation_used=False,
        cell_specific_validation_completed=False,age_sex_adjustment_completed=False,
        interpretation='Recovered-nucleus mixture is not tissue bulk and not an independent adjusted human replication.'))
    return pd.concat(results,ignore_index=True),'partial_question_unadjusted_human_context'


def decoded(a):return np.array([x.decode('utf-8') if isinstance(x,bytes) else str(x) for x in a],dtype=object)


def h5_frame(group):
    """Read old and modern AnnData categorical encodings without anndata."""
    result={}
    for k,node in group.items():
        if k=='__categories':continue
        if isinstance(node,h5py.Group):
            need({'codes','categories'}<=set(node),f'Unsupported categorical encoding: {k}')
            codes=node['codes'][:];cats=decoded(node['categories'][:])
            need(((codes>=-1)&(codes<len(cats))).all(),'Invalid category code')
            result[k]=np.array([cats[c] if c>=0 else '' for c in codes],dtype=object)
        elif '__categories' in group and k in group['__categories']:
            codes=node[:];cats=decoded(group['__categories'][k][:])
            need(((codes>=-1)&(codes<len(cats))).all(),'Invalid legacy category code')
            result[k]=np.array([cats[c] if c>=0 else '' for c in codes],dtype=object)
        else:
            a=node[:];result[k]=decoded(a) if a.dtype.kind in 'OSU' else a
    return pd.DataFrame(result)


def h5_pseudobulk(path,out):
    with h5py.File(path,'r') as f:
        obs=h5_frame(f['obs']);var=h5_frame(f['var'])
        need({'sampleID','atscore','final_celltype','apoe','trem2'}<=set(obs),'Unexpected donor/cell metadata schema')
        need('genename' in var,'Missing explicit gene symbols')
        need(obs[['sampleID','atscore','final_celltype','apoe','trem2']].ne('').all().all(),'Missing donor metadata')
        need(set(obs.atscore)<= {'A-T-','A+T-','A+T+'},'Unexpected pathology labels')
        for field in ['atscore','apoe','trem2']:
            need((obs.groupby('sampleID')[field].nunique()==1).all(),f'Donor metadata not constant: {field}')
        obs.groupby(['sampleID','final_celltype']).size().rename('n_cells').reset_index().to_csv(out/'all_donor_cell_counts.csv',index=False)
        donors=obs[['sampleID','atscore','apoe','trem2']].drop_duplicates().set_index('sampleID')
        donors.to_csv(out/'donor_metadata.csv')
        pairs=list(zip(obs.sampleID,obs.final_celltype));keys=sorted(set(pairs));ki={k:i for i,k in enumerate(keys)}
        codes=np.array([ki[k] for k in pairs]);sums=np.zeros((len(keys),len(var)),dtype=np.float64)
        node=f['X'];need(isinstance(node,h5py.Group),'Expected documented raw count CSR X')
        encoding=node.attrs.get('encoding-type','');encoding=encoding.decode() if isinstance(encoding,bytes) else encoding
        need(encoding=='csr_matrix','X must be CSR; no automatic normalized-layer substitution')
        shape=tuple(node.attrs['shape']);need(shape==(len(obs),len(var)),'X metadata dimensions disagree')
        writej(out/'h5_schema.json',dict(obs_fields=list(obs),var_fields=list(var),shape=list(map(int,shape)),selected_matrix='X',encoding=encoding))
        ptr=node['indptr'][:]
        for start in range(0,len(obs),4096):
            stop=min(start+4096,len(obs));lo,hi=int(ptr[start]),int(ptr[stop])
            data=node['data'][lo:hi];idx=node['indices'][lo:hi]
            need(np.isfinite(data).all() and (data>=0).all() and (data==np.floor(data)).all(),'X is not integer UMI counts')
            block=sp.csr_matrix((data,idx,ptr[start:stop+1]-lo),shape=(stop-start,len(var)))
            group=sp.csr_matrix((np.ones(stop-start),(codes[start:stop],np.arange(stop-start))),shape=(len(keys),stop-start))
            sums+=(group@block).toarray()
        meta=[];columns={}
        for j,(donor,cell) in enumerate(keys):
            n=int((codes==j).sum());d=donors.loc[donor]
            if n<50:continue
            uid=donor+'__'+cell;columns[uid]=sums[j]
            need(d.apoe in ['E3/E3','E3/E4','E4/E4'] and d.trem2 in ['WT','R47H'],'Unexpected APOE/TREM2 coding')
            meta.append(dict(unit_id=uid,donor=donor,cell_type=cell,n_cells=n,pathology=d.atscore,
                case=int(d.atscore=='A+T+'),apoe4=d.apoe.count('E4'),trem2=int(d.trem2=='R47H')))
        x=collapse_symbols(pd.DataFrame(columns),var.genename.astype(str).tolist(),'sum')
    return x,pd.DataFrame(meta).set_index('unit_id'),donors


def run269(out,dl,fixed,args):
    es=source(dl,'GSE243292',out)
    dl.get(geo_root('GSE243292')+'suppl/GSE243292_processed_data_readme.txt')
    gz=dl.get(geo_root('GSE243292')+'suppl/GSE243292_ADsnRNAseq_GEO.h5ad.gz')
    h5=dl.root/'GSE243292_ADsnRNAseq_GEO.h5ad';rec=h5.with_suffix('.decompressed.json')
    if h5.exists():need(rec.exists() and readj(rec)['sha256']==sha(h5) and readj(rec)['compressed_sha256']==sha(gz),'Decompressed cache mismatch')
    else:
        tmp=h5.with_suffix('.part')
        with gzip.open(gz,'rb') as src,tmp.open('wb') as dst:shutil.copyfileobj(src,dst,8*1024**2)
        os.replace(tmp,h5);writej(rec,dict(sha256=sha(h5),compressed_sha256=sha(gz)))
    x,meta,donors=h5_pseudobulk(h5,out)
    # Source sample numbers are explicitly documented; no guessed barcode-to-donor mapping.
    expected={}
    for e in es:
        if e['kind']!='SAMPLE':continue
        title=first(e,'Sample_title');m=re.search(r'sample (\d+)$',title)
        need(m is not None,'Unknown GEO sample title')
        expected[m[1]]=characteristics(e)['amyloid tau_(at)_pathology_score']
    need(set(donors.index)==set(expected),'GEO/H5AD donor set mismatch')
    need(all(donors.loc[d,'atscore']==v for d,v in expected.items()),'GEO/H5AD pathology mismatch')
    meta.to_csv(out/'all_sample_manifest.csv')
    results=[]
    need(meta.cell_type.nunique()==7,'Expected seven author cell types; source changed, revise plan before fitting')
    for i,cell in enumerate(sorted(meta.cell_type.unique())):
        m=meta[(meta.cell_type==cell)&meta.pathology.isin(['A+T-','A+T+'])].copy()
        need(m.donor.is_unique,'Repeated donor within cell type')
        results.append(fit_task(out,f'cell_{i:02d}_AD_vs_PA',x,m,'counts','~case','case',fixed,'human',args.rscript))
        terms=[v for v in ['apoe4','trem2'] if m[v].nunique()>1]
        formula='~'+'+'.join(terms+['case'])
        results.append(fit_task(out,f'cell_{i:02d}_AD_vs_PA_genetic_adjusted',x,m,'counts',formula,'case',fixed,'human',args.rscript,False))
    return pd.concat(results,ignore_index=True),'bounded_pathology_context_not_tau_causality'


def adjust_results(result):
    result=result.copy();p=pd.to_numeric(result.P_Value,errors='coerce')
    result['q_context_237']=np.nan;result['q_direction_focus']=np.nan;result['q_direction_237']=np.nan
    for _,ix in result.groupby('task',sort=False).groups.items():
        q=bh(p.loc[ix].fillna(1));result.loc[ix,'q_context_237']=np.where(p.loc[ix].notna(),q,np.nan)
    for primary in [True,False]:
        ix=result.index[result.primary_task==primary]
        if len(ix):
            q=bh(p.loc[ix].fillna(1));result.loc[ix,'q_direction_237']=np.where(p.loc[ix].notna(),q,np.nan)
        focus=result.index[(result.primary_task==primary)&result.human_gene_symbol.isin(FOCUS)]
        if len(focus):
            q=bh(p.loc[focus].fillna(1));result.loc[focus,'q_direction_focus']=np.where(p.loc[focus].notna(),q,np.nan)
    result['interpretation']='targeted_exploration_not_independent_validation_certification'
    return result
