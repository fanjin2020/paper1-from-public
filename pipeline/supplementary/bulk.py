"""GSE48350 public matrix and metadata parsing."""
from collections import Counter, defaultdict
import csv, gzip, math, re
from .common import need, digest, csvout as write_csv
from pipeline.io import dump
URLS = {
 'metadata': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350&targ=all&view=brief&form=text',
 'matrix': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/GSE48350_series_matrix.txt.gz',
 'platform': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&view=full&form=text',
}

REGIONS = {'entorhinal cortex': 'entorhinal_cortex', 'superior frontal gyrus': 'frontal_cortex'}

def read_csv(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def one(fields, key):
    values = fields.get(key, [])
    need(len(values)==1, f'Expected exactly one {key}, got {values}')
    return values[0]

def parse_metadata(path, out):
    samples, series, current = {}, defaultdict(list), None
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('^SAMPLE = '):
            gsm = line.split(' = ',1)[1]
            need(gsm not in samples, 'Duplicate GSM in metadata')
            current = defaultdict(list); samples[gsm] = current
        elif line.startswith('^'):
            current = None
        elif ' = ' in line:
            key, value = line.split(' = ',1)
            if key.startswith('!Series_'):
                series[key].append(value)
            if current is not None and key.startswith('!Sample_'):
                current[key].append(value)
    need(set(series['!Series_sample_id'])==set(samples), 'Series membership mismatch')
    rows, declarations = [], []
    for gsm, fields in samples.items():
        title = one(fields,'!Sample_title')
        chars = defaultdict(list)
        for value in fields.get('!Sample_characteristics_ch1',[]):
            need(': ' in value, f'Unparsed characteristic {gsm}: {value}')
            key, text = value.split(': ',1); chars[key.lower()].append(text)
        region = one(chars,'brain region').lower()
        sex = one(chars,'gender').lower()
        age = float(one(chars,'age (yrs)'))
        need(math.isfinite(age) and 0<age<120 and sex in {'female','male'}, f'Invalid age/sex {gsm}')
        individual = chars.get('individual',[])
        ctl = re.search(r'_indiv(\d+(?:-\d+)?)$',title,re.I)
        ad = re.search(r'_AD_(\d+)$',title)
        if ctl:
            group, ident, basis = 'control', ctl.group(1), 'title_indiv_and_individual_characteristic'
            need(len(individual)==1 and re.fullmatch(r'\s*'+re.escape(ident)+r'(?:,\s*[^,]+)?\s*',individual[0]),
                 f'Control title/individual conflict {gsm}')
            need(any('cognitively intact' in v.lower() for v in fields.get('!Sample_description',[])),
                 f'Control description not confirmed {gsm}')
        elif ad:
            group, ident, basis = 'AD', ad.group(1), 'published_AD_title_label_not_independently_certified_donor'
            need(one(fields,'!Sample_source_name_ch1').endswith('_AD'),f'AD source conflict {gsm}')
        else:
            raise ValueError(f'Cannot establish published group label: {gsm}, {title}')
        # Group namespacing preserves potentially overlapping numbering; no donor equivalence is assumed.
        record = dict(sample_id=gsm,title=title,group=group,age=age,sex=sex,
            source_region=region,region=REGIONS.get(region,''),source_person_label=ident,
            individual_characteristic=';'.join(individual),
            recorded_unit=f'{group}:{ident}',unit_basis=basis,
            platform=one(fields,'!Sample_platform_id'),
            source_series=';'.join(fields.get('!Sample_series_id',[])),
            analysis_include=region in REGIONS,
            exclusion_reason='' if region in REGIONS else 'outside_two_target_regions')
        need(record['platform']=='GPL570', 'Unexpected platform')
        rows.append(record)
        for key in ['!Sample_data_processing','!Sample_description','!Sample_source_name_ch1']:
            declarations.extend(dict(sample_id=gsm,field=key,text=v) for v in fields.get(key,[]))
    selected = [x for x in rows if x['analysis_include']]
    keys = [(x['region'],x['recorded_unit']) for x in selected]
    need(len(keys)==len(set(keys)), 'Repeated source unit within brain region; resolve replicates before fitting')
    units = defaultdict(set)
    for x in rows:
        units[x['recorded_unit']].add((x['age'],x['sex']))
    need(all(len(v)==1 for v in units.values()), 'Conflicting age/sex for a source person label')
    write_csv(out/'187_all_sample_register.csv',rows,list(rows[0]))
    write_csv(out/'187_processing_declarations.csv',declarations,['sample_id','field','text'])
    write_csv(out/'187_selected_samples.csv',selected,list(rows[0]))
    dump(out/'187_series_declarations.json',dict(series))
    return rows, selected

def annotation(path, out):
    rows, active, header = [], False, None
    with path.open(encoding='utf-8-sig',newline='') as f:
        for line in f:
            if line.startswith('!platform_table_begin'):
                active = True; continue
            if line.startswith('!platform_table_end'):
                break
            if not active:
                continue
            values = next(csv.reader([line],delimiter='\t'))
            if header is None:
                header = values
                need({'ID','Gene Symbol'}<=set(header),'GPL570 missing ID/Gene Symbol columns')
                continue
            need(len(values)==len(header),'Malformed GPL570 annotation row')
            d = dict(zip(header,values)); raw = d['Gene Symbol'].strip()
            symbols = sorted(set(s.strip() for s in raw.split('///') if s.strip() not in {'','---'}))
            status = 'unique_symbol' if len(symbols)==1 else 'ambiguous_or_unmapped'
            rows.append(dict(probe_id=d['ID'],raw_symbol=raw,
                             gene_symbol=symbols[0] if status=='unique_symbol' else '',status=status))
    need(rows and len({r['probe_id'] for r in rows})==len(rows),'Invalid annotation table')
    write_csv(out/'187_probe_annotation.csv',rows,list(rows[0]))

def extract_matrix(path, rows, selected, out):
    expected = {x['sample_id'] for x in rows}
    chosen = [x['sample_id'] for x in selected]
    active, header, seen, ended = False, None, set(), False
    with gzip.open(path,'rt',encoding='utf-8-sig',newline='') as src, \
            gzip.open(out/'187_selected_probe_matrix.csv.gz','wt',encoding='utf-8',newline='') as target:
        w = csv.writer(target)
        for line in src:
            if line.startswith('!series_matrix_table_begin'):
                active = True; continue
            if line.startswith('!series_matrix_table_end'):
                ended = True; break
            if not active:
                continue
            values = next(csv.reader([line],delimiter='\t'))
            if header is None:
                header = values
                need(header[0]=='ID_REF' and set(header[1:])==expected,'Expression/metadata sample set mismatch')
                need(len(set(header))==len(header),'Duplicate matrix sample')
                indices = [header.index(s) for s in chosen]
                w.writerow(['probe_id']+chosen); continue
            need(len(values)==len(header) and values[0] not in seen,'Duplicate/malformed expression row')
            seen.add(values[0]); w.writerow([values[0]]+[values[i] for i in indices])
    need(ended and len(seen)>10000,'Truncated/missing expression table')
    return len(seen)
