import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from xml.etree import ElementTree as E
from zipfile import ZipFile

from pipeline import publication_tables as P
from pipeline import workbook as X
from pipeline.io import sha
from pipeline.result_collection import validate_block_mapping

S = {'s':X.NS}


def cell_value(cell):
    if cell.get('t')=='inlineStr':
        return ''.join(cell.itertext())
    return cell.findtext('s:v',default='',namespaces=S)


class WorkbookTests(unittest.TestCase):
    def fixture(self, root, override=None):
        records, mapping = [],[]
        for block in sorted(P.EXPECTED_BLOCKS):
            sid = int(block[1:].split('.')[0]);name=f'S{sid:02}/{block}.csv'
            path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
            header=['gene','value','missing','literal']
            data=[['00123','0.125','','=not_a_formula'],['GFAP','-1e-10','','False']]
            if override and block in override:
                header,data=override[block]
            with path.open('w',encoding='utf-8',newline='') as f:
                writer=csv.writer(f);writer.writerow(header);writer.writerows(data)
            records.append(dict(section=f'S{sid:02}',blocks=block,file=name,rows=len(data),
                bytes=path.stat().st_size,sha256=sha(path),description='fixture'))
            mapping.append(dict(block=block,file=name,description='fixture'))
        self.write_index(root,records)
        with (root/'block_sources.csv').open('w',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(mapping[0]));writer.writeheader();writer.writerows(mapping)
        return records,mapping

    def write_index(self,root,records):
        with (root/'source_index.csv').open('w',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)

    def test_schema_and_legacy_collection(self):
        self.assertEqual(sum(P.BLOCK_COUNTS),112)
        mapping=[{'block':b} for b in P.EXPECTED_BLOCKS]
        self.assertEqual(validate_block_mapping(mapping),112)
        additions={'S6.6','S10.2'}|{f'S20.{i}' for i in range(13,18)}
        self.assertEqual(validate_block_mapping([r for r in mapping if r['block'] not in additions]),105)
        for bad in (mapping[:-1],mapping+[mapping[0]],mapping[:-1]+[{'block':'S99.1'}]):
            with self.assertRaises(ValueError):validate_block_mapping(bad)

    def test_public_sources_and_fixed_release(self):
        data=P.resource_catalog()
        self.assertEqual(len(data),12)
        self.assertEqual(list(data.columns),['resource','analysis_role','relationship','public_source','input_specification'])
        self.assertTrue(data.public_source.str.startswith('https://').all())
        self.assertIn('2026-06-22',data.set_index('resource').loc['SEA-AD','input_specification'])
        self.assertNotIn('download.json',data.to_csv(index=False))

    def test_job_expansion_is_lossless(self):
        jobs=[dict(id=f'job{i}',phase='paired',start=i,end=i+2,spec={'seed':i,'a':[None,True,'quoted']}) for i in range(4)]
        before=json.dumps(jobs,sort_keys=True)
        frame=P.simulation_jobs({'jobs':jobs})
        rebuilt=[dict(id=r.job_id,phase=r.phase,start=r.start,end=r.end,spec=json.loads(r.specification)) for r in frame.itertuples(index=False)]
        self.assertEqual(json.dumps(rebuilt,sort_keys=True),before)
        self.assertEqual(json.dumps(jobs,sort_keys=True),before)
        for bad in ([],jobs+[jobs[0]],[dict(jobs[0],extra=1)]):
            with self.assertRaises(ValueError):P.simulation_jobs({'jobs':bad})

    def test_display_preserves_numbers_and_urls(self):
        header=['path','type','value']
        self.assertEqual(X.display_values('S17.7',header,['/source','str','"/home/u/run/mouse_units.csv"'])[-1],'"mouse_units.csv"')
        self.assertEqual(X.display_values('S20.7',header,['/source','str','"C:\\run\\mouse_units.csv"'])[-1],'"mouse_units.csv"')
        url=['/url','str','"https://example.org/a/b"']
        self.assertEqual(X.display_values('S17.7',header,url),url)
        self.assertEqual(X.display_values('S20.7',header,['/seed','int','20261015'])[-1],'20261015')
        self.assertEqual(X.display_values('S17.7',header,['/jobs','array','x'*100000])[-1],X.JOB_NOTE)

    def test_workbook_catalog_and_all_cells(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture(root)
            path=X.generate(root)
            manifest=json.loads((root/'workbook_manifest.json').read_text())
            self.assertEqual(manifest['blocks'],112)
            self.assertEqual(sha(path),manifest['sha256'])
            with ZipFile(path) as z:
                self.assertIsNone(z.testzip())
                book=E.fromstring(z.read('xl/workbook.xml'))
                self.assertEqual([n.get('name') for n in book.find('s:sheets',S)],[f'S{i:02}' for i in range(1,21)])
                for sid in range(1,21):
                    sheet=E.fromstring(z.read(f'xl/worksheets/sheet{sid}.xml'))
                    cells={c.get('r'):c for c in sheet.findall('.//s:c',S)}
                    for entry in (e for e in manifest['sections'] if e['sheet']==sid):
                        h=entry['header_row']
                        self.assertEqual([cell_value(cells[f'{col}{h}']) for col in 'ABCD'],['gene','value','missing','literal'])
                        self.assertEqual([cell_value(cells[f'{col}{h+1}']) for col in 'ABCD'],['00123','0.125','','=not_a_formula'])
                        self.assertEqual(cells[f'A{h+1}'].get('t'),'inlineStr')
                        self.assertIsNone(cells[f'B{h+1}'].get('t'))
                        self.assertEqual(cell_value(cells[f'D{h+2}']),'0')
                    self.assertEqual(sheet.findall('.//s:f',S),[])
                guide=E.fromstring(z.read('xl/worksheets/sheet1.xml'))
                values={c.get('r'):cell_value(c) for c in guide.findall('.//s:c',S)}
                for row,entry in enumerate(manifest['sections'],6):
                    self.assertEqual(values[f'A{row}'],entry['block'])
                    self.assertEqual(int(values[f'D{row}']),entry['header_row'])

    def test_long_protocol_is_displayed_without_truncating_source(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);value=json.dumps([{'a':'x'*1000}]*300)
            self.fixture(root,{'S17.7':(['path','type','value'],[['/jobs','array',value]])})
            source=root/'S17/S17.7.csv';before=sha(source)
            path=X.generate(root)
            self.assertEqual(sha(source),before)
            with ZipFile(path) as z:
                self.assertIn(X.JOB_NOTE,z.read('xl/worksheets/sheet17.xml').decode())

    def test_changed_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture(root)
            (root/'S01/S1.1.csv').write_text('changed')
            with self.assertRaisesRegex(ValueError,'Changed workbook source'):X.generate(root)

    def test_bad_count_or_overlong_cell_never_replaces_good_workbook(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);records,_=self.fixture(root)
            target=X.generate(root);before=sha(target)
            records[0]['rows']=999;self.write_index(root,records)
            with self.assertRaisesRegex(ValueError,'row count'):X.generate(root)
            self.assertEqual(sha(target),before)
            self.assertFalse((root/'Supplementary_Tables.xlsx.tmp').exists())
        with self.assertRaisesRegex(ValueError,'exceeds limit'):X.cell_xml('A1','x'*32768)

    def test_gzip_and_large_csv_limit_restored(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'data.csv.gz'
            with gzip.open(path,'wt',encoding='utf-8') as f:f.write('value\n'+'x'*200000+'\n')
            old=csv.field_size_limit(131072)
            try:
                with X.rows(path) as reader:self.assertEqual(len(list(reader)[1][0]),200000)
                self.assertEqual(csv.field_size_limit(),131072)
            finally:csv.field_size_limit(old)

    def test_path_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ('../x.csv','C:/x.csv','/tmp/x.csv','..\\x.csv'):
                with self.assertRaises(ValueError):X.source_path(Path(d),name)


if __name__=='__main__':unittest.main()
