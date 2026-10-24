import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import collect_results


class CollectionTests(unittest.TestCase):
    def fixture(self,root):
        (root/'reporting').mkdir(parents=True);(root/'figures').mkdir()
        (root/'reporting/status.json').write_text(json.dumps({'status':'RESULT_EXPORTS_COMPLETE','run_id':'test'}))
        (root/'state.json').write_text(json.dumps({'status':'INTEGRATED_ANALYSES_COMPLETE','run_id':'test'}))
        for name in ['software.json','fresh_run.json','paper_coverage.json']:(root/name).write_text('{}')
        data=b'gene,p\nTEST,0.5\n';(root/'reporting/test.csv').write_bytes(data)
        (root/'reporting/source_index.csv').write_text('file,bytes,sha256\ntest.csv,'+str(len(data))+','+hashlib.sha256(data).hexdigest()+'\n')
        for name in ['matched/results/analysis_summary.json','supplementary/python.json',
                     'supplementary/status.json','benchmark/run_state.json']:
            path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('{}')
        (root/'core/raw').mkdir(parents=True);(root/'core/raw/counts.h5').write_text('NOT_AN_EXPORT')

    def test_collection_never_copies_raw(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'run';out=Path(d)/'export';self.fixture(root)
            self.assertEqual(collect_results.main(['--run',str(root),'--output',str(out)]),0)
            self.assertFalse((out/'core').exists())
            self.assertTrue((out/'reporting/test.csv').exists())
            with self.assertRaisesRegex(ValueError,'empty'):
                collect_results.main(['--run',str(root),'--output',str(out)])

    def test_collection_rejects_edited_result(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'run';self.fixture(root)
            (root/'reporting/test.csv').write_text('changed')
            with self.assertRaisesRegex(ValueError,'changed'):
                collect_results.main(['--run',str(root),'--output',str(Path(d)/'export')])


if __name__=='__main__':unittest.main()
