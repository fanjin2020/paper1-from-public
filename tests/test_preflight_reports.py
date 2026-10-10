"""验证预检查失败的明细会随审阅包一起保留。"""
from pathlib import Path
import tempfile
import unittest
import zipfile

from run import bundle


class PreflightReportTests(unittest.TestCase):
    """不运行R，用小型诊断文件检查审阅包的收录规则。"""

    def test_smoke_details_are_bundled(self):
        """预检查的CSV和警告文本不能只留在服务器日志目录。"""
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            logs = work/'logs'
            logs.mkdir()
            names = ['dream_retry_smoke_table.csv', 'dream_retry_smoke_diagnostics.csv',
                     'dream_retry_smoke_retry_attempts.csv', 'dream_retry_smoke_messages.txt']
            for name in names:
                (logs/name).write_text('test detail\n', encoding='utf-8')
            bundle(work)
            with zipfile.ZipFile(work/'review_bundle.zip') as archive:
                for name in names:
                    self.assertEqual(archive.read(f'logs/{name}'), (logs/name).read_bytes())
