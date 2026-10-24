"""Explicit manuscript coverage; unavailable workflows cannot become successful analyses."""
import json
from pathlib import Path

from .io import require

ROOT = Path(__file__).resolve().parents[1]


def report(root=ROOT):
    data = json.loads((Path(root) / 'inputs/paper_coverage.json').read_text(encoding='utf-8'))
    missing = [item for item in data['analyses'] if item['status'] != 'implemented']
    return {'complete_paper': not missing, 'analyses': data['analyses'],
            'missing': missing, 'input_boundary': data['input_boundary']}


def require_complete(root=ROOT):
    data = report(root)
    require(data['complete_paper'], 'FULL-PAPER REPRODUCTION NOT READY. Missing independently runnable workflows:\n' +
            '\n'.join(f"- {item['id']}: {item['description']}" for item in data['missing']) +
            '\nNo analysis or download was started. A primary-only run is not a full-paper reproduction.')
    return data
