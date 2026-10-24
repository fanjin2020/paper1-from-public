"""Shared I/O and public-download validation for supplementary analyses."""
import csv
from datetime import datetime, timezone
import json
import os
import re
from pathlib import Path
import shutil
from urllib.parse import unquote, urlparse

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from pipeline.io import download, dump, sha

ROOT = Path(__file__).resolve().parents[2]
RDIR = ROOT / 'analysis' / 'supplementary'
digest = sha


class Blocked(ValueError):
    """An explicitly specified data or estimability condition is not satisfied."""


def need(ok, message):
    if not ok:
        raise Blocked(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def readj(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def csvout(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row)) or ['status']
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def save(path, rows, columns=None):
    frame = rows.copy() if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows, columns=columns)
    frame.to_csv(path, index=False)
    return frame


def reserved_bh(values):
    p = np.asarray(values, dtype=float)
    need(p.ndim == 1 and not np.isinf(p).any(), 'Invalid probability vector')
    valid = np.isfinite(p)
    need(((p[valid] >= 0) & (p[valid] <= 1)).all(), 'Invalid probability')
    if not len(p):
        return p.copy()
    adjusted = false_discovery_control(np.where(valid, p, 1.0), method='bh')
    adjusted[~valid] = np.nan
    return adjusted



def discovery_direction_comparison(frame):
    """Interaction signs are not comparisons with the discovery main effect."""
    z = frame.copy()
    interaction = z.coefficient.astype('string').fillna('').str.contains(':', regex=False)
    direction = z.discovery_direction.astype('string').str.strip().str.lower()
    available = np.isfinite(z.logFC) & np.isfinite(z.P_Value)
    known = direction.isin(['up', 'down'])
    comparable = ~interaction & available & known
    z['direction_comparison_status'] = 'comparable'
    z.loc[~known, 'direction_comparison_status'] = 'unknown_discovery_direction'
    z.loc[~available, 'direction_comparison_status'] = 'unavailable_estimate'
    z.loc[interaction, 'direction_comparison_status'] = 'not_applicable_interaction'
    z['same_as_discovery_direction'] = pd.Series(pd.NA, index=z.index, dtype='boolean')
    same = (direction.eq('up') & z.logFC.gt(0)) | (direction.eq('down') & z.logFC.lt(0))
    z.loc[comparable, 'same_as_discovery_direction'] = same.loc[comparable]
    return z


def permitted_download(url):
    if url in {
        'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350&targ=all&view=brief&form=text',
        'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&view=full&form=text',
    }:
        return True
    p = urlparse(url)
    if p.scheme != 'https' or p.query or p.fragment:
        return False
    if p.netloc == 'ftp.ncbi.nlm.nih.gov':
        if re.fullmatch(r'/geo/samples/GSM\d+nnn/GSM\d+/suppl/[^/]+\.(?:h5|gz)', p.path):
            return True
        accessions = ('GSE1297', 'GSE48350', 'GSE168137', 'GSE157827', 'GSE92926',
                      'GSE134151', 'GSE243292', 'GSE167490', 'GSE222007', 'GSE118553')
        return (p.path.startswith('/geo/series/') and any('/' + x + '/' in p.path for x in accessions)
                or p.path in ('/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz',
                              '/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz'))
    return url in {
        'https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt',
        'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt',
        'https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt',
    }


class Download:
    """Import only the requested URL's verified download, never a prepared matrix."""
    def __init__(self, root, caches=(), offline=False):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.caches = [Path(p).expanduser().resolve() for p in caches]
        for cache in self.caches:
            need(cache.is_dir(), 'Missing download cache: ' + str(cache))
            need(not self.root.is_relative_to(cache) and not cache.is_relative_to(self.root),
                 'Download cache must be separate from the new supplementary output')
        self.offline = offline
        self.receipts = []
        self._index = None

    def cache_index(self):
        if self._index is not None:
            return self._index
        self._index = {}
        for cache in self.caches:
            for suffix in ('.download.json', '.receipt.json'):
                for path in cache.rglob('*' + suffix):
                    need(path.resolve().is_relative_to(cache), 'Receipt escapes download cache')
                    try:
                        info = readj(path)
                    except (OSError, ValueError):
                        continue
                    url = info.get('url') if isinstance(info, dict) else None
                    if isinstance(url, str) and permitted_download(url):
                        payload = path.with_name(path.name.removesuffix(suffix))
                        if payload.is_file():
                            self._index.setdefault(url, []).append((cache, payload, info))
        return self._index

    def get(self, url, name=None):
        need(permitted_download(url), 'Unregistered public input URL: ' + url)
        name = name or unquote(urlparse(url).path.rsplit('/', 1)[1])
        need(name and Path(name).name == name and name not in ('.', '..') and '\\' not in name,
             'Unsafe download name')
        target = self.root / name
        if not target.exists():
            candidates = self.cache_index().get(url, [])
            hashes = set()
            chosen = None
            for cache, path, info in candidates:
                need(path.resolve().is_relative_to(cache), 'Payload escapes download cache')
                actual = sha(path)
                need(info.get('sha256') == actual and info.get('bytes', path.stat().st_size) == path.stat().st_size,
                     'Damaged or incomplete raw download receipt: ' + str(path))
                hashes.add(actual)
                chosen = path, dict(info, bytes=path.stat().st_size)
            need(len(hashes) <= 1, 'Conflicting public download versions: ' + url)
            if chosen is not None:
                path, info = chosen
                try:
                    os.link(path, target)
                except OSError:
                    shutil.copy2(path, target)
                dump(target.with_name(name + '.download.json'), info)
                print('REUSE PUBLIC DOWNLOAD', path, flush=True)
        need(not self.offline or target.exists(), 'Offline public download missing: ' + url)
        download(url, target)
        self.receipts.append(readj(target.with_name(name + '.download.json')))
        return target


class CurrentFiles:
    """Access artifacts generated within this supplementary stage only."""
    def __init__(self, root):
        self.root = Path(root).resolve()

    def get(self, name):
        path = (self.root / name).resolve()
        need(path.is_relative_to(self.root) and path.is_file(), 'Missing current-run artifact: ' + name)
        return path
