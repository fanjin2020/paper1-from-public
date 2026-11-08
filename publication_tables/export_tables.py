"""Export checksum-verified publication tables without fitting models."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def verify():
    manifest=json.loads((ROOT/'manifest.json').read_text(encoding='utf-8'))
    for row in manifest:
        data=(ROOT/row['file']).read_bytes()
        if hashlib.sha256(data).hexdigest()!=row['sha256']:
            raise ValueError('Checksum mismatch: '+row['file'])
        values=list(csv.reader(io.StringIO(gzip.decompress(data).decode('utf-8'))))
        if len(values)!=row['rows']+1 or any(len(r)!=row['columns'] for r in values):
            raise ValueError('Table dimensions differ: '+row['section'])
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    manifest=verify()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Output must be new or empty')
    args.output.mkdir(parents=True,exist_ok=True)
    for row in manifest:
        (args.output/(row['section']+'.csv')).write_bytes(
            gzip.decompress((ROOT/row['file']).read_bytes()))
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Exported',len(manifest),'sections to',args.output.resolve())


if __name__=='__main__':main()
