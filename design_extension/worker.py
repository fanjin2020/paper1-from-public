import json
import sys
from pathlib import Path
from .study import gaussian_job

if __name__ == '__main__':
    job = json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))
    gaussian_job(job, Path(job['output']))
