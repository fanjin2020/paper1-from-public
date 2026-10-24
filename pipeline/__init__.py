"""Public downloads, count preparation, models and result export."""
import json
import os
from pathlib import Path
import sys

# Optional existing packages are appended, never put ahead of the selected environment.
for _entry in json.loads(os.environ.get('PAPER1_OPTIONAL_PACKAGE_PATHS', '[]')):
    _path = Path(_entry).expanduser().resolve()
    if not _path.is_dir():
        raise RuntimeError('Existing optional-package path not found: ' + str(_path))
    if str(_path) not in sys.path:
        sys.path.append(str(_path))
