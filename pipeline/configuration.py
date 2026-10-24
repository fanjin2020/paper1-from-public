"""Validate release configuration before creating an analysis directory."""
import json
from pathlib import Path

from pipeline.io import require

VERSION = '2.1.0'
ROOT = Path(__file__).resolve().parents[1]


def load_config(path=None, workers=None):
    path = Path(path) if path is not None else ROOT / 'config.json'
    try:
        config = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError) as exc:
        raise ValueError('Cannot read configuration ' + str(path) + ': ' + str(exc)) from exc
    require(isinstance(config, dict), 'Configuration must be a JSON object: ' + str(path))
    require(config.get('analysis_version') == VERSION,
            f'Configuration version mismatch in {path}: expected {VERSION}, '
            f'found {config.get("analysis_version")!r}. Synchronize the source and config.json '
            'from the same release; no analysis has started.')
    if workers is not None:
        config['workers'] = workers
    require(type(config.get('workers')) is int and config['workers'] >= 1,
            'Configuration workers must be a positive integer')
    require(type(config.get('calibration_repeats')) is int and config['calibration_repeats'] >= 100,
            'Configuration calibration_repeats must be an integer of at least 100')
    require(isinstance(config.get('ebayes_robust'), bool)
            and isinstance(config.get('dream_ebayes_robust'), bool),
            'Configuration must specify boolean ebayes_robust and dream_ebayes_robust')
    return config
