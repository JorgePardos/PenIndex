"""Run provenance: records exactly how a set of output files was produced,
so results can be traced back to the exact config and package versions that
made them - the kind of record a manuscript's methods section (and a
reviewer) implicitly expects to be reconstructable, matching the "data
availability" expectations set by the source literature (e.g. Echeverria &
Alvarez, Chem. Sci. 2023: raw data for every figure provided in the ESI).
"""

import json
import platform
import sys
from datetime import datetime, timezone
from typing import Any, Dict, List

import penindex


def build_manifest(config: Dict[str, Any], config_path: str, output_files: List[str]) -> Dict[str, Any]:
    """Builds a run manifest: UTC timestamp, penindex/Python/platform
    versions, the versions of the key scientific packages involved, the
    exact resolved config used, and the list of files this run wrote."""
    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'penindex_version': penindex.__version__,
        'python_version': sys.version,
        'platform': platform.platform(),
        'package_versions': _key_package_versions(),
        'config_path': config_path,
        'config': config,
        'output_files': output_files,
    }


def _key_package_versions() -> Dict[str, str]:
    versions = {}
    for module_name in ('numpy', 'pandas', 'scipy', 'MDAnalysis', 'matplotlib', 'yaml'):
        try:
            module = __import__(module_name)
            versions[module_name] = getattr(module, '__version__', 'unknown')
        except ImportError:
            versions[module_name] = 'not installed'
    return versions


def write_manifest(manifest: Dict[str, Any], out_file: str) -> None:
    with open(out_file, 'w') as f:
        json.dump(manifest, f, indent=2, default=str)
