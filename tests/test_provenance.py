import json

from penindex import provenance


def test_build_manifest_contains_expected_keys():
    config = {'system': {'topology': 'earp.pdb'}, 'analysis': {'hbond_enabled': True}}
    manifest = provenance.build_manifest(config, "config.yaml", ["out.csv", "out.pml"])

    assert manifest['config'] == config
    assert manifest['config_path'] == "config.yaml"
    assert manifest['output_files'] == ["out.csv", "out.pml"]
    assert 'timestamp_utc' in manifest
    assert 'penindex_version' in manifest
    assert 'package_versions' in manifest
    for package in ('numpy', 'pandas', 'scipy', 'MDAnalysis', 'matplotlib', 'yaml'):
        assert package in manifest['package_versions']


def test_write_manifest_produces_valid_json(tmp_path):
    manifest = provenance.build_manifest({'a': 1}, "config.yaml", ["out.csv"])
    out_file = tmp_path / "run_manifest.json"

    provenance.write_manifest(manifest, str(out_file))

    with open(out_file) as f:
        loaded = json.load(f)
    assert loaded['output_files'] == ["out.csv"]
