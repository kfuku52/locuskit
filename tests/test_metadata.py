import re
import tomllib
from pathlib import Path

from kffractbias import __version__

ROOT = Path(__file__).resolve().parents[1]


def _declared_version(path: Path) -> str:
    match = re.search(
        r"^version\s*[:=]\s*[\"']?([^\"'\s]+)", path.read_text(encoding="utf-8"), re.MULTILINE
    )
    assert match is not None
    return match.group(1)


def test_version_metadata_stays_synchronized() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    expected = project["version"]
    assert _declared_version(ROOT / "CITATION.cff") == expected
    assert __version__ == expected
    lock = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))
    (locked_project,) = (
        package
        for package in lock["package"]
        if package["name"] == project["name"] and package["source"] == {"editable": "."}
    )
    assert locked_project["version"] == expected


def test_compatibility_metadata_tracks_locuskit() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    shim = tomllib.loads((ROOT / "compat/kffractbias/pyproject.toml").read_text())["project"]
    version = project["version"]
    assert shim["version"] == version
    assert shim["dependencies"] == [f"locuskit=={version}"]
    for extra, requirements in shim["optional-dependencies"].items():
        assert extra in project["optional-dependencies"]
        assert requirements == [f"locuskit[{extra}]=={version}"]
