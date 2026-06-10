import os
import sys

import pytest


def pytest_addoption(parser):
    parser.addoption("--cores")


@pytest.fixture(scope="session", autouse=True)
def cores(request):
    return request.config.getoption("--cores")


@pytest.fixture(scope="session")
def capcruncher_test_bin(tmp_path_factory):
    """Inject a ``capcruncher`` shim that resolves to this checkout.

    Snakemake rules shell out to ``capcruncher``, so tests need a PATH entry
    that points at the local source rather than any installed entry point.
    All other pipeline tools (flash2, multiqc, split, …) are taken from the
    pixi environment PATH unchanged.
    """

    repo_root = os.path.dirname(__file__)
    bin_dir = tmp_path_factory.mktemp("capcruncher_test_bin")

    capcruncher = bin_dir / "capcruncher"
    capcruncher.write_text(
        f"""#!/usr/bin/env python
import sys

sys.path.insert(0, {repo_root!r})

from capcruncher.cli import cli

sys.exit(cli())
""",
        encoding="utf-8",
    )

    capcruncher.chmod(0o755)

    return bin_dir


@pytest.fixture(scope="session")
def capcruncher_subprocess_env(capcruncher_test_bin):
    repo_root = os.path.dirname(__file__)
    pythonpath = os.environ.get("PYTHONPATH")
    if pythonpath:
        pythonpath = f"{repo_root}{os.pathsep}{pythonpath}"
    else:
        pythonpath = repo_root

    env_bin = os.path.dirname(sys.executable)
    base_path = os.environ.get("PATH", "")
    path = os.pathsep.join(
        filter(None, [str(capcruncher_test_bin), env_bin, base_path])
    )

    return {
        **os.environ,
        "PATH": path,
        "PYTHONPATH": pythonpath,
    }


class MockFastqRecord:
    """Testing class used to supply a pysam FastqProxy like object"""

    def __init__(self, name, sequence, quality):
        self.name = name
        self.sequence = sequence
        self.quality = quality
        self.comment = ""

    def __repr__(self) -> str:
        return "|".join([self.name, self.sequence, "+", self.quality])


class MockFastaRecord:
    """Testing class used to supply a pysam FastqProxy like object"""

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __repr__(self) -> str:
        return f">{self.name}\n{self.sequence}\n"
