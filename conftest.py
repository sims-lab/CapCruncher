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
    """Executables used by subprocess-based pipeline tests.

    The pipeline shells out to ``capcruncher`` from Snakemake rules, so tests need
    a PATH entry that resolves to this checkout rather than any installed conda
    entry point.
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

    gzcat = bin_dir / "gzcat"
    gzcat.write_text(
        """#!/usr/bin/env python
import gzip
import shutil
import sys

with gzip.open(sys.argv[1], "rb") as source:
    shutil.copyfileobj(source, sys.stdout.buffer)
""",
        encoding="utf-8",
    )

    gsplit = bin_dir / "gsplit"
    gsplit.write_text(
        """#!/usr/bin/env python
import gzip
import pathlib
import sys

args = sys.argv[1:]
n_lines = int(args[args.index("-l") + 1])
suffix = ""
for arg in args:
    if arg.startswith("--additional-suffix="):
        suffix = arg.split("=", 1)[1]
prefix = args[-1]

lines = sys.stdin.buffer.readlines()
for part, offset in enumerate(range(0, len(lines), n_lines)):
    output = f"{prefix}{part:02d}{suffix}.gz"
    pathlib.Path(output).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wb") as handle:
        handle.writelines(lines[offset : offset + n_lines])
""",
        encoding="utf-8",
    )

    flash = bin_dir / "flash"
    flash.write_text(
        """#!/usr/bin/env python
import subprocess
import sys

args = sys.argv[1:]
# flash2 does not support --compress-prog-args; strip it and its value
filtered = []
skip_next = False
for arg in args:
    if skip_next:
        skip_next = False
        continue
    if arg == "--compress-prog-args":
        skip_next = True
        continue
    filtered.append(arg)

result = subprocess.run(["flash2", *filtered], check=False)
sys.exit(result.returncode)
""",
        encoding="utf-8",
    )

    multiqc = bin_dir / "multiqc"
    multiqc.write_text(
        """#!/usr/bin/env python
import pathlib
import sys

args = sys.argv[1:]
outdir = pathlib.Path(".")
name = "multiqc_report.html"
for idx, arg in enumerate(args):
    if arg in {"-o", "--outdir"}:
        outdir = pathlib.Path(args[idx + 1])
    elif arg == "-n":
        name = args[idx + 1]

outdir.mkdir(parents=True, exist_ok=True)
(outdir / name).write_text("<html></html>\\n", encoding="utf-8")

data_dir = outdir / "multiqc_data"
data_dir.mkdir(exist_ok=True)
(data_dir / "multiqc_cutadapt.txt").write_text(
    "Sample\\tr_processed\\tr_written\\tr_with_adapters\\n"
    "SAMPLE-A_REP1_part0_1\\t10\\t9\\t1\\n"
    "SAMPLE-A_REP1_part0_2\\t10\\t8\\t2\\n",
    encoding="utf-8",
)
(data_dir / "multiqc_flash_combo_stats.txt").write_text(
    "Sample\\tcombopairs\\tuncombopairs\\n"
    "SAMPLE-A_REP1_part0\\t7\\t3\\n",
    encoding="utf-8",
)
(data_dir / "multiqc_bowtie2.txt").write_text("Sample\\n", encoding="utf-8")
""",
        encoding="utf-8",
    )

    for executable in [capcruncher, flash, gzcat, gsplit, multiqc]:
        executable.chmod(0o755)

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
    path = os.pathsep.join(filter(None, [str(capcruncher_test_bin), env_bin, base_path]))

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
