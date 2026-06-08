import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_pyproject_uses_modern_license_metadata():
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())

    assert pyproject["project"]["license"] == "GPL-3.0-only"
    assert pyproject["project"]["license-files"] == ["LICENSE"]
    assert "setuptools>=77,<80" in pyproject["build-system"]["requires"]


def test_source_distribution_manifest_prunes_non_runtime_trees():
    manifest = (REPO_ROOT / "MANIFEST.in").read_text(encoding="utf-8")

    assert "graft capcruncher" in manifest
    for directory in ("prune .github", "prune docs", "prune tests"):
        assert directory in manifest
    for root_file in ("exclude pixi.lock", "exclude uv.lock", "exclude Dockerfile"):
        assert root_file in manifest
    for pattern in ("global-exclude *.py[cod]", "global-exclude __pycache__"):
        assert pattern in manifest


def test_python_optional_features_have_aggregate_extra():
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    optional_dependencies = pyproject["project"]["optional-dependencies"]
    all_extra = set(optional_dependencies["all"])

    for extra in ("full", "config", "plot", "hub", "hpc", "differential"):
        assert extra in optional_dependencies
        assert set(optional_dependencies[extra]) <= all_extra


def test_root_environment_is_pyranges1_only():
    environment = (REPO_ROOT / "environment.yml").read_text(encoding="utf-8")

    assert "pyranges1>=1.3,<2" in environment
    assert "pyranges" + ">=1.0,<2" not in environment
    assert "polars>=1.39,<1.42" in environment
    assert "pyarrow>=24,<25" in environment


def test_critical_dependency_bounds_are_aligned_across_manifests():
    manifest_texts = {
        "pyproject.toml": (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"),
        "environment.yml": (REPO_ROOT / "environment.yml").read_text(encoding="utf-8"),
        "pixi.toml": (REPO_ROOT / "pixi.toml").read_text(encoding="utf-8"),
        "workflow environment": (
            REPO_ROOT / "capcruncher/pipeline/workflow/envs/environment.yml"
        ).read_text(encoding="utf-8"),
    }

    required_bounds = {
        "capcruncher-tools": (">=0.2.4", "<0.3.0"),
        "polars": (">=1.39", "<1.42"),
        "pyarrow": (">=24", "<25"),
        "pyranges1": (">=1.3", "<2"),
        "snakemake": (">=9.21", "<10"),
    }

    for manifest_name, manifest_text in manifest_texts.items():
        for package_name, bounds in required_bounds.items():
            assert package_name in manifest_text, (
                f"{package_name} missing from {manifest_name}"
            )
            for bound in bounds:
                assert bound in manifest_text, (
                    f"{package_name} {bound} missing from {manifest_name}"
                )

    # matplotlib and pyyaml appear in environment.yml and pixi.toml but not in
    # pyproject.toml core deps or the workflow environment
    env_and_pixi = {
        k: v for k, v in manifest_texts.items() if k in ("environment.yml", "pixi.toml")
    }
    for manifest_name, manifest_text in env_and_pixi.items():
        assert "matplotlib" in manifest_text, f"matplotlib missing from {manifest_name}"
        assert ">=3.10.9" in manifest_text, (
            f"matplotlib >=3.10.9 bound missing from {manifest_name}"
        )
        assert "pyyaml" in manifest_text, f"pyyaml missing from {manifest_name}"
        assert ">=6" in manifest_text, f"pyyaml >=6 bound missing from {manifest_name}"


def test_install_docs_present_supported_routes_in_priority_order():
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    installation = (REPO_ROOT / "docs/installation.md").read_text(encoding="utf-8")

    assert readme.index("mamba create -n capcruncher") < readme.index("apptainer exec")
    assert readme.index("apptainer exec") < readme.index("docker run")
    assert readme.index("docker run") < readme.index("pip install capcruncher")
    assert "Pixi is used for development and CI reproducibility" in readme

    assert installation.index("## Recommended Native Install") < installation.index(
        "## Highly Recommended: Containers"
    )
    assert installation.index("## Developer Install") < installation.index(
        "## Detailed Conda Setup"
    )
    assert "Pure pip does not install native pipeline tools" in installation


def test_docker_build_context_and_smoke_contract_are_runtime_scoped():
    dockerfile = (REPO_ROOT / "Dockerfile").read_text(encoding="utf-8")
    dockerignore = (REPO_ROOT / ".dockerignore").read_text(encoding="utf-8")
    container_workflow = (
        REPO_ROOT / ".github/workflows/container-build.yml"
    ).read_text(encoding="utf-8")

    assert "COPY --chown=$MAMBA_USER:$MAMBA_USER . ." not in dockerfile
    assert (
        "COPY --chown=$MAMBA_USER:$MAMBA_USER capcruncher ./capcruncher" in dockerfile
    )
    for excluded_path in ("tests", "docs", "pixi.lock", "uv.lock"):
        assert excluded_path in dockerignore

    assert "docker run --rm capcruncher:smoke-test --help" in container_workflow
    assert "--entrypoint apptainer" in container_workflow
    assert "--entrypoint quarto" not in container_workflow


def test_install_method_ci_covers_documented_routes():
    workflow = (REPO_ROOT / ".github/workflows/install-methods.yml").read_text(
        encoding="utf-8"
    )

    assert "name: Python wheel install" in workflow
    assert "uv pip install dist/*.whl" in workflow
    assert "name: Fallback conda environment" in workflow
    assert "environment-file: environment.yml" in workflow
    assert "uv pip install --no-deps dist/*.whl" in workflow
    assert "name: Docker install smoke" in workflow
    assert "docker build -t capcruncher:install-smoke ." in workflow
    assert "--entrypoint apptainer" in workflow
    assert "name: Apptainer install smoke" in workflow
    assert "apptainer exec capcruncher.sif capcruncher --version" in workflow
    # every smoke test must gate on --version before --help
    assert "capcruncher --version" in workflow
