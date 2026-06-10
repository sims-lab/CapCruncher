from __future__ import annotations

import importlib
import importlib.metadata
import re

CAPCRUNCHER_TOOLS_DISTRIBUTION = "capcruncher-tools"
CAPCRUNCHER_TOOLS_MODULE = "capcruncher_tools"
CAPCRUNCHER_TOOLS_REQUIREMENT = ">=0.2.4,<0.3.0"


class DependencyVersionError(RuntimeError):
    """Raised when an installed runtime dependency is outside the supported range."""


def _module_path(module_name: str) -> str:
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError:
        return "<not importable>"

    return str(getattr(module, "__file__", "<unknown>"))


def _version_satisfies(version: str, requirement: str) -> bool:
    try:
        from packaging.specifiers import SpecifierSet
        from packaging.version import Version
    except ModuleNotFoundError:
        version_match = re.match(r"^(\d+)\.(\d+)\.(\d+)", version)
        if requirement != CAPCRUNCHER_TOOLS_REQUIREMENT or not version_match:
            return False
        release = tuple(int(part) for part in version_match.groups())
        return (0, 2, 4) <= release < (0, 3, 0)

    return Version(version) in SpecifierSet(requirement)


def require_capcruncher_tools() -> str:
    """Return the installed capcruncher-tools version if it is supported."""
    module_path = _module_path(CAPCRUNCHER_TOOLS_MODULE)

    try:
        installed_version = importlib.metadata.version(CAPCRUNCHER_TOOLS_DISTRIBUTION)
    except importlib.metadata.PackageNotFoundError as exc:
        raise DependencyVersionError(
            f"{CAPCRUNCHER_TOOLS_DISTRIBUTION} is required "
            f"({CAPCRUNCHER_TOOLS_REQUIREMENT}) but is not installed. "
            f"Imported module path: {module_path}"
        ) from exc

    if not _version_satisfies(installed_version, CAPCRUNCHER_TOOLS_REQUIREMENT):
        raise DependencyVersionError(
            f"{CAPCRUNCHER_TOOLS_DISTRIBUTION} {CAPCRUNCHER_TOOLS_REQUIREMENT} "
            f"is required, but version {installed_version} is installed. "
            f"Imported module path: {module_path}"
        )

    return installed_version
