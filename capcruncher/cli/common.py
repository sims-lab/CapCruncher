from importlib import import_module
from typing import Any

HELP_SETTINGS = {"help_option_names": ["-h", "--help"]}


def run_imported(import_path: str, *args: Any, **kwargs: Any) -> None:
    """Import and run a command implementation on demand."""
    module_name, function_name = import_path.rsplit(":", 1)
    module = import_module(module_name)
    getattr(module, function_name)(*args, **kwargs)
