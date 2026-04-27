from __future__ import annotations

from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from typing import Any

import pandas as pd
from loguru import logger

from capcruncher.types import Executor


def iter_count_results(
    viewpoints: Iterable[str],
    count_kwargs: dict[str, Any],
    executor: str,
    n_cores: int,
) -> Iterable[tuple[str, pd.DataFrame]]:
    from capcruncher_tools.count import count_viewpoint_pixels

    if executor == Executor.LOCAL.value:
        for viewpoint in viewpoints:
            yield count_viewpoint_pixels(viewpoint=viewpoint, **count_kwargs)
        return

    if executor == Executor.PROCESS.value:
        process_kwargs: dict[str, Any] = {"max_workers": n_cores}
        try:
            process_kwargs["mp_context"] = get_context("fork")
        except ValueError:
            pass
        try:
            with ProcessPoolExecutor(**process_kwargs) as pool:
                futures = [
                    pool.submit(
                        count_viewpoint_pixels, viewpoint=viewpoint, **count_kwargs
                    )
                    for viewpoint in viewpoints
                ]
                for future in as_completed(futures):
                    yield future.result()
        except PermissionError:
            logger.warning(
                "Process executor is unavailable in this environment; "
                "falling back to local execution"
            )
            for viewpoint in viewpoints:
                yield count_viewpoint_pixels(viewpoint=viewpoint, **count_kwargs)
        return

    if executor == Executor.RAY.value:
        try:
            import ray
        except ImportError as exc:
            raise RuntimeError(
                "Ray executor requested but ray is not installed. "
                "Install capcruncher-tools[ray]."
            ) from exc

        ray.init(num_cpus=n_cores, ignore_reinit_error=True)
        remote_count = ray.remote(count_viewpoint_pixels)
        futures = [
            remote_count.remote(viewpoint=viewpoint, **count_kwargs)
            for viewpoint in viewpoints
        ]
        while futures:
            completed, futures = ray.wait(futures)
            for future in completed:
                yield ray.get(future)
        return

    raise ValueError(f"Unknown executor: {executor}")
