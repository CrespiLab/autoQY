"""Bounded, ordered Monte Carlo execution with Windows-safe worker processes."""

from collections import deque
from contextlib import contextmanager
from copy import deepcopy
import os


def available_cpu_count():
    detector = getattr(os, 'process_cpu_count', os.cpu_count)
    count = detector()
    if hasattr(os, 'sched_getaffinity'):
        try:
            count = min(count or 1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return max(1, count or 1)


def cpu_worker_limits():
    cpus = available_cpu_count()
    maximum = max(1, min(61, cpus - 1))
    suggested = min(4, maximum, max(2, cpus // 2))
    return cpus, suggested, maximum


def resolve_workers(requested):
    if isinstance(requested, bool) or not isinstance(requested, int) or not 1 <= requested <= 61:
        raise ValueError('Monte Carlo workers must be an integer from 1 to 61')
    cpus, _, maximum = cpu_worker_limits()
    return min(requested, maximum), cpus


def _initialize_worker():
    # Load both NumPy and SciPy BLAS libraries before setting their limits.
    import numpy
    import scipy.linalg
    from threadpoolctl import threadpool_limits

    global _worker_thread_limit
    _worker_thread_limit = threadpool_limits(limits=1)


@contextmanager
def process_pool(workers):
    if workers == 1:
        yield None
        return
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    try:
        import threadpoolctl
    except ImportError as error:
        raise RuntimeError('Parallel support is unavailable in this installation. '
                           'Update AutoQY, or turn off parallel simulations.') from error
    executor = ProcessPoolExecutor(max_workers=workers,
                                   mp_context=multiprocessing.get_context('spawn'),
                                   initializer=_initialize_worker)
    try:
        yield executor
    finally:
        executor.shutdown(wait=True, cancel_futures=True)


def trial_results(start, stop, workers, prepare, rng, executor):
    """Keep at most workers draws in memory and yield a resumable ordered prefix.

    All randomness stays in the parent. Each saved state belongs to its own
    draw, even though later inputs may already be prepared. The caller bounds
    each block at the next precision check so early stopping does no extra fits.
    """
    from concurrent.futures.process import BrokenProcessPool
    import numpy as np
    from .monte_carlo_draw import fit_trial

    errors = (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError)
    pending = deque()
    next_iteration = start

    def enqueue():
        nonlocal next_iteration
        iteration = next_iteration
        next_iteration += 1
        future, record, failure = None, None, None
        try:
            trial = prepare()
            if executor is None:
                record = fit_trial(trial)
            else:
                future = executor.submit(fit_trial, trial)
        except BrokenProcessPool:
            raise
        except errors as error:
            failure = f'Simulation {iteration + 1}: {error}'
        pending.append((iteration, deepcopy(rng.bit_generator.state), future, record, failure))

    try:
        for _ in range(min(workers, stop - start)):
            enqueue()
        while pending:
            iteration, state, future, record, failure = pending.popleft()
            if future is not None:
                try:
                    record = future.result()
                except BrokenProcessPool:
                    raise
                except errors as error:
                    failure = f'Simulation {iteration + 1}: {error}'
            yield iteration + 1, record, failure, state
            if next_iteration < stop:
                enqueue()
    except BrokenProcessPool as error:
        raise RuntimeError('Parallel workers stopped unexpectedly. Completed simulations are saved. '
                           'Turn off parallel simulations and run again to continue.') from error
    finally:
        for _, _, future, _, _ in pending:
            if future is not None:
                future.cancel()
