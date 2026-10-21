"""Durable, input-specific Monte Carlo draws, without executable serialization."""

from dataclasses import fields
import hashlib
from io import BytesIO
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np
import scipy


def checkpoint_key(data, reactant, product, error_metric, seed):
    digest = hashlib.sha256()

    def add(value):
        if isinstance(value, np.ndarray):
            array = np.ascontiguousarray(value)
            add((array.dtype.str, array.shape))
            digest.update(array.tobytes())
        elif isinstance(value, (tuple, list)):
            digest.update(b'[')
            for item in value:
                add(item)
            digest.update(b']')
        else:
            encoded = repr(value).encode('utf-8')
            digest.update(str(len(encoded)).encode('ascii') + b':' + encoded)

    add(('AutoQY Monte Carlo checkpoint 1', sys.version, np.__version__, scipy.__version__,
         error_metric, seed))
    for field in fields(data):
        add(field.name)
        add(getattr(data, field.name))
    for envelope in (reactant, product):
        add(hashlib.sha256(Path(envelope.source_path).read_bytes()).hexdigest())
    # Keep saved draws separate when numerical or sampling implementations change.
    for name in ('epsilon_monte_carlo.py', 'kinetics.py', 'pipeline.py', 'spectra.py',
                 'epsilon_uncertainty.py', 'monte_carlo_checkpoint.py'):
        add(hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest())
    return digest.hexdigest()


def _pack(record):
    if record is None:
        return None
    stream = BytesIO()
    np.savez_compressed(stream, **{f'r{i}': np.asarray(value) for i, value in enumerate(record)})
    return stream.getvalue()


def _unpack(payload):
    with np.load(BytesIO(payload), allow_pickle=False) as arrays:
        return (arrays['r0'], arrays['r1'], tuple(bool(x) for x in arrays['r2']),
                float(arrays['r3']), arrays['r4'], arrays['r5'], arrays['r6'], float(arrays['r7']))


class MonteCarloCheckpoint:
    def __init__(self, directory, key):
        self.path = Path(directory) / f'{key}.sqlite3'
        self.key = key
        self.connection = self.guard = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        try:
            # A separate SQLite transaction holds an OS-managed run lock while
            # the draw database commits each result. Crashes release the lock.
            self.guard = sqlite3.connect(str(self.path) + '.lock', timeout=0)
            self.guard.execute('BEGIN EXCLUSIVE')
            self.connection = sqlite3.connect(self.path)
            self.connection.execute('PRAGMA synchronous=FULL')
            self.connection.execute('CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)')
            self.connection.execute('CREATE TABLE IF NOT EXISTS draws '
                                    '(iteration INTEGER PRIMARY KEY, result BLOB, error TEXT, rng_state TEXT NOT NULL)')
            with self.connection:
                self.connection.execute('INSERT OR IGNORE INTO metadata VALUES (?, ?)', ('fingerprint', self.key))
            key = self.connection.execute('SELECT value FROM metadata WHERE key=?', ('fingerprint',)).fetchone()
            if key != (self.key,):
                raise ValueError('The saved simulations do not match this analysis. Choose a different output folder.')
            return self
        except sqlite3.OperationalError as error:
            self.__exit__(None, None, None)
            if 'locked' in str(error).lower():
                raise ValueError('These Monte Carlo simulations are already running in another analysis. '
                                 'Wait for that run to finish, then try again.') from error
            raise
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *args):
        if self.connection is not None:
            self.connection.close()
        if self.guard is not None:
            self.guard.close()

    def load(self, maximum):
        records, failures, state = [], [], None
        attempted = 0
        try:
            for iteration, payload, error, rng_state in self.connection.execute(
                    'SELECT iteration, result, error, rng_state FROM draws '
                    'WHERE iteration <= ? ORDER BY iteration', (maximum,)):
                attempted += 1
                if iteration != attempted:
                    raise ValueError('Saved simulation numbers are incomplete')
                if error is None:
                    records.append(_unpack(payload))
                else:
                    failures.append(error)
                state = json.loads(rng_state)
        except (ValueError, KeyError, TypeError, sqlite3.DatabaseError) as error:
            raise ValueError('The saved Monte Carlo simulations could not be read. '
                             'Choose a different output folder to start a fresh calculation.') from error
        return records, failures, attempted, state

    def append(self, iteration, record, error, rng_state):
        with self.connection:
            self.connection.execute('INSERT INTO draws VALUES (?, ?, ?, ?)',
                                    (iteration, _pack(record), error, json.dumps(rng_state)))
