"""Durable, input-specific Monte Carlo draws, without executable serialization."""

import ast
from dataclasses import fields
import hashlib
from io import BytesIO
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np
import scipy


def _calculation_signature(directory=None):
    """Hash sampling and fitting code independently of stopping and reporting.

    AST serialization ignores line endings, formatting, and source locations.
    The draw's model preparation and try/except body preserve the random draw
    order and fit inputs. Changes to those invalidate saved draws.
    """
    directory = Path(directory) if directory is not None else Path(__file__).parent
    tree = ast.parse((directory / 'epsilon_monte_carlo.py').read_text(encoding='utf-8'))
    parts = []
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            parts.append(node)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            if not (isinstance(node, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == 'SCOPE' for t in node.targets)):
                parts.append(node)
        elif isinstance(node, ast.FunctionDef):
            if node.name == '_run_monte_carlo':
                parts.extend(item for item in node.body if isinstance(item, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == 'rng' for t in item.targets))
                start = next(i for i, item in enumerate(node.body)
                             if isinstance(item, ast.Assign) and any(
                                 isinstance(t, ast.Name) and t.id == 'nominal' for t in item.targets))
                end = next(i for i, item in enumerate(node.body)
                           if isinstance(item, ast.Assign) and any(
                               isinstance(t, ast.Name) and t.id == 'rows' for t in item.targets))
                parts.extend(node.body[start:end + 1])
                loop = next(item for item in node.body if isinstance(item, ast.For))
                parts.append(next(item for item in loop.body if isinstance(item, ast.Try)))
            elif node.name not in ('_endpoint_error', '_summary_precision', 'run_monte_carlo'):
                parts.append(node)
    epsilon_tree = ast.parse((directory / 'epsilon_uncertainty.py').read_text(encoding='utf-8'))
    parts.extend(node for node in epsilon_tree.body
                 if isinstance(node, (ast.FunctionDef, ast.Import, ast.ImportFrom)))
    code = '\n'.join(ast.dump(node, include_attributes=False) for node in parts)
    for name in ('kinetics.py', 'pipeline.py', 'spectra.py'):
        code += ast.dump(ast.parse((directory / name).read_text(encoding='utf-8')),
                         include_attributes=False)
    return hashlib.sha256(code.encode('utf-8')).hexdigest()


def _fingerprint(data, reactant, product, error_metric, seed, version, code_hashes):
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

    add((f'AutoQY Monte Carlo checkpoint {version}', sys.version, np.__version__, scipy.__version__,
         error_metric, seed))
    for field in fields(data):
        add(field.name)
        add(getattr(data, field.name))
    for envelope in (reactant, product):
        add(hashlib.sha256(Path(envelope.source_path).read_bytes()).hexdigest())
    for code_hash in code_hashes:
        add(code_hash)
    return digest.hexdigest()


def checkpoint_key(data, reactant, product, error_metric, seed):
    return _fingerprint(data, reactant, product, error_metric, seed, 2,
                        (_calculation_signature(),))


def legacy_checkpoint_keys(data, reactant, product, error_metric, seed):
    from .monte_carlo_legacy import CALCULATION_SIGNATURE, SOURCE_PROFILES

    # Only explicitly verified, numerically identical legacy implementations
    # may migrate. Inputs, seed, Python, NumPy and SciPy must still match exactly.
    if _calculation_signature() != CALCULATION_SIGNATURE:
        return ()
    return tuple(_fingerprint(data, reactant, product, error_metric, seed, 1, profile)
                 for profile in SOURCE_PROFILES)


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
    def __init__(self, directory, key, legacy_keys=()):
        self.path = Path(directory) / f'{key}.sqlite3'
        self.key = key
        self.legacy_keys = tuple(legacy_keys)
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
            if not self.connection.execute('SELECT 1 FROM draws LIMIT 1').fetchone():
                self._import_legacy()
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

    def _import_legacy(self):
        for key in self.legacy_keys:
            if key == self.key or not self.path.with_name(f'{key}.sqlite3').exists():
                continue
            # Lock and verify the original database, then copy all rows atomically.
            # Keep its bytes intact, including every saved draw and RNG state.
            with MonteCarloCheckpoint(self.path.parent, key) as previous:
                with self.connection:
                    self.connection.executemany('INSERT INTO draws VALUES (?, ?, ?, ?)',
                        previous.connection.execute('SELECT iteration, result, error, rng_state '
                                                    'FROM draws ORDER BY iteration'))
                    self.connection.execute('INSERT OR REPLACE INTO metadata VALUES (?, ?)',
                                            ('migrated_from', key))
            break

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
