"""Profile imports, app construction, and first HTTP requests without a browser.

Run from the repository: python tools/profile_gui_startup.py --output startup.json
Use a fresh Python process for every sample; compare the same GUI and environment.
This does not measure browser launch, rendering, or the first analysis operation.
"""

from argparse import ArgumentParser
import cProfile
import importlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import pstats
import sys
import time


def main(argv=None):
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--gui', choices=('analysis', 'smoother', 'power'), default='analysis')
    parser.add_argument('--output', type=Path, default=Path('startup.json'))
    args = parser.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    profiler = cProfile.Profile()
    profiler.enable()
    start = time.perf_counter()
    module = importlib.import_module(f'autoqy_core.tools.{args.gui}_gui')
    imported = time.perf_counter()
    app = module.create_app()
    created = time.perf_counter()
    client = app.server.test_client()
    responses = {path: client.get(path).status_code for path in
                 ('/', '/_dash-layout', '/_dash-dependencies')}
    ready = time.perf_counter()
    profiler.disable()
    result = {
        'gui': args.gui,
        'import_seconds': imported - start,
        'create_app_seconds': created - imported,
        'first_requests_seconds': ready - created,
        'total_seconds': ready - start,
        'responses': responses,
        'matplotlib_pyplot_loaded': 'matplotlib.pyplot' in sys.modules,
        'scipy_loaded': 'scipy' in sys.modules,
        'module_count': len(sys.modules),
        'python': sys.version,
        'platform': platform.platform(),
        'dependencies': {name: version(name) for name in
                         ('numpy', 'scipy', 'pandas', 'matplotlib', 'dash', 'plotly')},
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    with args.output.with_suffix('.profile.txt').open('w', encoding='utf-8') as stream:
        pstats.Stats(profiler, stream=stream).sort_stats('cumulative').print_stats(80)
    print(json.dumps(result, indent=2))
    return 0 if all(status == 200 for status in responses.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
