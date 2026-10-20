"""GUI-independent AutoQY calculations, imported only when requested.

Keep the public API without loading numerical solvers or Matplotlib merely
to launch a GUI or import a lightweight utility.
"""

from importlib import import_module

_EXPORT_MODULES = {
    'AnalysisConfig': 'config',
    'ConfigError': 'config',
    'load_config': 'config',
    'validate_config': 'config',
    'AnalysisInput': 'pipeline',
    'AnalysisResult': 'pipeline',
    'run_analysis_pipeline': 'pipeline',
    'run_concentration_analysis': 'pipeline',
    'format_value_uncertainty': 'output',
    'write_detailed_data': 'output',
    'write_results': 'output',
    'write_figure': 'plotting',
    'PowerBaseline': 'power',
    'PowerResult': 'power',
    'PowerTraceResult': 'power',
    'baseline_power_trace': 'power',
    'combine_power_traces': 'power',
    'process_power_configuration': 'power',
    'run_power_analysis': 'power',
    'RunOutput': 'runner',
    'run_analysis': 'runner',
}

__all__ = list(_EXPORT_MODULES)


def __getattr__(name):
    module = _EXPORT_MODULES.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f".{module}", __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
