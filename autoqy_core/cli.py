"""Terminal interface for AutoQY Core."""

from argparse import ArgumentParser
import sys

from .config import load_config
from .output import format_value_uncertainty
from .power import run_power_analysis
from .runner import run_analysis


def main(argv=None):
    parser = ArgumentParser(prog="autoqy-core")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="validate a configuration without running")
    validate.add_argument("config")
    run = commands.add_parser("run", help="run an analysis")
    run.add_argument("config")
    run.add_argument("--output-directory")
    power = commands.add_parser("power", help="process optical power-monitor traces")
    power.add_argument("config")
    power.add_argument("--output")
    power_gui = commands.add_parser("power-gui", help="open the power-treatment GUI")
    power_gui.add_argument("--host", default="127.0.0.1")
    power_gui.add_argument("--port", default=8050, type=int)
    power_gui.add_argument("--no-open", action="store_true")
    smoother_gui = commands.add_parser("smoother-gui", help="open the spectral treatment GUI")
    smoother_gui.add_argument("--host", default="127.0.0.1")
    smoother_gui.add_argument("--port", default=8051, type=int)
    smoother_gui.add_argument("--no-open", action="store_true")
    analysis_gui = commands.add_parser("analysis-gui", help="open the analysis builder GUI")
    analysis_gui.add_argument("--host", default="127.0.0.1")
    analysis_gui.add_argument("--port", default=8052, type=int)
    analysis_gui.add_argument("--no-open", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command in {"power-gui", "smoother-gui", "analysis-gui"}:
            if args.command == "power-gui":
                from .tools.power_gui import run_server
            elif args.command == "smoother-gui":
                from .tools.smoother_gui import run_server
            else:
                from .tools.analysis_gui import run_server
            run_server(args.host, args.port, not args.no_open)
            return 0
        if args.command == "power":
            result = run_power_analysis(args.config, args.output)
            formatted = format_value_uncertainty(result.power_mw, result.power_error_mw)
            print(f"Power: {formatted[0]} +/- {formatted[1]} mW")
            return 0
        config = load_config(args.config)
        if args.command == "validate":
            print(f"Valid: {config.source}")
            return 0
        output = run_analysis(config, args.output_directory)
        values = output.result.yield_fit.values * 100
        errors = output.result.yield_errors * 100
        forward = format_value_uncertainty(values[0], errors[0], two_digit_threshold=2)
        backward = format_value_uncertainty(values[1], errors[1], two_digit_threshold=2)
        print(f"Phi_R->P: {forward[0]} +/- {forward[1]}%")
        print(f"Phi_P->R: {backward[0]} +/- {backward[1]}%")
        for path in output.files:
            print(path)
        return 0
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
