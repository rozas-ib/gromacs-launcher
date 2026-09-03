#!/usr/bin/env python3

import argparse

from launcher.analysis import run_analysis


def cli_main():
    parser = argparse.ArgumentParser(description="Run post-simulation analyses for configured systems")
    parser.add_argument("config", nargs="?", default="config.toml")
    args = parser.parse_args()
    try:
        paths = run_analysis(args.config)
    except (KeyError, OSError, ValueError) as exc:
        parser.error(str(exc))
    for path in paths:
        print(f"Analysis report written to {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(cli_main())
