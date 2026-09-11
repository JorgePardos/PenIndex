#!/usr/bin/env python3
"""Backward-compatible entry point.

The implementation has moved to the penindex/ package (see penindex/core.py,
penindex/cli.py, and the project README for the new layout). This file is
kept so `python hbond_analyzer.py -c config.yaml` keeps working unchanged.
"""

from penindex.cli import main

if __name__ == "__main__":
    main()
