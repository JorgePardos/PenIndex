"""Exceptions raised by penindex for expected, user-facing failure modes
(bad config, unreadable structure file, empty selections, etc.).

Only cli.py is expected to catch these and turn them into a clean exit
message; anything else (a genuine bug) should propagate as a normal
traceback rather than being swallowed.
"""


class PenIndexError(Exception):
    """Base class for all expected PenIndex failures."""


class ConfigError(PenIndexError):
    """The YAML configuration file is missing or malformed."""


class StructureLoadError(PenIndexError):
    """The topology/trajectory/QM output could not be loaded."""
