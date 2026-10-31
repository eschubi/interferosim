"""Minimal hello-world example for the Interferosim package."""

from __future__ import annotations

import argparse
from typing import Sequence


def hello_world(name: str = "world") -> str:
    """Return a friendly greeting for the given name."""
    return f"Hello, {name}!"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="interferosim-hello",
        description="Print a simple hello-world message.",
    )
    parser.add_argument(
        "name",
        nargs="?",
        default="world",
        help="name to greet (default: world)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    print(hello_world(args.name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
