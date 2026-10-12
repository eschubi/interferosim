# src/interferosim/sky/cli.py -- command line tool "interferosim-sky" (interferosim.sky.cli)
"""Show or save the map of a sky model file.

Usage::

    interferosim-sky configs/skies/aas6_demo.json                 # window
    interferosim-sky configs/skies/aas6_demo.json -o demo.png     # file
    interferosim-sky configs/skies/aas6_demo.json --overview -o demo.png
    interferosim-sky --demo -o demo.png                           # built-in demo

Version 0.2.0 (2026-10-12): first implementation.
"""

from __future__ import annotations

import argparse
from typing import List, Optional

from .models import SkyModel


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="interferosim-sky",
                                description="Show the map of an Interferosim sky model (JSON).")
    p.add_argument("file", nargs="?", help="sky model file (sim_sky.json format)")
    p.add_argument("--demo", action="store_true", help="use the built-in demo field")
    p.add_argument("-o", "--out", help="save the map to this file instead of showing it")
    p.add_argument("--overview", action="store_true", help="full field plus zoom")
    p.add_argument("--fov", type=float, default=None, help="field of view [deg]")
    p.add_argument("--npix", type=int, default=512, help="pixels per side (default 512)")
    p.add_argument("--scale", choices=["log", "linear", "asinh"], default="log")
    p.add_argument("--list", action="store_true", help="print the component list")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.demo:
        sky = SkyModel.demo()
    elif args.file:
        sky = SkyModel.from_json(args.file)
    else:
        build_parser().error("give a sky model file or --demo")
        return 2

    if args.list:
        print(f"{sky.name or 'sky model'}: centre {sky.center.label} "
              f"(RA {sky.center.ra_deg}, Dec {sky.center.dec_deg})")
        for c in sky:
            print(f"  {c.kind:6s} {c.name:<20s} dx {c.dx_deg:+9.5f}  dy {c.dy_deg:+9.5f}"
                  f"  {c.flux_jy:10.1f} Jy")
        print(f"  total {sky.total_flux_jy:.1f} Jy")

    import matplotlib
    if args.out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..plotting import plot_sky_map, plot_sky_overview
    if args.overview:
        fig, _ = plot_sky_overview(sky, npix=args.npix, scale=args.scale)
    else:
        fig, _, _ = plot_sky_map(sky, fov_deg=args.fov, npix=args.npix, scale=args.scale)
    if args.out:
        fig.savefig(args.out, dpi=150, bbox_inches="tight")
        print(f"saved {args.out}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
