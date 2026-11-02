"""Copy the figures that the top-level README.md embeds from ../figures (git-ignored, rendered
by the scripts here) to ../assets (tracked). Run it after re-rendering a figure that the README
shows.

    uv run python readme_assets.py
"""

import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGURES = HERE.parent / "figures"
ASSETS = HERE.parent / "assets"

# the figures README.md embeds, in the order they appear there
NAMES = ["hos_fig01", "hos_fig05", "soh_fig03", "hos_fig03", "hos_fig08", "hos_fig09", "hso_fig01"]

ASSETS.mkdir(exist_ok=True)
total = 0
for name in NAMES:
    src, dst = FIGURES / f"{name}.png", ASSETS / f"{name}.png"
    shutil.copyfile(src, dst)
    total += dst.stat().st_size
    print(f"{dst.relative_to(HERE.parent)}  {dst.stat().st_size / 1e3:.0f} kB")
print(f"{len(NAMES)} files, {total / 1e6:.2f} MB")
