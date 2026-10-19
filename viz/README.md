# viz: plotting and visualisation for HoyngDynamo

Python side of the project (dev/DECISIONS.md D13). Reads the HDF5 output of the Julia
package and makes the figures; it is never needed to run a simulation.

Managed by [uv](https://docs.astral.sh/uv/); Python 3.12, matplotlib, numpy, h5py.

```sh
cd viz
uv sync                      # create .venv from uv.lock
uv run python <script>.py    # run a figure script in the environment
```

Layout (to come, build step 6): one script per paper figure named after it
(`hos_fig01.py`, …), writing to `../figures/`, plus a shared style module and a reader
for the HDF5 run files.
