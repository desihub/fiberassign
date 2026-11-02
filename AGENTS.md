# AGENTS.md

This file provides guidance to AI coding agents when working with code in this repository. Claude Code reads it via the `@AGENTS.md` import in `CLAUDE.md`.

## What This Is

`fiberassign` is the DESI fiber assignment tool. It plans which of DESI's 5000 fiber positioners (across 10 petals) observe which astronomical targets during each telescope pointing (~10,000 pointings total, targeting ~50 million objects). The core algorithm is implemented in C++ for performance, exposed to Python via pybind11. Rendered documentation is at https://fiberassign.readthedocs.io/en/latest/ (built from `doc/`).

## Commands

**Build the C++ extension in-place (required before importing or testing):**
```bash
CC=g++ CXX=g++ CFLAGS="-fPIC" LDSHARED="g++ -shared" python setup.py build_ext --inplace
```

**Run all tests:**
```bash
pytest
```

**Run a single test file:**
```bash
pytest py/fiberassign/test/test_assign.py
```

**Build documentation:**
```bash
cd doc && make html
```

## Architecture

The package is a C++/Python hybrid with a clear separation:

**C++ core (`src/`)** — performance-critical algorithms compiled into the `_internal` Python extension via pybind11 (bundled at `src/pybind11/`). The six compiled source files are:
- `assign.cpp` — main assignment algorithm and collision detection
- `hardware.cpp` — positioner geometry, coordinate transforms, exclusion zones
- `targets.cpp` — target data structures and management
- `utils.cpp` — logging, timing, utility functions
- `tiles.cpp` — tile/pointing data management
- `_pyfiberassign.cpp` — all pybind11 bindings (the bridge layer)

`kdTree.h` provides header-only KD-tree spatial indexing used for nearest-neighbor target lookups.

**Python layer (`py/fiberassign/`)** — high-level workflow, FITS I/O, and scripting:
- `assign.py` — top-level assignment API and FITS output
- `hardware.py` — loads hardware configuration from desimodel
- `targets.py` — reads target catalogs and classifies by type (QSO, LRG, ELG, BGS, MWS, sky, standards)
- `fba_launch_io.py` — the main survey launch I/O (largest module at ~4000 lines)
- `fba_patch_io.py` — I/O for patching existing fiberassign outputs (`fba_patch`)
- `fba_rerun_io.py` — I/O for rerunning assignment on existing tiles (`fba_rerun*`)
- `qa.py` / `vis.py` — quality assurance metrics and visualization
- `gfa.py`, `stucksky.py`, `tiles.py`, `fba_tertiary_io.py`, `scripts/`: GFA camera I/O, stuck-sky handling, tile helpers, tertiary-catalog I/O, and the implementations behind the `bin/` entry points

**Scripts (`bin/`)** — 22 entry-point scripts. The main workflow is:
1. `fba_run` — run assignment, producing per-tile FITS files
2. `fba_merge_results` — merge raw output with original target catalogs
3. `fba_run_qa` — compute per-tile QA metrics (writes `qa.json`)
4. `fba_plot` — generate PDF plots of tile assignments

## Key Dependencies

The package requires the DESI ecosystem:
- **desimodel** — instrument configuration; must set `$DESIMODEL` env var pointing to a data directory. The full data comes via svn from `https://desi.lbl.gov/svn/code/desimodel/`, but for running the tests the small test branch is sufficient (this is what CI does, `.github/workflows/test.yml:70`): `svn export https://desi.lbl.gov/svn/code/desimodel/branches/test-0.19/data $DESIMODEL/data`
- **desitarget** — target classification bitmasks
- **desimeter** — metrology/positioning calculations
- **desiutil** — DESI utility functions; `desiutil.setup.get_version` is called at build time

Standard scientific stack: numpy, scipy, astropy, fitsio, matplotlib, healpy.

## Assignment Algorithm (`assign_unused`, `src/assign.cpp:369`)

The core per-tile assignment is **target-centric**, not fiber-centric. A common misconception is that fibers are sorted by the priority of the best target they can reach; that is wrong. The actual loop structure:

1. Targets available for the tile are sorted by priority (highest to lowest).
2. For each target in that order, the available fibers that can physically reach it are iterated closest-to-furthest.
3. The first fiber that passes collision checks (`ok_to_assign`) is assigned.

Fibers have no global ordering — they are only ordered *within the context of a specific target*. `ok_to_assign` rejects a fiber if the positioner is stuck/broken, if moving it would collide with a stuck/broken neighbor, if it would collide with a currently-assigned neighbor, or if it would hit a GFA/petal edge. The "nearest possible fiber" means the nearest one that clears all these checks.

## Calibration Force Algorithm (`assign_force`, `src/assign.cpp:763`)

`assign_force` runs after `assign_unused` to guarantee minimum calibration target counts (sky, standards) per petal and slitblock. It is willing to displace science targets to meet quotas.

Per-tile loop:
1. Collects all fibers currently assigned to pure-science targets (dual science+standard targets are never bumped).
2. Sorts them by **inverse priority** (lowest priority first) — least valuable science targets get bumped first.
3. Calls `tile_available` restricted to those locations to find reachable calibration targets.
4. For each science target (lowest priority first), skips if the petal/slitblock already meets its quota.
5. Finds the highest-priority available calibration target at that location that passes `ok_to_assign`.
6. Calls `reassign_science_target(t+1, tstop, ...)` to attempt relocation of the bumped science target to a **future tile** — while the target is still technically assigned (before `unassign_tileloc`).
7. Unconditionally unassigns the science target; assigns the calibration target; if a future slot was found, places the science target there.

Notable behaviors:
- No early exit once all petals meet their quota — the loop keeps iterating but `petal_count_max` prevents over-assigning.
- Displaced science targets can only go to future tiles (`t+1` onward), so targets bumped from the last tile are permanently lost.
- Tile processing order matters: earlier tiles have more future tiles as fallback for displaced science targets.

## Rundate-Dependent Behavior

Several behaviors are keyed off the `rundate` (an ISO UTC string passed through the launch workflow) via `py/fiberassign/data/cutoff-dates.yaml`. Helper functions in `py/fiberassign/utils.py` (`get_date_cutoff`, `get_obsdate`, `get_fba_use_fabs`, `get_whether_to_reorder_mtl`, `get_whether_to_use_np_concatenate`) compare the rundate against cutoff dates and return era-appropriate values, so a rerun of an old date reproduces historical (sometimes intentionally buggy) behavior. When adding a behavior switch, add a key there rather than hardcoding a date.

## Key Conventions

- Packaging metadata (name, dependencies, package/data discovery, pytest config) lives entirely in `pyproject.toml`'s `[project]`/`[tool.setuptools]`/`[tool.pytest.ini_options]` tables. There is no `setup.cfg`. `setup.py` only retains what can't be expressed declaratively: the pybind11/C++ `Extension` build (`BuildExt`, compiler-flag probing), the `RealClean` command, and the `src/_version.h` stamping step.
- Build artifacts land in `build/` (ignored by git; removed by `setup.py clean`) and the compiled extension `_internal*.so` appears inside the `py/fiberassign/` package directory when built `--inplace`.
- `old/` holds legacy code and, like `bin`, `doc`, `etc`, and `src`, is excluded from pytest collection via `norecursedirs` in `pyproject.toml` — do not treat it as current behavior.
- The test suite writes its outputs under `test_fiberassign_output/` at the repo root (untracked and not gitignored) — do not commit it.
- The C++ extension writes a `src/_version.h` file during `setup.py` execution — this is auto-generated and should not be edited.
- CI (`.github/workflows/test.yml`) has three jobs: unit tests on Ubuntu with Python 3.10–3.13, each pinned to a numpy/astropy range; a coverage job (`pytest --cov`, plus lcov over `build/temp*/src/` for the C++ code); and a docs job (`sphinx-build --keep-going -b html doc doc/_build/html`). The desimodel test data is fetched from DESI's svn server during CI.
- Tests in `py/fiberassign/test/` depend on `$DESIMODEL` being set; many tests are integration-level and require the full DESI data stack.
