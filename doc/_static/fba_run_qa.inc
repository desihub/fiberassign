.. code-block:: console

    usage: fba_run_qa [-h] --dir DIR [--prefix PREFIX] [--split] [--out OUT]
                      [--footprint FOOTPRINT] [--tiles TILES] [--rundate RUNDATE]

    options:
      -h, --help            show this help message and exit
      --dir DIR             Directory containing fiberassign outputs.
      --prefix PREFIX       Prefix of each file (before the <tile>.fits).
      --split               Results are in tile prefix directories.
      --out OUT             Output file for QA. Default is a file in directory
                            containing the fiberassign output.
      --footprint FOOTPRINT
                            Optional FITS file defining the footprint. If not
                            specified, the default footprint from desimodel is
                            used.
      --tiles TILES         Optional text file containing a subset of the tile IDs
                            to use in the footprint, one ID per line. Default uses
                            all tiles in the footprint.
      --rundate RUNDATE     Optional date to simulate for this run of fiber
                            assignment, used to load the correct focalplane
                            properties and state from desimodel. Default uses the
                            current date. Format is YYYY-MM-DDTHH:mm:ss in UTC
                            time.
