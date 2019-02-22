.. code-block:: console

    usage: fba_run_qa [-h] --dir DIR [--prefix PREFIX] [--split] [--out OUT]
                      [--footprint FOOTPRINT] [--tiles TILES]
                      [--positioners POSITIONERS] [--status STATUS]
                      [--rundate RUNDATE]

    optional arguments:
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
      --positioners POSITIONERS
                            Optional FITS file describing the fiber positioner
                            locations. Default uses the file from desimodel.
      --status STATUS       Optional fiber status file in astropy ECSV format.
                            Default treats all fibers as good.
      --rundate RUNDATE     Optional date to simulate for this run of fiber
                            assignment, used with the fiber status file to
                            determine which fibers currently have problems.
                            Default uses the current date. Format is YYYY-MM-DD or
                            YYYY-MM-DDTHH:mm:ss in UTC time.
