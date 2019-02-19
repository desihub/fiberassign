.. code-block:: console

    usage: fba_run [-h] --targets TARGETS [TARGETS ...] [--gfafile GFAFILE]
                   [--footprint FOOTPRINT] [--tiles TILES]
                   [--positioners POSITIONERS] [--status STATUS]
                   [--rundate RUNDATE] [--dir DIR] [--prefix PREFIX] [--split]
                   [--standards_per_petal STANDARDS_PER_PETAL]
                   [--sky_per_petal SKY_PER_PETAL] [--write_all_targets]
                   [--overwrite] [--sciencemask SCIENCEMASK] [--stdmask STDMASK]
                   [--skymask SKYMASK] [--safemask SAFEMASK]
                   [--excludemask EXCLUDEMASK]

    optional arguments:
      -h, --help            show this help message and exit
      --targets TARGETS [TARGETS ...]
                            Input file with targets of any type. This argument can
                            be specified multiple times (for example if standards
                            / skies / science targets are in different files). By
                            default, the DESI_TARGET column and bitfield values
                            defined in desitarget are used to determine the type
                            of each target. Each filename may be optionally
                            followed by comma separated values that will be passed
                            to the arguments of the load_target_file() function.
                            See docs for details.
      --gfafile GFAFILE     Optional GFA targets FITS file
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
                            Default uses the current date. Format is YYYY-MM-
                            DDTHH:mm:ss in UTC time.
      --dir DIR             Output directory.
      --prefix PREFIX       Prefix of each file (before the <tile>.fits).
      --split               Split output into tile prefix directories.
      --standards_per_petal STANDARDS_PER_PETAL
                            Required number of standards per petal
      --sky_per_petal SKY_PER_PETAL
                            Required number of sky targets per petal
      --write_all_targets   When writing target properties, write data for all
                            available targets, not just those which are assigned.
                            This is convenient, but increases the write time and
                            the file size.
      --overwrite           Overwrite any pre-existing output files
      --sciencemask SCIENCEMASK
                            Default DESI_TARGET mask to use for science targets
      --stdmask STDMASK     Default DESI_TARGET mask to use for stdstar targets
      --skymask SKYMASK     Default DESI_TARGET mask to use for sky targets
      --safemask SAFEMASK   Default DESI_TARGET mask to use for safe backup
                            targets
      --excludemask EXCLUDEMASK
                            Default DESI_TARGET mask to exclude from any
                            assignments
