.. code-block:: console

    usage: fba_fulltest [-h] --targets TARGETS [TARGETS ...] [--sky SKY [SKY ...]]
                        [--gfafile GFAFILE] [--footprint FOOTPRINT]
                        [--tiles TILES] [--out OUT] [--plotpetals PLOTPETALS]
                        [--positioners POSITIONERS] [--status STATUS]
                        [--rundate RUNDATE]
                        [--standards_per_petal STANDARDS_PER_PETAL]
                        [--sky_per_petal SKY_PER_PETAL]

    optional arguments:
      -h, --help            show this help message and exit
      --targets TARGETS [TARGETS ...]
                            Input file with targets of any type. This argument can
                            be specified multiple times (for example if standards
                            / skies / science targets are in different files). By
                            default, the '--mask_column' (default
                            DESI_TARGET)column and bitfield values defined in
                            desitarget are used to determine the type of each
                            target. Each filename may be optionally followed by
                            comma and then one of the strings 'science',
                            'standard', 'sky' or 'safe' to force all targets in
                            that file to be treated as a fixed target type.
      --sky SKY [SKY ...]   Input file with sky or 'bad sky' targets. This option
                            exists in order to treat main-survey sky target files
                            as valid for other survey types. If you are running a
                            main survey assignment, you can just pass the sky file
                            to the --targets list.
      --gfafile GFAFILE     Optional GFA targets FITS file
      --footprint FOOTPRINT
                            Optional FITS file defining the footprint. If not
                            specified, the default footprint from desimodel is
                            used.
      --tiles TILES         Optional text file containing a subset of the tile IDs
                            to use in the footprint, one ID per line. Default uses
                            all tiles in the footprint.
      --out OUT             Top-level output directory.
      --plotpetals PLOTPETALS
                            Comma-separated list of petals to plot (default is all
                            petals)
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
      --standards_per_petal STANDARDS_PER_PETAL
                            Required number of standards per petal
      --sky_per_petal SKY_PER_PETAL
                            Required number of sky targets per petal
