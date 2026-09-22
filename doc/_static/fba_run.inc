.. code-block:: console

    usage: fba_run [-h] --targets TARGETS [TARGETS ...] [--sky SKY [SKY ...]]
                   [--gfafile GFAFILE] [--footprint FOOTPRINT] [--tiles TILES]
                   [--rundate RUNDATE] [--obsdate OBSDATE] [--ha HA]
                   [--fieldrot FIELDROT] [--fieldrot_corr FIELDROT_CORR]
                   [--dir DIR] [--prefix PREFIX] [--split]
                   [--standards_per_petal STANDARDS_PER_PETAL]
                   [--sky_per_petal SKY_PER_PETAL]
                   [--sky_per_slitblock SKY_PER_SLITBLOCK]
                   [--margin-pos MARGIN_POS] [--margin-petal MARGIN_PETAL]
                   [--margin-gfa MARGIN_GFA] [--write_all_targets] [--overwrite]
                   [--mask_column MASK_COLUMN] [--sciencemask SCIENCEMASK]
                   [--stdmask STDMASK] [--gaia_stdmask GAIA_STDMASK]
                   [--skymask SKYMASK] [--safemask SAFEMASK]
                   [--excludemask EXCLUDEMASK] [--by_tile] [--no_redistribute]
                   [--no_zero_obsremain] [--lookup_sky_source {ls,gaia}]
                   [--fafns_for_stucksky FAFNS_FOR_STUCKSKY]
                   [--fba_use_fabs FBA_USE_FABS]

    options:
      -h, --help            show this help message and exit
      --targets TARGETS [TARGETS ...]
                            Input file with targets of any type. This argument
                            can be specified multiple times (for example if
                            standards / skies / science targets are in different
                            files). By default, the '--mask_column' (default
                            DESI_TARGET) column and bitfield values defined in
                            desitarget are used to determine the type of each
                            target. Each filename may be optionally followed by
                            comma and then one of the strings 'science',
                            'standard', 'sky' or 'safe' to force all targets in
                            that file to be treated as a fixed target type.
      --sky SKY [SKY ...]   Input file with sky or 'bad sky' targets. This option
                            exists in order to treat main-survey sky target files
                            as valid for other survey types. If you are running
                            a main survey assignment, you can just pass the sky
                            file to the --targets list.
      --gfafile GFAFILE     Optional GFA targets FITS file
      --footprint FOOTPRINT
                            Optional FITS file defining the footprint. If not
                            specified, the default footprint from desimodel is
                            used.
      --tiles TILES         Optional text file containing a subset of the tile IDs
                            to use in the footprint, one ID per line.
                            Default uses all tiles in the footprint.
      --rundate RUNDATE     Optional date to simulate for this run of fiber
                            assignment, used to load the correct focalplane
                            properties and state from desimodel. Default uses the
                            current date. Format is YYYY-MM-DDTHH:mm:ss+-zz:zz.
      --obsdate OBSDATE     Plan field rotations for this date (YEARMMDD, or
                            ISO 8601 YEAR-MM-DD with or without time);
                            default=2022-07-01 or rundate+1yr, according to rundate
      --ha HA               Design for the given Hour Angle in degrees.
      --fieldrot FIELDROT   Override obsdate and use this field rotation for
                            all tiles (degrees counter clockwise in CS5)
      --fieldrot_corr FIELDROT_CORR
                            apply correction to computed fieldrot (default:
                            False if args.fieldrot is provided or if
                            args.rundate<rundate_cutoff, True else)
      --dir DIR             Output directory.
      --prefix PREFIX       Prefix of each file (before the <tile>.fits).
      --split               Split output into tile prefix directories.
      --standards_per_petal STANDARDS_PER_PETAL
                            Required number of standards per petal
      --sky_per_petal SKY_PER_PETAL
                            Required number of sky targets per petal
      --sky_per_slitblock SKY_PER_SLITBLOCK
                            Required number of sky targets per fiber slitblock
      --margin-pos MARGIN_POS, --margin_pos MARGIN_POS
                            Add margin (in mm) around positioner keep-out polygons
      --margin-petal MARGIN_PETAL, --margin_petal MARGIN_PETAL
                            Add margin (in mm) around petal-boundary keep-out polygons
      --margin-gfa MARGIN_GFA, --margin_gfa MARGIN_GFA
                            Add margin (in mm) around GFA keep-out polygons
      --write_all_targets   When writing target properties, write data for all
                            available targets, not just those which are assigned.
                            This is convenient, but increases the write time
                            and the file size.
      --overwrite           Overwrite any pre-existing output files
      --mask_column MASK_COLUMN
                            FITS column to use for applying target masks
      --sciencemask SCIENCEMASK
                            Default DESI_TARGET mask to use for science targets
      --stdmask STDMASK     Default DESI_TARGET mask to use for stdstar targets
      --gaia_stdmask GAIA_STDMASK
                            Default MWS_TARGET mask to use for Gaia stdstar targets
                            (for BACKUP program)
      --skymask SKYMASK     Default DESI_TARGET mask to use for sky targets
      --safemask SAFEMASK   Default DESI_TARGET mask to use for safe backup targets
      --excludemask EXCLUDEMASK
                            Default DESI_TARGET mask to exclude from any assignments
      --by_tile             Do assignment one tile at a time. This disables
      redistribution.
      --no_redistribute     Disable redistribution of science targets.
      --no_zero_obsremain   Disable oversubscription of science targets with
                            leftover fibers.
      --lookup_sky_source {ls,gaia}
                            Source for the look-up table for sky positions for
                            stuck fibers: 'ls': uses $SKYBRICKS_DIR;
                            'gaia': uses $SKYHEALPIXS_DIR (default=ls)
      --fafns_for_stucksky FAFNS_FOR_STUCKSKY
                            csv list of fiberassign-TILEID.fits.gz files;
                            if set, fiberassign will use this fiberassign file
                            to know which stuck fibers are usable for sky
      --fba_use_fabs FBA_USE_FABS
                            value to determine the cpp behavior,
                            see PR470 (default: based on the rundate)

