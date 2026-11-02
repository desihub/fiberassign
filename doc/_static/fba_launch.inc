.. code-block:: console

    usage: fba_launch [-h] --outdir OUTDIR --tileid TILEID [--forcetileid {y,n}]
                      --tilera TILERA --tiledec TILEDEC [--ha HA]
                      --survey {sv1,sv2,sv3,main}
                      --program {DARK,DARK1B,BRIGHT,BRIGHT1B,BACKUP}
                      [--goaltype {DARK,BRIGHT,BACKUP}] [--goaltime GOALTIME]
                      [--sbprof {ELG,BGS,PSF,FLT}] [--mintfrac MINTFRAC]
                      [--rundate RUNDATE] [--mtltime MTLTIME]
                      [--standards_per_petal STANDARDS_PER_PETAL]
                      [--sky_per_petal SKY_PER_PETAL]
                      [--sky_per_slitblock SKY_PER_SLITBLOCK]
                      [--margin-pos MARGIN_POS] [--margin-petal MARGIN_PETAL]
                      [--margin-gfa MARGIN_GFA] [--dr DR] [--gaiadr GAIADR]
                      --dtver DTVER [--pmcorr PMCORR]
                      [--pmtime_utc_str PMTIME_UTC_STR] [--doclean DOCLEAN]
                      [--log-stdout] [--hdr_survey HDR_SURVEY]
                      [--hdr_faprgrm HDR_FAPRGRM] [--steps STEPS]
                      [--nosteps NOSTEPS] [--svntiledir SVNTILEDIR]
                      [--worldreadable] [--lookup_sky_source {ls,gaia}]
                      [--targ_std_only] [--too_tile]
                      [--custom_too_file CUSTOM_TOO_FILE]
                      [--custom_too_development] [--fba_use_fabs FBA_USE_FABS]
                      [--try_too_all] [--skybricks_dir SKYBRICKS_DIR]

    options:
      -h, --help            show this help message and exit
      --outdir OUTDIR       output directory
      --tileid TILEID       output TILEID (e.g., 63142)
      --forcetileid {y,n}   y/n, if y, allows to generate tileids which already
                            exist in svn (default=n)
      --tilera TILERA       tile centre ra
      --tiledec TILEDEC     tile centre dec
      --ha HA               Design for the given Hour Angle in degrees.
      --survey {sv1,sv2,sv3,main}
                            survey
      --program {DARK,DARK1B,BRIGHT,BRIGHT1B,BACKUP}
                            program
      --goaltype {DARK,BRIGHT,BACKUP}
                            goaltype (default=args.program, without '1B')
      --goaltime GOALTIME   goal effective time (default for main: DARK=1000.0,
                            DARK1B=1000.0, BRIGHT=180.0, BRIGHT1B=180.0,
                            BACKUP=60.0)
      --sbprof {ELG,BGS,PSF,FLT}
                            surface brightness profile to be used for computing
                            time -> efftime relation (default for main: DARK=ELG,
                            DARK1B=ELG, BRIGHT=BGS, BRIGHT1B=BGS, BACKUP=PSF)
      --mintfrac MINTFRAC   minimum exposure time fraction needed for this tile to
                            be considered done (default for main: 0.85
      --rundate RUNDATE     yyyy-mm-ddThh:mm:ss+00:00 rundate for focalplane with
                            UTC timezone formatting (default=current UTC time)
      --mtltime MTLTIME     yyyy-mm-ddThh:mm:ss+00:00 MTL isodate, with UTC
                            timezone formatting (default=latest MTL timestamp for
                            args.program if available, if not current UTC time)
      --standards_per_petal STANDARDS_PER_PETAL
                            required number of standards per petal (default=10)
      --sky_per_petal SKY_PER_PETAL
                            required number of sky targets per petal (default=40)
      --sky_per_slitblock SKY_PER_SLITBLOCK
                            Required number of sky targets per fiber slitblock
                            (default=1)
      --margin-pos, --margin_pos MARGIN_POS
                            Add margin (in mm) around positioner keep-out polygons
                            (default: 0.05)
      --margin-petal, --margin_petal MARGIN_PETAL
                            Add margin (in mm) around petal-boundary keep-out
                            polygons (default: 0.4)
      --margin-gfa, --margin_gfa MARGIN_GFA
                            Add margin (in mm) around GFA keep-out polygons
                            (default: 0.4)
      --dr DR               comma-separated list of legacypipe dr
                            (default=dr9,dr11). If multiple are passed, it must be
                            the same length as `--gaiadr`.
      --gaiadr GAIADR       comma-separated list of gaia dr
                            (default=gaiadr2,gaiadr3). If multiple are passed, it
                            must be the same length as `--dr`. The FIRST value is
                            used for gaia ref epochs for GFAs.
      --dtver DTVER         comma-separated list of desitarget catalogue version.
                            If a single value is passed, it is used for ALL data
                            releases. If multiple are passed, exactly ONE must be
                            passed for EACH legacypipe dr.
      --pmcorr PMCORR       apply proper-motion-correction before fiber
                            assignment? (y/n) (default=n)
      --pmtime_utc_str PMTIME_UTC_STR
                            yyyy-mm-ddThh:mm:ss+00:00, UTC time use to compute new
                            coordinates after applying proper motion since
                            REF_EPOCH (default=current UTC time)
      --doclean DOCLEAN     delete TILEID-{tiles,sky,std,gfa,targ,scnd,too}.fits
                            files (y/n)
      --log-stdout, --log_stdout
                            log to stdout instead of redirecting to a file
      --hdr_survey HDR_SURVEY
                            SURVEY keyword in fiberassign-TILEID.fits.gz file
                            (default=args.survey); FAFLAVOR will be
                            {args.hdr_survey}{args.hdr_faprgrm}
      --hdr_faprgrm HDR_FAPRGRM
                            FAPRGRM keyword in fiberassign-TILEID.fits.gz file
                            (default=args.program.lower()); FAFLAVOR will be
                            {args.hdr_survey}{args.hdr_faprgrm}
      --steps STEPS         comma-separated list of steps to accomplish, amongst
                            ['tiles', 'sky', 'gfa', 'targ', 'scnd', 'too', 'fa',
                            'zip', 'move', 'qa'] (defaults to
                            tiles,sky,gfa,targ,scnd,too,fa,zip,move,qa if
                            args.nosteps=None)
      --nosteps NOSTEPS     comma-separated list of steps *not* to accomplish,
                            amongst ['tiles', 'sky', 'gfa', 'targ', 'scnd', 'too',
                            'fa', 'zip', 'move', 'qa'] (default=None); cannot be
                            used if some args.steps values are also provided
      --svntiledir SVNTILEDIR
                            SVN tile directory; fiberassign-* files created should
                            be copied here.
      --worldreadable       Permissions should be world readable.
      --lookup_sky_source {ls,gaia}
                            'ls' or 'gaia' (default='gaia' if main/backup, 'ls'
                            otherwise)
      --targ_std_only       only select standard stars for TILEID-targ.fits (e.g.,
                            for ToO-dedicated tiles); can only be used with
                            args.survey='main'
      --too_tile            do not cut on TOO_TYPE, hence keeping both
                            TOO_TYPE="FIBER" *and* TOO_TYPE="TILE" for ToO
                            dedicated tiles,
      --custom_too_file CUSTOM_TOO_FILE
                            full path to a custom ToO file, or a comma-separated
                            list of custom ToO files, for tertiary program or
                            development work, which overrides the official one
                            (default=None)
      --custom_too_development
                            is this for development? (allows args.custom_too_file
                            to be outside of $DESI_SURVEYOPS)
      --fba_use_fabs FBA_USE_FABS
                            value to determine the cpp behavior, see PR470
                            (default: based on the rundate)
      --try_too_all         if true, look for a file named ToO-fiber-all.ecsv
                            rather than the default ToO-fiber.ecsv. If not found,
                            fallback to ToO-fiber.ecsv. Defaults to False.
      --skybricks_dir SKYBRICKS_DIR
                            if passed, this will override the environment variable
                            $SKYBRICKS_DIR, which is the location used for
                            skybricks lookup for stuck fibers.
