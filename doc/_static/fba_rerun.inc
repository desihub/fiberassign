.. code-block:: console

    usage: fba_rerun [-h] --infiberassign INFIBERASSIGN --outdir OUTDIR
                     [--doclean DOCLEAN] [--dtver DTVER] [--steps STEPS]
                     [--nosteps NOSTEPS] [--no_dvsro] [--force_try_too_all]
                     [--dont_override_skybricks]

    Code to rerun a tile design with a call to fba_launch.
    - SV3: will fail, as fiberassign/master is not backward-compatible for SV3; one should use fiberassign/fba_rerun_io.py functions;
    - Main: for 2021 pre-shutdown tiles designed with dtver=1.0.0, use --dtver 1.1.1 (see [desi-survey 2989] from 2021, Sept. 8th).

    options:
      -h, --help            show this help message and exit
      --infiberassign INFIBERASSIGN
                            full path to fiberassign-TILEID.fits.gz file we want to re-run
      --outdir OUTDIR       output directory
      --doclean DOCLEAN     delete TILEID-{tiles,sky,std,gfa,targ,scnd,too}.fits files? (y/n; default=n); overrides the one in fiberassign-TILEID.fits.gz
      --dtver DTVER         desitarget catalogue version; if not None, overrides the one in fiberassign-TILEID.fits.gz
      --steps STEPS         comma-separated list of steps to accomplish, amongst ['tiles', 'sky', 'gfa', 'targ', 'scnd', 'too', 'fa', 'zip', 'move', 'qa'] (defaults to tiles,sky,gfa,targ,scnd,too,fa,zip,move,qa if args.nosteps=None)
      --nosteps NOSTEPS     comma-separated list of steps *not* to accomplish, amongst ['tiles', 'sky', 'gfa', 'targ', 'scnd', 'too', 'fa', 'zip', 'move', 'qa'] (default=None); cannot be used if some args.steps values are also provided
      --no_dvsro            do not replace '/global/cfs/desi/' by '/dvs_ro/cfs/desi/' in $DESI_TARGET, $SKYBRICKS_DIR, $SKYHEALPIXS_DIR
      --force_try_too_all   Override the value of try_too_all in the infiberassign header and use try_too_all regardless.
      --dont_override_skybricks
                            If passed, do not automatically override the $SKYBRICKS_DIR environment variable. Otherwise, try override with the value in the fiberassign header.
