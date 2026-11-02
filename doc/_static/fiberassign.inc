.. code-block:: console

    usage: fiberassign [-h] --mtl MTL --sky SKY [--stdstar STDSTAR]
                       [--fibstatusfile FIBSTATUSFILE] [--footprint FOOTPRINT]
                       [--positioners POSITIONERS] [--surveytiles SURVEYTILES]
                       [--outdir OUTDIR] [--sciencemask SCIENCEMASK]
                       [--starmask STARMASK] [--stdmask STDMASK]
                       [--skymask SKYMASK] [--safemask SAFEMASK]
                       [--excludemask EXCLUDEMASK] [--rundate RUNDATE]
                       [--gfafile GFAFILE] [--nstarpetal NSTARPETAL]
                       [--nskypetal NSKYPETAL] [--overwrite] [--version]

    options:
      -h, --help            show this help message and exit
      --mtl MTL             input targets (FITS file)
      --sky SKY             input sky positions (FITS file)
      --stdstar STDSTAR     input std stars (FITS file)
      --fibstatusfile FIBSTATUSFILE
                            list of positioners and its status (ECSV or txt file)
      --footprint FOOTPRINT
                            list of tiles defining the footprint (FITS file)
      --positioners POSITIONERS
                            list of positioners on the focal plane (FITS file)
      --surveytiles SURVEYTILES
                            set of tiles to run fiberassign on (text file)
      --outdir OUTDIR       output directory (default = ./)
      --sciencemask SCIENCEMASK
                            Default DESI_TARGET mask to use for science targets
      --starmask STARMASK   DEPRECATED: use --stdmask instead
      --stdmask STDMASK     Default DESI_TARGET mask to use for stdstar targets
      --skymask SKYMASK     Default DESI_TARGET mask to use for sky targets
      --safemask SAFEMASK   Default DESI_TARGET mask to use for safe backup
                            targets
      --excludemask EXCLUDEMASK
                            Default DESI_TARGET mask to exclude from any
                            assignments
      --rundate RUNDATE     run date [YYYY-MM-DD]
      --gfafile GFAFILE     GFA file (FITS tile)
      --nstarpetal NSTARPETAL
                            number of standard stars per petal (default=10)
      --nskypetal NSKYPETAL
                            number of sky fibers per petal (default=40)
      --overwrite           overwrite pre-existing output files
      --version             Print code version and exit
