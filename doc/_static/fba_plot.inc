.. code-block:: console

    usage: fba_plot [-h] [--out_dir OUT_DIR] [--petals PETALS] [--simple]
                    [--serial] files [files ...]

    positional arguments:
      files

    options:
      -h, --help         show this help message and exit
      --out_dir OUT_DIR  Output directory for all plots.
      --petals PETALS    Comma-separated list of petals to plot (default is all
                         petals)
      --simple           Plot simple lines for the positioner arms and no
                         exclusion zones. This is less informative, but faster.
      --serial           Disable the use of multiprocessing. Needed by some
                         unit tests to avoid issues with matplotlib.
