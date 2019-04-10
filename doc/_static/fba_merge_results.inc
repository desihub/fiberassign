.. code-block:: console

    usage: fba_merge_results [-h] --targets TARGETS [TARGETS ...] --dir DIR
                             [--prefix PREFIX] [--split] [--out OUT]
                             [--out_prefix OUT_PREFIX] [--out_split]
                             [--columns COLUMNS] [--skip_raw]

    optional arguments:
      -h, --help            show this help message and exit
      --targets TARGETS [TARGETS ...]
                            Input file with targets of any type. This argument can
                            be specified multiple times (for example if standards
                            / skies / science targets are in different files).
      --sky SKY [SKY ...]   Input file with sky or 'bad sky' targets. This option
                        exists in order to treat main-survey sky target files
                        as valid for other survey types. If you are running a
                        main survey assignment, you can just pass the sky file
                        to the --targets list.

      --dir DIR             Directory containing fiberassign results.
      --prefix PREFIX       Prefix of each file (before the <tile>.fits).
      --split               Results are in tile prefix directories.
      --out OUT             Output directory for the merged files. Default is the
                            directory containing the fiberassign output.
      --out_prefix OUT_PREFIX
                            Prefix of each output file.
      --out_split           Split output into tile prefix directories.
      --columns COLUMNS     Override the column names of target data to be copied
                            from the target files into the fiber assignment files.
                            This should be a comma-separated list.
      --skip_raw            Do not copy the raw fiberassign HDUs to the merged
                            output.
