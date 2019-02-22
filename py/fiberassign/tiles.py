# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.tiles
=====================

Functions for loading the tiles

"""
from __future__ import absolute_import, division, print_function

import numpy as np

import desimodel

from ._internal import Tiles


def load_tiles(tiles_file=None, select=None):
    """Load tiles from a file.

    Load tile data either from the specified file or from the default provided
    by desimodel.  Optionally select a subset of tile IDs.

    Args:
        tiles_file (str):  Optional path to a FITS format footprint file.
        select (list):  List of tile IDs to keep when reading the file.

    Returns:
        (Tiles):  A Tiles object.

    """
    # Read in the tile information
    if tiles_file is None:
        tiles_file = desimodel.io.findfile('footprint/desi-tiles.fits')
    tiles_data = desimodel.io.load_tiles(tilesfile=tiles_file, cache=False)

    keeprows = np.arange(len(tiles_data))
    if select is not None:
        keeprows = np.where(np.isin(tiles_data["TILEID"], select))

    tls = Tiles(tiles_data["TILEID"][keeprows], tiles_data["RA"][keeprows],
                tiles_data["DEC"][keeprows],
                tiles_data["OBSCONDITIONS"][keeprows])

    return tls
