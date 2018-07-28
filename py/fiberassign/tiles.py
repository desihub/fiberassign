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


def load_tiles(hw, tiles_file=None, select=None):
    # Read in the tile information
    if tiles_file is None:
        tiles_file = desimodel.io.findfile('footprint/desi-tiles.fits')
    tiles_data = desimodel.io.load_tiles(tilesfile=tiles_file, cache=False)

    keeprows = np.arange(len(tiles_data))
    if select is not None:
        keeprows = np.where(tiles_data["TILEID"] in select)

    tls = Tiles(hw, tiles_data["TILEID"][keeprows], tiles_data["RA"][keeprows],
                tiles_data["DEC"][keeprows],
                tiles_data["OBSCONDITIONS"][keeprows])

    return tls
