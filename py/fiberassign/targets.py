# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.targets
=====================

Functions for loading the target list

"""
from __future__ import absolute_import, division, print_function

import sys
import numpy as np

import fitsio

# FIXME:  If / when SV bit names diverge from main survey names, we
# should import the SV bitmasks here.
# AR duplicating what s required for sv2;
# AR should definitely be re-written to handle sv3, etc
# AR and duplicating for sv3...

from desitarget.targetmask import desi_mask, mws_mask

from desitarget.cmx.cmx_targetmask import cmx_mask

from desitarget.sv1.sv1_targetmask import desi_mask as sv1_mask
from desitarget.sv1.sv1_targetmask import scnd_mask as sv1_scnd_mask
# AR
from desitarget.sv2.sv2_targetmask import desi_mask as sv2_mask
from desitarget.sv3.sv3_targetmask import desi_mask as sv3_mask

from desitarget.targets import main_cmx_or_sv

from .utils import Logger, Timer
from .hardware import radec2xy, cs52xy

from ._internal import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                        TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                        TARGET_TYPE_SUPPSKY,
                        Target, Targets, TargetsAvailable,
                        LocationsAvailable)

from fiberassign.utils import assert_isoformat_utc, get_date_cutoff

from datetime import datetime
from astropy.time import Time

import desimodel.io

class TargetTagalong(object):
    '''
    This class holds data from the targeting input files that we want
    to propagate to the output fiberassign files, and that are not
    needed by the C++ layer.
    '''
    def __init__(self, columns, outnames={}, aliases={}):
        '''
        Create a new tag-along object.

        Args:
        *columns*: list of strings: the column names that will be saved.
        *outnames*: dict, string to string: mapping from 'columns' to the name
                    the column will be given in the output file; None to omit
                    from the output file.
        *aliases*: dict, string to string: for get_for_ids(), column aliases.
        '''
        self.columns = columns
        self.outnames = outnames
        self.aliases = aliases
        # Internally, we store one tuple for each targeting file read
        # (to avoid manipulating/reformatting the arrays too much),
        # where each tuple starts with the TARGETID of the targets, followed
        # by the data arrays for each column in *columns*.
        self.data = []

    def get_default(self, column):
        '''
        Returns the default value to return for a given *column*, or
        None if not set.
        '''
        return None

    def get_output_name(self, column):
        '''
        Returns the column name to use in the output file for the
        given input column name.
        '''
        return self.outnames.get(column, column)

    def add_data(self, targetids, tabledata, fake={}):
        '''
        Stores data from an input (targeting file) table.

        Arguments:
        *targetids*: numpy array of TARGETID values that must be in the same
                     order as the data arrays in *tabledata*.
        *tabledata*: numpy record-array / table from which this tagalong's
                     *columns* will be read.
        *fake*: dict from string column name to numpy array, containing column
                     data that will be used in place of reading from *tabledata*.
        '''
        tgarrays = [targetids]
        for k in self.columns:
            if k in fake:
                assert(len(fake[k]) == len(targetids))
                tgarrays.append(fake[k])
            else:
                assert(len(tabledata[k]) == len(targetids))
                tgarrays.append(tabledata[k])
        self.data.append(tgarrays)

    def set_data(self, targetids, tabledata):
        '''
        Sets *ALL* rows of the given *tabledata* object to defaults,
        and then fills in values for the given targetids, if they are found.
        '''
        # Set defaults, and grab output arrays
        outarrs = []
        for c in self.columns:
            defval = self.get_default(c)
            outname = self.get_output_name(c)
            if outname is None:
                # We're omitting this column from the output
                outarrs.append(None)
            else:
                outarr = tabledata[outname]
                if defval is not None:
                    outarr[:] = defval
                outarrs.append(outarr)
        # Build output targetid-to-index map
        outmap = dict([(tid,i) for i,tid in enumerate(targetids)])
        # Go through my many data arrays
        for thedata in self.data:
            # TARGETIDs are the first element in the tuple
            tids = thedata[0]
            # Search for output array indices for these targetids
            outinds = np.array([outmap.get(tid, -1) for tid in tids])
            # Keep only the indices of targetids that were found
            ininds = np.flatnonzero(outinds >= 0)
            outinds = outinds[ininds]
            for outarr,inarr in zip(outarrs, thedata[1:]):
                if outarr is None:
                    continue
                outarr[outinds] = inarr[ininds]

    def get_for_ids(self, targetids, names):
        '''
        Fetch arrays for the given columns names and given targetids.
        '''
        # Create output arrays
        outarrs = []
        colinds = []
        for name in names:
            name = self.aliases.get(name, name)
            ic = self.columns.index(name)
            # Look at the data saved for my first dataset to determine
            # the output type.
            dtype = self.data[0][ic+1].dtype
            outarrs.append(np.zeros(len(targetids), dtype))
            colinds.append(ic+1)
        # Build output targetid-to-index map
        outmap = dict([(tid,i) for i,tid in enumerate(targetids)])
        # Go through my many data arrays
        for thedata in self.data:
            tids = thedata[0]
            # Search for output array indices for these targetids
            outinds = np.array([outmap.get(tid, -1) for tid in tids])
            ininds = np.flatnonzero(outinds >= 0)
            outinds = outinds[ininds]
            for outarr,ic in zip(outarrs, colinds):
                outarr[outinds] = thedata[ic][ininds]
        return outarrs

def create_tagalong(plate_radec=True):
    cols = [
        'TARGET_RA',
        'TARGET_DEC',
        'OBSCOND',
        'FA_TARGET',
        ]

    if plate_radec:
        cols.extend([
            'PLATE_RA',
            'PLATE_DEC',
            # 'PLATE_REF_EPOCH',
            ])
        # If PLATE_{RA,DEC} exist in the input target tables, use those
        # when converting RA,DEC to focal-plane coords.
        aliases = {
            'RA':  'PLATE_RA',
            'DEC': 'PLATE_DEC',
            }
    else:
        aliases = {
            'RA':  'TARGET_RA',
            'DEC': 'TARGET_DEC',
            }

    # (OBSCOND doesn't appear in all the fiberassign output HDUs,
    # so we handle it specially)
    return TargetTagalong(cols, outnames={'OBSCOND':None}, aliases=aliases)

def str_to_target_type(input):
    if input == "science":
        return TARGET_TYPE_SCIENCE
    elif input == "sky":
        return TARGET_TYPE_SKY
    elif input == "suppsky":
        return TARGET_TYPE_SUPPSKY
    elif input == "standard":
        return TARGET_TYPE_STANDARD
    elif input == "safe":
        return TARGET_TYPE_SAFE
    else:
        raise ValueError("unknown target type '{}'".format(input))
    return None


def default_main_sciencemask():
    """Returns default mask of bits for science targets in main survey.

    Notes:
        20250221: enable "LGE"
    """
    sciencemask = 0
    sciencemask |= desi_mask["LRG"].mask
    sciencemask |= desi_mask["LGE"].mask
    sciencemask |= desi_mask["ELG"].mask
    sciencemask |= desi_mask["QSO"].mask
    sciencemask |= desi_mask["BGS_ANY"].mask
    sciencemask |= desi_mask["MWS_ANY"].mask
    if "SCND_ANY" in desi_mask.names():
        sciencemask |= desi_mask["SCND_ANY"].mask
    else:
        sciencemask |= desi_mask["SECONDARY_ANY"].mask
    return sciencemask


def default_main_stdmask(rundate=None):
    """Returns default mask of bits for standards in main survey.

    Args:
        rundate (optional, defaults to None): yyyy-mm-ddThh:mm:ss+00:00 rundate
            for focalplane with UTC timezone formatting (string)

    Notes:
        20210930 : we make the default behavior to discard STD_WD;
                    if rundate < 2021-10-01T19:00:00+00:00, we include STD_WD in the stdmask
                    to preserve backward-compatibility.
    """
    log = Logger.get()

    stdmask = 0
    stdmask |= desi_mask["STD_FAINT"].mask
    stdmask |= desi_mask["STD_BRIGHT"].mask
    # AR cutoff rundate for:
    # AR - including or not STD_WD
    rundate_cutoff = get_date_cutoff("rundate", "std_wd")
    rundate_mjd_cutoff = Time(datetime.strptime(rundate_cutoff, "%Y-%m-%dT%H:%M:%S%z")).mjd
    use_wd = False
    if rundate is not None:
        if not assert_isoformat_utc(rundate):
            log.info("provided rundate={} does not follow the expected formatting (yyyy-mm-ddThh:mm:ss+00:00); exiting".format(rundate))
            sys.exit(1)
        rundate_mjd = Time(datetime.strptime(rundate, "%Y-%m-%dT%H:%M:%S%z")).mjd
        if rundate_mjd < rundate_mjd_cutoff:
            use_wd = True
    if use_wd:
        stdmask |= desi_mask["STD_WD"].mask
        log.info("rundate = {} is before rundate_cutoff_std_wd = {}, so STD_WD are counted as TARGET_TYPE_STANDARD".format(rundate, rundate_cutoff))
    else:
        log.info("rundate = {} is after rundate_cutoff_std_wd = {}, so STD_WD are *not* counted as TARGET_TYPE_STANDARD".format(rundate, rundate_cutoff))
    return stdmask

def default_main_gaia_stdmask():
    """Returns default mask of bits for Gaia standards in main survey.
        Purposely not include GAIA_STD_WD, as we plan to remove WD from
        the standards counting requirement.
    """
    gaia_stdmask = 0
    gaia_stdmask |= mws_mask["GAIA_STD_FAINT"].mask
    gaia_stdmask |= mws_mask["GAIA_STD_BRIGHT"].mask
    return gaia_stdmask

def default_main_skymask():
    """Returns default mask of bits for sky targets in main survey.
    """
    skymask = 0
    skymask |= desi_mask["SKY"].mask
    return skymask


def default_main_suppskymask():
    """Returns default mask of bits for suppsky targets in main survey.
    """
    suppskymask = 0
    suppskymask |= desi_mask["SUPP_SKY"].mask
    return suppskymask


def default_main_safemask():
    """Returns default mask of bits for 'safe' targets in main survey.

    Note: these are targets of last resort; they are safe locations where
    we won't saturate the detector, but aren't good for anything else.
    """
    safemask = 0
    safemask |= desi_mask["BAD_SKY"].mask
    return safemask


def default_main_excludemask():
    """Returns default mask of bits for main survey targets to NOT observe.
    """
    excludemask = 0
    # Exclude BRIGHT_OBJECT and IN_BRIGHT_OBJECT, but not NEAR_BRIGHT_OBJECT
    excludemask |= desi_mask.BRIGHT_OBJECT
    excludemask |= desi_mask.IN_BRIGHT_OBJECT
    return excludemask


def default_sv3_sciencemask():
    """Returns default mask of bits for science targets in SV1 survey.
    """
    sciencemask = 0
    sciencemask |= sv3_mask["LRG"].mask
    sciencemask |= sv3_mask["ELG"].mask
    sciencemask |= sv3_mask["QSO"].mask
    sciencemask |= sv3_mask["BGS_ANY"].mask
    sciencemask |= sv3_mask["MWS_ANY"].mask
    sciencemask |= sv3_mask["SCND_ANY"].mask
    return sciencemask

def default_sv3_stdmask():
    """Returns default mask of bits for standards in SV1 survey.
    """
    stdmask = 0
    stdmask |= sv3_mask["STD_FAINT"].mask
    stdmask |= sv3_mask["STD_WD"].mask
    stdmask |= sv3_mask["STD_BRIGHT"].mask
    return stdmask

def default_sv3_skymask():
    """Returns default mask of bits for sky targets in SV1 survey.
    """
    skymask = 0
    skymask |= sv3_mask["SKY"].mask
    return skymask

def default_sv3_suppskymask():
    """Returns default mask of bits for suppsky targets in SV1 survey.
    """
    suppskymask = 0
    suppskymask |= sv3_mask["SUPP_SKY"].mask
    return suppskymask

def default_sv3_safemask():
    """Returns default mask of bits for 'safe' targets in SV1 survey.

    Note: these are targets of last resort; they are safe locations where
    we won't saturate the detector, but aren't good for anything else.
    """
    safemask = 0
    safemask |= sv3_mask["BAD_SKY"].mask
    return safemask

def default_sv3_excludemask():
    """Returns default mask of bits for SV1 survey targets to NOT observe.
    """
    excludemask = 0
    # Exclude BRIGHT_OBJECT and IN_BRIGHT_OBJECT, but not NEAR_BRIGHT_OBJECT
    excludemask |= sv3_mask.BRIGHT_OBJECT
    excludemask |= sv3_mask.IN_BRIGHT_OBJECT
    return excludemask


def default_sv2_sciencemask():
    """Returns default mask of bits for science targets in SV1 survey.
    """
    sciencemask = 0
    sciencemask |= sv2_mask["LRG"].mask
    sciencemask |= sv2_mask["ELG"].mask
    sciencemask |= sv2_mask["QSO"].mask
    sciencemask |= sv2_mask["BGS_ANY"].mask
    sciencemask |= sv2_mask["MWS_ANY"].mask
    sciencemask |= sv2_mask["SCND_ANY"].mask
    return sciencemask

def default_sv2_stdmask():
    """Returns default mask of bits for standards in SV1 survey.
    """
    stdmask = 0
    stdmask |= sv2_mask["STD_FAINT"].mask
    stdmask |= sv2_mask["STD_WD"].mask
    stdmask |= sv2_mask["STD_BRIGHT"].mask
    return stdmask


def default_sv2_skymask():
    """Returns default mask of bits for sky targets in SV1 survey.
    """
    skymask = 0
    skymask |= sv2_mask["SKY"].mask
    return skymask

def default_sv2_suppskymask():
    """Returns default mask of bits for suppsky targets in SV1 survey.
    """
    suppskymask = 0
    suppskymask |= sv2_mask["SUPP_SKY"].mask
    return suppskymask


def default_sv2_safemask():
    """Returns default mask of bits for 'safe' targets in SV1 survey.

    Note: these are targets of last resort; they are safe locations where
    we won't saturate the detector, but aren't good for anything else.
    """
    safemask = 0
    safemask |= sv2_mask["BAD_SKY"].mask
    return safemask


def default_sv2_excludemask():
    """Returns default mask of bits for SV1 survey targets to NOT observe.
    """
    excludemask = 0
    # Exclude BRIGHT_OBJECT and IN_BRIGHT_OBJECT, but not NEAR_BRIGHT_OBJECT
    excludemask |= sv2_mask.BRIGHT_OBJECT
    excludemask |= sv2_mask.IN_BRIGHT_OBJECT
    return excludemask

def default_sv1_sciencemask():
    """Returns default mask of bits for science targets in SV1 survey.
    """
    sciencemask = 0
    sciencemask |= sv1_mask["LRG"].mask
    sciencemask |= sv1_mask["ELG"].mask
    sciencemask |= sv1_mask["QSO"].mask
    sciencemask |= sv1_mask["BGS_ANY"].mask
    sciencemask |= sv1_mask["MWS_ANY"].mask

    if "SCND_ANY" in desi_mask.names():
        sciencemask |= desi_mask["SCND_ANY"].mask
    else:
        sciencemask |= desi_mask["SECONDARY_ANY"].mask

    return sciencemask


def default_sv1_stdmask():
    """Returns default mask of bits for standards in SV1 survey.
    """
    stdmask = 0
    stdmask |= sv1_mask["STD_FAINT"].mask
    stdmask |= sv1_mask["STD_WD"].mask
    stdmask |= sv1_mask["STD_BRIGHT"].mask
    return stdmask


def default_sv1_skymask():
    """Returns default mask of bits for sky targets in SV1 survey.
    """
    skymask = 0
    skymask |= sv1_mask["SKY"].mask
    return skymask


def default_sv1_suppskymask():
    """Returns default mask of bits for suppsky targets in SV1 survey.
    """
    suppskymask = 0
    suppskymask |= sv1_mask["SUPP_SKY"].mask
    return suppskymask


def default_sv1_safemask():
    """Returns default mask of bits for 'safe' targets in SV1 survey.

    Note: these are targets of last resort; they are safe locations where
    we won't saturate the detector, but aren't good for anything else.
    """
    safemask = 0
    safemask |= sv1_mask["BAD_SKY"].mask
    return safemask


def default_sv1_excludemask():
    """Returns default mask of bits for SV1 survey targets to NOT observe.
    """
    excludemask = 0
    # Exclude BRIGHT_OBJECT and IN_BRIGHT_OBJECT, but not NEAR_BRIGHT_OBJECT
    excludemask |= sv1_mask.BRIGHT_OBJECT
    excludemask |= sv1_mask.IN_BRIGHT_OBJECT
    return excludemask


def default_cmx_sciencemask():
    """Returns default mask of bits for science targets in CMX survey.
    """
    sciencemask = 0
    sciencemask |= cmx_mask["STD_GAIA"].mask
    sciencemask |= cmx_mask["SV0_STD_BRIGHT"].mask
    sciencemask |= cmx_mask["STD_TEST"].mask
    sciencemask |= cmx_mask["STD_CALSPEC"].mask
    sciencemask |= cmx_mask["STD_DITHER"].mask
    sciencemask |= cmx_mask["STD_FAINT"].mask
    sciencemask |= cmx_mask["SV0_BGS"].mask
    sciencemask |= cmx_mask["SV0_MWS"].mask
    sciencemask |= cmx_mask["SV0_LRG"].mask
    sciencemask |= cmx_mask["SV0_ELG"].mask
    sciencemask |= cmx_mask["SV0_QSO"].mask
    sciencemask |= cmx_mask["SV0_WD"].mask
    sciencemask |= cmx_mask["BACKUP_BRIGHT"].mask
    sciencemask |= cmx_mask["BACKUP_FAINT"].mask
    sciencemask |= cmx_mask["M31_STD_BRIGHT"].mask
    sciencemask |= cmx_mask["M31_H2PN"].mask
    sciencemask |= cmx_mask["M31_GC"].mask
    sciencemask |= cmx_mask["M31_VAR"].mask
    sciencemask |= cmx_mask["M31_QSO"].mask
    sciencemask |= cmx_mask["M31_BSPL"].mask
    sciencemask |= cmx_mask["M31_M31cen"].mask
    sciencemask |= cmx_mask["M31_M31out"].mask
    sciencemask |= cmx_mask["ORI_STD_BRIGHT"].mask
    sciencemask |= cmx_mask["ORI_QSO"].mask
    sciencemask |= cmx_mask["ORI_ORI"].mask
    sciencemask |= cmx_mask["ORI_HA"].mask
    # NEW bits for SV0- March 2020- desitarget 0.37.0
    sciencemask |= cmx_mask["SV0_QSO_Z5"].mask
    sciencemask |= cmx_mask["SV0_MWS_CLUSTER"].mask
    sciencemask |= cmx_mask["SV0_MWS_CLUSTER_VERYBRIGHT"].mask



    sciencemask |= cmx_mask["ROS_STD_BRIGHT"].mask
    sciencemask |= cmx_mask["ROS_QSO"].mask
    sciencemask |= cmx_mask["ROS_ROSM17"].mask
    sciencemask |= cmx_mask["ROS_ROS1"].mask
    sciencemask |= cmx_mask["ROS_HA"].mask
    sciencemask |= cmx_mask["ROS_ROS2"].mask
    sciencemask |= cmx_mask["M33_STD_BRIGHT"].mask
    sciencemask |= cmx_mask["M33_H2PN"].mask
    sciencemask |= cmx_mask["M33_GC"].mask
    sciencemask |= cmx_mask["M33_QSO"].mask
    sciencemask |= cmx_mask["M33_M33cen"].mask
    sciencemask |= cmx_mask["M33_M33out"].mask
    sciencemask |= cmx_mask["MINI_SV_LRG"].mask
    sciencemask |= cmx_mask["MINI_SV_ELG"].mask
    sciencemask |= cmx_mask["MINI_SV_QSO"].mask
    sciencemask |= cmx_mask["MINI_SV_BGS_BRIGHT"].mask


    return sciencemask


def default_cmx_stdmask():
    """Returns default mask of bits for standards in CMX survey.
    """
    # Nothing in a CMX file is currently treated as a "standard".  The
    # objects are all things which should be assigned as science targets.
    stdmask = 0
    stdmask |= cmx_mask["STD_FAINT"].mask
    stdmask |= cmx_mask["STD_BRIGHT"].mask

    return stdmask


def default_cmx_skymask():
    """Returns default mask of bits for sky targets in CMX survey.
    """
    skymask = 0
    skymask |= cmx_mask["SKY"].mask
    return skymask


def default_cmx_suppskymask():
    """Returns default mask of bits for suppsky targets in CMX survey.
    """
    suppskymask = 0
    suppskymask |= cmx_mask["SUPP_SKY"].mask
    return suppskymask


def default_cmx_safemask():
    """Returns default mask of bits for 'safe' targets in CMX survey.

    Note: these are targets of last resort; they are safe locations where
    we won't saturate the detector, but aren't good for anything else.
    """
    safemask = 0
    safemask |= cmx_mask["BAD_SKY"].mask
    return safemask


def default_cmx_excludemask():
    """Returns default mask of bits for CMX survey targets to NOT observe.
    """
    excludemask = 0
    return excludemask


def desi_target_type(desi_target, mws_target, sciencemask, stdmask,
                     skymask, suppskymask, safemask, excludemask, gaia_stdmask):
    """Determine fiber assign type from the data column.

    Args:
        desi_target (iterable):  Scalar or array-like integer values.
        mws_target (iterable): Scalar or array-like integer values.
        sciencemask (int):  Integer value to bitwise-and when checking for
            science targets.
        stdmask (int):  Integer value to bitwise-and when checking for
            standards targets.
        skymask (int):  Integer value to bitwise-and when checking for
            sky targets.
        suppskymask (int):  Integer value to bitwise-and when checking for
            suppsky targets.
        safemask (int):  Integer value to bitwise-and when checking for
            safe targets.
        excludemask (int):  Integer value to bitwise-and when checking for
            targets to exclude.
        gaia_stdmask (int):  Integer value to bitwise-and when checking for
            Gaia standards targets.

    Returns:
        (array):  The fiberassign target types.

    """
    # print('sciencemask {}'.format(sciencemask))
    # print('stdmask     {}'.format(stdmask))
    # print('gaia_stdmask  {}'.format(gaia_stdmask))
    # print('skymask     {}'.format(skymask))
    # print('safemask    {}'.format(safemask))
    # print('excludemask {}'.format(excludemask))

    if np.isscalar(desi_target):
        ttype = 0
        if desi_target & sciencemask != 0:
            ttype |= TARGET_TYPE_SCIENCE
        if desi_target & stdmask != 0:
            ttype |= TARGET_TYPE_STANDARD
        if mws_target & gaia_stdmask != 0:
            ttype |= TARGET_TYPE_STANDARD
        if desi_target & skymask != 0:
            ttype |= TARGET_TYPE_SKY
        if desi_target & suppskymask != 0:
            ttype |= TARGET_TYPE_SUPPSKY
        if desi_target & safemask != 0:
            ttype |= TARGET_TYPE_SAFE
        if desi_target & excludemask != 0:
            ttype = 0
    else:
        desi_target = np.asarray(desi_target)
        mws_target = np.asarray(mws_target)
        ttype = np.zeros(len(desi_target), dtype=np.uint8)
        ttype[desi_target & sciencemask != 0] |= TARGET_TYPE_SCIENCE
        ttype[desi_target & stdmask != 0] |= TARGET_TYPE_STANDARD
        ttype[mws_target  & gaia_stdmask != 0] |= TARGET_TYPE_STANDARD
        ttype[desi_target & skymask != 0] |= TARGET_TYPE_SKY
        ttype[desi_target & suppskymask != 0] |= TARGET_TYPE_SUPPSKY
        ttype[desi_target & safemask != 0] |= TARGET_TYPE_SAFE
        ttype[desi_target & excludemask != 0] = 0

    return ttype


def default_survey_target_masks(survey, rundate=None):
    """Return the default masks for the survey.

    Args:
        survey (str): The survey name.
        rundate (optional, defaults to None): yyyy-mm-ddThh:mm:ss+00:00 rundate
            for focalplane with UTC timezone formatting (string)

    Returns:
        (tuple): The science mask, standard mask, sky mask, suppsky mask,
            safe mask, exclude mask, and gaia standard mask for the data.

    Notes:
        20210930 : include rundate argument, for default_main_stdmask().
    """
    sciencemask = None
    stdmask = None
    skymask = None
    suppskymask = None
    safemask = None
    excludemask = None
    gaia_stdmask = None
    if survey == "main":
        sciencemask = default_main_sciencemask()
        stdmask = default_main_stdmask(rundate=rundate)
        skymask = default_main_skymask()
        suppskymask = default_main_suppskymask()
        safemask = default_main_safemask()
        excludemask = default_main_excludemask()
        gaia_stdmask = default_main_gaia_stdmask()
    elif survey == "cmx":
        sciencemask = default_cmx_sciencemask()
        stdmask = default_cmx_stdmask()
        skymask = default_cmx_skymask()
        suppskymask = default_cmx_suppskymask()
        safemask = default_cmx_safemask()
        excludemask = default_cmx_excludemask()
    elif survey == "sv1":
        sciencemask = default_sv1_sciencemask()
        stdmask = default_sv1_stdmask()
        skymask = default_sv1_skymask()
        suppskymask = default_sv1_suppskymask()
        safemask = default_sv1_safemask()
        excludemask = default_sv1_excludemask()
    # AR duplicating for sv2...
    elif survey == "sv2":
        sciencemask = default_sv2_sciencemask()
        stdmask = default_sv2_stdmask()
        skymask = default_sv2_skymask()
        suppskymask = default_sv2_suppskymask()
        safemask = default_sv2_safemask()
        excludemask = default_sv2_excludemask()
    # AR duplicating for sv3...
    elif survey == "sv3":
        sciencemask = default_sv3_sciencemask()
        stdmask = default_sv3_stdmask()
        skymask = default_sv3_skymask()
        suppskymask = default_sv3_suppskymask()
        safemask = default_sv3_safemask()
        excludemask = default_sv3_excludemask()

    return (sciencemask, stdmask, skymask, suppskymask, safemask, excludemask, gaia_stdmask)


def default_target_masks(data, rundate=None):
    """Return the column name and default mask values for the data table.

    This identifies the type of target data and returns the defaults for
    the program type.

    Args:
        data (Table):  A Table or recarray.
        rundate (optional, defaults to None): yyyy-mm-ddThh:mm:ss+00:00 rundate for
            focalplane with UTC timezone formatting (string)

    Returns:
        (tuple):  The survey, column name, science mask, standard mask,
            sky mask, suppsky mask, safe mask, exclude mask, and
            gaia standard mask for the data.

    Notes:
        20210930 : include rundate argument, for default_main_stdmask().
    """
    col = None
    filecols, filemasks, filesurvey = main_cmx_or_sv(data)
    if filesurvey == "main":
        col = "DESI_TARGET"
    elif filesurvey == "cmx":
        col = filecols[0]
    elif filesurvey == "sv1":
        col = "SV1_DESI_TARGET"
    elif filesurvey == "sv2":
        col = "SV2_DESI_TARGET"
    elif filesurvey == "sv3":
        col = "SV3_DESI_TARGET"
    sciencemask, stdmask, skymask, suppskymask, safemask, excludemask, gaia_stdmask = \
        default_survey_target_masks(filesurvey, rundate=rundate)
    return (filesurvey, col, sciencemask, stdmask, skymask, suppskymask,
            safemask, excludemask, gaia_stdmask)


def append_target_table(tgs, tagalong, tgdata, survey, typeforce, typecol,
                        sciencemask,
                        stdmask, skymask, suppskymask, safemask, excludemask, gaia_stdmask):
    """Append a target recarray / table to a Targets object.

    This function is used to take a slice of targets table (as read from a
    file) and extract the columns containing properties which are stored
    internally in a Targets object.  These targets and their properties are
    added to the Targets object.

    Args:
        tgs (Targets): The targets object to modify.
        tagalong (TargetTagalong): a data structure that carries RA,Dec, and other information for targets from the targeting files, to be written to the fiberassign outputs.
        tgdata (Table): The table or recarray containing the input data.
        survey (str):  The survey type.
        typeforce (int): If not None, all targets are considered to be this
            type.
        typecol (str): The name of the column to use for bitmask operations.
        sciencemask (int):  Integer value to bitwise-and when checking for
            science targets.
        stdmask (int):  Integer value to bitwise-and when checking for
            standards targets.
        skymask (int):  Integer value to bitwise-and when checking for
            sky targets.
        suppskymask (int):  Integer value to bitwise-and when checking for
            suppsky targets.
        safemask (int):  Integer value to bitwise-and when checking for
            safe targets.
        excludemask (int):  Integer value to bitwise-and when checking for
            targets to exclude.
        gaia_stdmask (int):  Integer value to bitwise-and when checking for
            Gaia standards targets.

    Returns:
        None

    """

    log = Logger.get()

    validtypes = [
        TARGET_TYPE_SCIENCE,
        TARGET_TYPE_SKY,
        TARGET_TYPE_SUPPSKY,
        TARGET_TYPE_STANDARD,
        TARGET_TYPE_SAFE
    ]
    if typeforce is not None:
        if typeforce not in validtypes:
            raise RuntimeError("Cannot force objects to be an invalid type")
    # Create buffers for column data
    nrows = len(tgdata["TARGETID"][:])
    # Arrays needed by C++
    d_targetid = np.zeros(nrows, dtype=np.int64)
    d_type = np.zeros(nrows, dtype=np.uint8)
    d_nobs = np.zeros(nrows, dtype=np.int32)
    d_prior = np.zeros(nrows, dtype=np.int32)
    d_subprior = np.zeros(nrows, dtype=np.float64)
    # Arrays that have special handling
    d_ra = np.zeros(nrows, dtype=np.float64)
    d_dec = np.zeros(nrows, dtype=np.float64)
    d_bits = np.zeros(nrows, dtype=np.int64)
    d_obscond = np.zeros(nrows, dtype=np.int32)

    d_targetid[:] = tgdata["TARGETID"][:]
    if "TARGET_RA" in tgdata.dtype.names:
        d_ra[:] = tgdata["TARGET_RA"][:]
    else:
        d_ra[:] = tgdata["RA"][:]
    if "TARGET_DEC" in tgdata.dtype.names:
        d_dec[:] = tgdata["TARGET_DEC"][:]
    else:
        d_dec[:] = tgdata["DEC"][:]

    if typeforce is not None:
        d_type[:] = typeforce
        # In this case we leave the targets bits at zero since we are
        # forcibly assigning a type.  In this case, the target bits cannot
        # be used to determine anything about the object for QA, etc.
    else:
        if typecol == "FA_TYPE":
            # We are using the pre-established target categories.
            d_type[:] = tgdata["FA_TYPE"][:]
            d_bits[:] = tgdata["FA_TARGET"][:]
        else:
            d_bits[:] = tgdata[typecol][:]
            # AR trying to protect against cases where the *MWS_TARGET column
            # AR    is not present; that does not happen for the main
            # AR    survey, but imagining cases where fiberassign is
            # AR    run on some special targets set (e.g. only secondaries)
            mws_targets = 0 * tgdata[typecol]
            if typecol.replace("DESI", "MWS") in tgdata.dtype.names:
                mws_targets = tgdata[typecol.replace("DESI", "MWS")]
            d_type[:] = desi_target_type(
                tgdata[typecol], mws_targets, sciencemask, stdmask, skymask, suppskymask,
                safemask, excludemask, gaia_stdmask)

    if "OBSCONDITIONS" in tgdata.dtype.fields:
        d_obscond[:] = tgdata["OBSCONDITIONS"][:]
    else:
        # Set obs conditions mask to be all bits
        d_obscond[:] = np.invert(np.zeros(nrows, dtype=np.int32))

    if "NUMOBS_MORE" in tgdata.dtype.fields:
        d_nobs[:] = tgdata["NUMOBS_MORE"][:]
    elif "NUMOBS_INIT" in tgdata.dtype.fields:
        d_nobs[:] = tgdata["NUMOBS_INIT"][:]
    else:
        d_nobs[:] = np.zeros(nrows, dtype=np.int32)

    if "PRIORITY" in tgdata.dtype.fields:
        d_prior[:] = tgdata["PRIORITY"][:]
    elif "PRIORITY_INIT" in tgdata.dtype.fields:
        d_prior[:] = tgdata["PRIORITY_INIT"][:]
    else:
        d_prior[:] = np.zeros(nrows, dtype=np.int32)

    if "SUBPRIORITY" in tgdata.dtype.fields:
        d_subprior[:] = tgdata["SUBPRIORITY"][:]
    else:
        d_subprior[:] = np.zeros(nrows, dtype=np.float64)

    fake = {'TARGET_RA': d_ra,
            'TARGET_DEC': d_dec,
            'FA_TARGET': d_bits,
            'OBSCOND': d_obscond}
    if not 'PLATE_RA' in tgdata.dtype.fields:
        log.info('Warning: no PLATE_RA, PLATE_DEC in target file; using RA,DEC or TARGET_RA,DEC')
        fake.update({'PLATE_RA': d_ra,
                     'PLATE_DEC': d_dec,})
    tagalong.add_data(d_targetid, tgdata, fake=fake)

    # Append the data to our targets list.  This will print a
    # warning if there are duplicate target IDs.
    tgs.append(survey, d_targetid, d_nobs, d_prior, d_subprior, d_type)
    return


def load_target_table(tgs, tagalong, tgdata, survey=None, typeforce=None, typecol=None,
                      sciencemask=None, stdmask=None, skymask=None,
                      suppskymask=None, safemask=None, excludemask=None, gaia_stdmask=None,
                      rundate=None):
    """Append targets from a table.

    Use the table data to append targets to the input Targets object.
    A subset of the columns in the file will be stored in each Target added
    to the Targets object.  Each target is classified into one or more of the
    4 types used internally in assignment (science, standard, sky, safe).

    This classification is controlled by applying bitmasks to the specified
    data column.  Alternatively, all targets in the file can be forced to one
    type.

    Args:
        tgs (Targets): The targets object on which to append this data.
        tagalong (TargetTagalong): a data structure that carries RA,Dec, and other information for targets from the targeting files, to be written to the fiberassign outputs.
        tgdata (Table): A table or recarray with the target properties.
        survey (str):  The survey type.  If None, query from columns.
        typeforce (int): If specified, it must equal one of the TARGET_TYPE_*
            values.  All targets read from the file will be assigned this type.
        typecol (str): Optional column to use for bitmask matching (default
            uses the result of main_cmx_or_sv from desitarget).
        sciencemask (int): Bitmask for classifying targets as science.
        stdmask (int): Bitmask for classifying targets as a standard.
        skymask (int): Bitmask for classifying targets as sky.
        suppskymask (int): Bitmask for classifying targets as suppsky.
        safemask (int): Bitmask for classifying targets as a safe location.
        excludemask (int): Bitmask for excluding targets.
        gaia_stdmask (int): Bitmask for classifying targets as a Gaia standard.
        rundate (optional, defaults to None): yyyy-mm-ddThh:mm:ss+00:00 rundate
            for focalplane with UTC timezone formatting (string)

    Returns:
        None

    Notes:
        20210930 : include rundate argument, for default_main_stdmask().
    """
    log = Logger.get()
    if "TARGETID" not in tgdata.dtype.names:
        msg = "TARGETID column is required"
        log.error(msg)
        raise RuntimeError(msg)
    if tgdata.dtype["TARGETID"].char != "l":
        msg = "TARGETID column should be int64"
        log.error(msg)
        raise RuntimeError(msg)
    if "PRIORITY" in tgdata.dtype.names:
        if tgdata.dtype["PRIORITY"].char not in ["i", "l"]:
            msg = "PRIORITY column should be an integer type"
            log.error(msg)
            raise RuntimeError(msg)
    if "SUBPRIORITY" not in tgdata.dtype.names:
        msg = "SUBPRIORITY column is required"
        log.error(msg)
        raise RuntimeError(msg)
    if tgdata.dtype["SUBPRIORITY"].char != "d":
        msg = "SUBPRIORITY column should be float64"
        log.error(msg)
        raise RuntimeError(msg)
    if "NUMOBS_MORE" in tgdata.dtype.names:
        if tgdata.dtype["NUMOBS_MORE"].char not in ["i", "l"]:
            msg = "NUMOBS_MORE column should be an integer type"
            log.error(msg)
            raise RuntimeError(msg)
    if "NUMOBS_INIT" in tgdata.dtype.names:
        if tgdata.dtype["NUMOBS_INIT"].char not in ["i", "l"]:
            msg = "NUMOBS_INIT column should be an integer type"
            log.error(msg)
            raise RuntimeError(msg)
    if "OBSCONDITIONS" not in tgdata.dtype.names:
        msg = "OBSCONDITIONS column is required"
        log.error(msg)
        raise RuntimeError(msg)
    if tgdata.dtype["OBSCONDITIONS"].char not in ["i", "l"]:
        msg = "OBSCONDITIONS column should be an integer type"
        log.error(msg)
        raise RuntimeError(msg)

    # Are we loading raw output?  If so, we require the survey key to get
    # the default masks.
    fsurvey = None
    fcol = None
    fsciencemask = None
    fstdmask = None
    fskymask = None
    fsuppskymask = None
    fsafemask = None
    fexcludemask = None
    fgaia_stdmask = None
    if typecol == "FA_TYPE":
        if survey is None:
            msg = "When loading raw fiberassign tables, the survey must be \
                specified"
            log.error(msg)
            raise RuntimeError(msg)
        fsciencemask, fstdmask, fskymask, fsuppskymask, fsafemask, \
            fexcludemask, fgaia_stdmask = default_survey_target_masks(survey, rundate=rundate)
    else:
        fsurvey, fcol, fsciencemask, fstdmask, fskymask, fsuppskymask, \
            fsafemask, fexcludemask, fgaia_stdmask = default_target_masks(tgdata, rundate=rundate)
        if fcol is None:
            # File could not be identified.  In this case, the user must
            # completely specify the bitmask and column to use.
            if typeforce is None:
                if (typecol is None) or (sciencemask is None) \
                        or (stdmask is None) or (skymask is None) \
                        or (suppskymask is None) or (safemask is None) \
                        or (excludemask is None) or (gaia_stdmask is None):
                    msg = "Unknown survey type.  To use this table, \
                        specify the column name and every bitmask."
                    log.error(msg)
                    raise RuntimeError(msg)

    if survey is None:
        survey = fsurvey
    if typecol is None:
        typecol = fcol
    if sciencemask is None:
        sciencemask = fsciencemask
    if stdmask is None:
        stdmask = fstdmask
    if skymask is None:
        skymask = fskymask
    if suppskymask is None:
        suppskymask = fsuppskymask
    if safemask is None:
        safemask = fsafemask
    if excludemask is None:
        excludemask = fexcludemask
    if gaia_stdmask is None:
        gaia_stdmask = fgaia_stdmask

    log.debug("Target table using survey '{}', column {}:".format(survey, typecol))
    if survey == "main":
        log.debug("  sciencemask {}".format(
            "|".join(desi_mask.names(sciencemask))))
        log.debug("  stdmask     {}".format(
            "|".join(desi_mask.names(stdmask))))
        log.debug("  skymask     {}".format(
            "|".join(desi_mask.names(skymask))))
        log.debug("  suppskymask     {}".format(
            "|".join(desi_mask.names(suppskymask))))
        log.debug("  safemask    {}".format(
            "|".join(desi_mask.names(safemask))))
        log.debug("  excludemask {}".format(
            "|".join(desi_mask.names(excludemask))))
        log.debug("  gaia_stdmask {}".format(
            "|".join(desi_mask.names(gaia_stdmask))))
    elif survey == "cmx":
        log.debug("  sciencemask {}".format(
            "|".join(cmx_mask.names(sciencemask))))
        log.debug("  stdmask     {}".format(
            "|".join(cmx_mask.names(stdmask))))
        log.debug("  skymask     {}".format(
            "|".join(cmx_mask.names(skymask))))
        log.debug("  suppskymask     {}".format(
            "|".join(cmx_mask.names(suppskymask))))
        log.debug("  safemask    {}".format(
            "|".join(cmx_mask.names(safemask))))
        log.debug("  excludemask {}".format(
            "|".join(cmx_mask.names(excludemask))))
    elif survey == "sv1":
        log.debug("  sciencemask {}".format(
            "|".join(sv1_mask.names(sciencemask))))
        log.debug("  stdmask     {}".format(
            "|".join(sv1_mask.names(stdmask))))
        log.debug("  skymask     {}".format(
            "|".join(sv1_mask.names(skymask))))
        log.debug("  suppskymask     {}".format(
            "|".join(sv1_mask.names(suppskymask))))
        log.debug("  safemask    {}".format(
            "|".join(sv1_mask.names(safemask))))
        log.debug("  excludemask {}".format(
            "|".join(sv1_mask.names(excludemask))))
    # AR adding sv2...
    elif survey == "sv2":
        log.debug("  sciencemask {}".format(
            "|".join(sv2_mask.names(sciencemask))))
        log.debug("  stdmask     {}".format(
            "|".join(sv2_mask.names(stdmask))))
        log.debug("  skymask     {}".format(
            "|".join(sv2_mask.names(skymask))))
        log.debug("  suppskymask     {}".format(
            "|".join(sv2_mask.names(suppskymask))))
        log.debug("  safemask    {}".format(
            "|".join(sv2_mask.names(safemask))))
        log.debug("  excludemask {}".format(
            "|".join(sv2_mask.names(excludemask))))
    # AR adding sv3...
    elif survey == "sv3":
        log.debug("  sciencemask {}".format(
            "|".join(sv3_mask.names(sciencemask))))
        log.debug("  stdmask     {}".format(
            "|".join(sv3_mask.names(stdmask))))
        log.debug("  skymask     {}".format(
            "|".join(sv3_mask.names(skymask))))
        log.debug("  suppskymask     {}".format(
            "|".join(sv3_mask.names(suppskymask))))
        log.debug("  safemask    {}".format(
            "|".join(sv3_mask.names(safemask))))
        log.debug("  excludemask {}".format(
            "|".join(sv3_mask.names(excludemask))))
    else:
        raise RuntimeError("unknown survey type, should never get here!")
    append_target_table(tgs, tagalong, tgdata, survey, typeforce, typecol, sciencemask,
                        stdmask, skymask, suppskymask, safemask, excludemask, gaia_stdmask)
    return


def load_target_file(tgs, tagalong, tfile, survey=None, typeforce=None, typecol=None,
                     sciencemask=None, stdmask=None, skymask=None,
                     suppskymask=None, safemask=None, excludemask=None,
                     rowbuffer=1000000, gaia_stdmask=None, rundate=None):
    """Append targets from a file.

    Read the specified file and append targets to the input Targets object.
    A subset of the columns in the file will be stored in each Target added
    to the Targets object.  Each target is classified into one or more of the
    4 types used internally in assignment (science, standard, sky, safe).

    This classification is controlled by applying bitmasks to the specified
    data column.  Alternatively, all targets in the file can be forced to one
    type.

    Args:
        tgs (Targets): The targets object on which to append this data.
        tagalong (TargetTagalong): a data structure that carries RA,Dec, and other information for targets from the targeting files, to be written to the fiberassign outputs.  A new one can be created using `fiberassign.targets.create_tagalong`.
        tfile (str): The path to the target catalog.
        survey (str):  The survey type.  If None, query from columns and
            the FITS header.
        typeforce (int): If specified, it must equal one of the TARGET_TYPE_*
            values.  All targets read from the file will be assigned this type.
        typecol (str): Optional column to use for bitmask matching (default
            uses the result of main_cmx_or_sv from desitarget).
        sciencemask (int): Bitmask for classifying targets as science.
        stdmask (int): Bitmask for classifying targets as a standard.
        skymask (int): Bitmask for classifying targets as sky.
        suppskymask (int): Bitmask for classifying targets as suppsky.
        safemask (int): Bitmask for classifying targets as a safe location.
        excludemask (int): Bitmask for excluding targets.
        rowbuffer (int): Optional number of rows to read at once when loading
            very large files.
        gaia_stdmask (int): Bitmask for classifying targets as a Gaia standard.
        rundate (optional, defaults to None): yyyy-mm-ddThh:mm:ss+00:00 rundate
            for focalplane with UTC timezone formatting (string)

    Returns:
        (str): The survey type.

    Notes:
        20210930 : include rundate argument, for default_main_stdmask().
    """
    tm = Timer()
    tm.start()

    log = Logger.get()

    # Open file
    fits = fitsio.FITS(tfile, mode="r")

    # Total number of rows
    nrows = fits[1].get_nrows()
    log.info("Target file {} has {} rows.  Reading in chunks of {}".format(tfile, nrows, rowbuffer))

    header = fits[1].read_header()
    if survey is None:
        if "FA_SURV" in header:
            survey = str(header["FA_SURV"]).rstrip()

    offset = 0
    n = rowbuffer
    while offset < nrows:
        if offset + n > nrows:
            n = nrows - offset
        data = fits[1].read(rows=np.arange(offset, offset + n, dtype=np.int64))
        log.debug("Target file {} read rows {} - {}".format(tfile, offset, offset + n - 1))
        load_target_table(tgs, tagalong, data, survey=survey,
                          typeforce=typeforce,
                          typecol=typecol,
                          sciencemask=sciencemask,
                          stdmask=stdmask,
                          skymask=skymask,
                          suppskymask=suppskymask,
                          safemask=safemask,
                          excludemask=excludemask,
                          gaia_stdmask=gaia_stdmask,
                          rundate=rundate)
        offset += n

    tm.stop()
    tm.report("Read target file {}".format(tfile))

    return survey

def targets_in_tiles(hw, tgs, tiles, tagalong):
    '''
    Returns tile_targetids, tile_x, tile_y,
    which are maps from tileid to numpy arrays.
    '''

    log = Logger.get()

    try:
        plate_scale = desimodel.io.load_platescale() # this loads the numbers in $DESIMODEL
        focalplane_radius_deg = plate_scale['theta'].max() 
    except:
        focalplane_radius_deg = hw.focalplane_radius_deg # this loads the numbers hardcoded in fiberassign
    
    log.info(f'Using {focalplane_radius_deg} as the focal plane radius in degrees')
    
    tile_targetids = {}
    tile_x = {}
    tile_y = {}
    tile_xy_cs5 = {}

    target_ids = tgs.ids()
    target_ra, target_dec, target_obscond = tagalong.get_for_ids(
        target_ids, ['RA', 'DEC', 'OBSCOND'])

    kd = _radec2kd(target_ra, target_dec)

    for (tile_id, tile_ra, tile_dec, tile_obscond, tile_ha, tile_obstheta,
         tile_obstime) in zip(
            tiles.id, tiles.ra, tiles.dec, tiles.obscond, tiles.obshourang,
            tiles.obstheta, tiles.obstime):

        log.info(f'Tile {tile_id} at RA,Dec {tile_ra},{tile_dec} obscond: {tile_obscond} HA: {tile_ha} obstime: {tile_obstime}')

        inds = _kd_query_radec(kd, tile_ra, tile_dec, focalplane_radius_deg)
        match = np.flatnonzero(target_obscond[inds] & tile_obscond)
        inds = inds[match]
        del match
        ras  = target_ra [inds]
        decs = target_dec[inds]
        tids = target_ids[inds]
        del inds
        log.info(f'Found {len(tids)} targets near tile and matching obscond')

        x, y = radec2xy(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta,
                        tile_ha, ras, decs, True)
        # Save CS5 mapping
        tile_xy_cs5[tile_id] = dict((tid,(xi,yi)) for tid,xi,yi in zip(tids, x, y))
        x, y = cs52xy(x, y)

        tile_targetids[tile_id] = tids
        tile_x[tile_id] = x
        tile_y[tile_id] = y

    return tile_targetids, tile_x, tile_y, tile_xy_cs5


def _radec2kd(ra, dec):
    """
    Creates a scipy KDTree from the given *ra*, *dec* arrays (in deg).
    """
    from scipy.spatial import KDTree
    xyz = _radec2xyz(ra, dec)
    return KDTree(xyz)


def _radec2xyz(ra, dec):
    """
    Converts arrays from *ra*, *dec* (in deg) to XYZ unit-sphere
    coordinates.
    """
    rr = np.deg2rad(ra)
    dd = np.deg2rad(dec)
    return np.vstack((np.cos(rr) * np.cos(dd),
                      np.sin(rr) * np.cos(dd),
                      np.sin(dd))).T

def _kd_query_radec(kd, ra, dec, radius_deg):
    searchrad = np.deg2rad(radius_deg)
    # Convert from radius to (tangent) distance on the unit sphere.
    searchrad = np.sqrt(2. * (1. - np.cos(searchrad)))
    xyz = _radec2xyz([ra], [dec])
    inds = kd.query_ball_point(xyz[0, :], searchrad)
    inds = np.array(inds)
    return inds
