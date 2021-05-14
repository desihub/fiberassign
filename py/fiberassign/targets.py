# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.targets
=====================

Functions for loading the target list

"""
from __future__ import absolute_import, division, print_function

import numpy as np

import fitsio

# FIXME:  If / when SV bit names diverge from main survey names, we
# should import the SV bitmasks here.
# AR duplicating what s required for sv2;
# AR should definitely be re-written to handle sv3, etc
# AR and duplicating for sv3...

from desitarget.targetmask import desi_mask

from desitarget.cmx.cmx_targetmask import cmx_mask

from desitarget.sv1.sv1_targetmask import desi_mask as sv1_mask
from desitarget.sv1.sv1_targetmask import scnd_mask as sv1_scnd_mask
# AR
from desitarget.sv2.sv2_targetmask import desi_mask as sv2_mask
from desitarget.sv3.sv3_targetmask import desi_mask as sv3_mask

from desitarget.targets import main_cmx_or_sv

from .utils import Logger, Timer
from .hardware import radec2xy

from ._internal import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                        TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                        TARGET_TYPE_SUPPSKY,
                        Target, Targets, TargetsAvailable,
                        LocationsAvailable)


class TargetTagalong(object):
    def __init__(self, columns, outnames={}):
        self.columns = columns
        self.outnames = outnames
        self.data = []

    def get_default(self, column):
        return None

    def get_output_name(self, column):
        return self.outnames.get(column, column)

    def add_data(self, targetids, tabledata, fake={}):
        tgarrays = [targetids]
        for k in self.columns:
            if k in fake:
                tgarrays.append(fake[k])
            else:
                 tgarrays.append(tabledata[k])
        self.data.append(tgarrays)

    def set_data(self, targetids, tabledata):
        '''
        Sets *ALL* rows of the given *tabledata* object to defaults,
        and then fills in values for the given targetids, if they are found.
        '''
        # Set defaults
        for c in self.columns:
            defval = self.get_default(c)
            if defval is None:
                continue
            outarr = tabledata[self.get_output_name(c)]
            outarr[:] = defval

        # Build output targetid-to-index map
        outmap = dict([(tid,i) for i,tid in enumerate(targetids)])
        # Put tagalong data into these output arrays
        outarrs = [tabledata[self.get_output_name(c)] for c in self.columns]
        # Go through my many data arrays
        for thedata in self.data:
            tids = thedata[0]
            # Search for output array indices for these targetids
            outinds = np.array([outmap.get(tid, -1) for tid in tids])
            ininds = np.flatnonzero(outinds >= 0)
            outinds = outinds[ininds]

            for outarr,inarr in zip(outarrs, thedata[1:]):
                outarr[outinds] = inarr[ininds]

    def get(self, name):
        '''
        Fetch appended array of the given name.
        '''
        arrs = []
        i = self.columns.index(name)
        for thedata in self.data:
            # +1 because 'thedata' starts with TARGETID.
            arrs.append(thedata[i+1])
        return np.hstack(arrs)

    def get_for_ids(self, targetids, names):
        '''
        Fetch arrays for the given names and given targetids.
        '''
        # Create output arrays
        outarrs = []
        colinds = []
        for name in names:
            ic = self.columns.index(name)
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
    """
    sciencemask = 0
    sciencemask |= desi_mask["LRG"].mask
    sciencemask |= desi_mask["ELG"].mask
    sciencemask |= desi_mask["QSO"].mask
    sciencemask |= desi_mask["BGS_ANY"].mask
    sciencemask |= desi_mask["MWS_ANY"].mask
    if "SCND_ANY" in desi_mask.names():
        sciencemask |= desi_mask["SCND_ANY"].mask
    else:
        sciencemask |= desi_mask["SECONDARY_ANY"].mask
    return sciencemask


def default_main_stdmask():
    """Returns default mask of bits for standards in main survey.
    """
    stdmask = 0
    stdmask |= desi_mask["STD_FAINT"].mask
    stdmask |= desi_mask["STD_WD"].mask
    stdmask |= desi_mask["STD_BRIGHT"].mask
    return stdmask


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


def desi_target_type(desi_target, sciencemask, stdmask,
                     skymask, suppskymask, safemask, excludemask):
    """Determine fiber assign type from the data column.

    Args:
        desi_target (iterable):  Scalar or array-like integer values.
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

    Returns:
        (array):  The fiberassign target types.

    """
    # print('sciencemask {}'.format(sciencemask))
    # print('stdmask     {}'.format(stdmask))
    # print('skymask     {}'.format(skymask))
    # print('safemask    {}'.format(safemask))
    # print('excludemask {}'.format(excludemask))

    if np.isscalar(desi_target):
        ttype = 0
        if desi_target & sciencemask != 0:
            ttype |= TARGET_TYPE_SCIENCE
        if desi_target & stdmask != 0:
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
        ttype = np.zeros(len(desi_target), dtype=np.uint8)
        ttype[desi_target & sciencemask != 0] |= TARGET_TYPE_SCIENCE
        ttype[desi_target & stdmask != 0] |= TARGET_TYPE_STANDARD
        ttype[desi_target & skymask != 0] |= TARGET_TYPE_SKY
        ttype[desi_target & suppskymask != 0] |= TARGET_TYPE_SUPPSKY
        ttype[desi_target & safemask != 0] |= TARGET_TYPE_SAFE
        ttype[desi_target & excludemask != 0] = 0

    return ttype


def default_survey_target_masks(survey):
    """Return the default masks for the survey.

    Args:
        survey (str): The survey name.

    Returns:
        (tuple): The science mask, standard mask, sky mask, suppsky mask,
            safe mask, and exclude mask for the data.

    """
    sciencemask = None
    stdmask = None
    skymask = None
    suppskymask = None
    safemask = None
    excludemask = None
    if survey == "main":
        sciencemask = default_main_sciencemask()
        stdmask = default_main_stdmask()
        skymask = default_main_skymask()
        suppskymask = default_main_suppskymask()
        safemask = default_main_safemask()
        excludemask = default_main_excludemask()
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

    return (sciencemask, stdmask, skymask, suppskymask, safemask, excludemask)


def default_target_masks(data):
    """Return the column name and default mask values for the data table.

    This identifies the type of target data and returns the defaults for
    the program type.

    Args:
        data (Table):  A Table or recarray.

    Returns:
        (tuple):  The survey, column name, science mask, standard mask,
            sky mask, suppsky mask, safe mask, and exclude mask for the data.

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
    sciencemask, stdmask, skymask, suppskymask, safemask, excludemask = \
        default_survey_target_masks(filesurvey)
    return (filesurvey, col, sciencemask, stdmask, skymask, suppskymask,
            safemask, excludemask)


def append_target_table(tgs, tgdata, survey, typeforce, typecol, sciencemask,
                        stdmask, skymask, suppskymask, safemask, excludemask,
                        tagalong=None):
    """Append a target recarray / table to a Targets object.

    This function is used to take a slice of targets table (as read from a
    file) and extract the columns containing properties which are stored
    internally in a Targets object.  These targets and their properties are
    added to the Targets object.

    Args:
        tgs (Targets): The targets object to modify.
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

    Returns:
        None

    """
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
    d_obscond = np.zeros(nrows, dtype=np.int32)
    d_targetid = np.zeros(nrows, dtype=np.int64)
    d_ra = np.zeros(nrows, dtype=np.float64)
    d_dec = np.zeros(nrows, dtype=np.float64)
    d_pra = np.zeros(nrows, dtype=np.float64)
    d_pdec = np.zeros(nrows, dtype=np.float64)
    d_pepoch = np.zeros(nrows, dtype=np.float64)
    d_type = np.zeros(nrows, dtype=np.uint8)
    d_nobs = np.zeros(nrows, dtype=np.int32)
    d_prior = np.zeros(nrows, dtype=np.int32)
    d_subprior = np.zeros(nrows, dtype=np.float64)
    d_targetid[:] = tgdata["TARGETID"][:]

    d_bits = np.zeros(nrows, dtype=np.int64)

    if "TARGET_RA" in tgdata.dtype.names:
        d_ra[:] = tgdata["TARGET_RA"][:]
    else:
        d_ra[:] = tgdata["RA"][:]
    if "TARGET_DEC" in tgdata.dtype.names:
        d_dec[:] = tgdata["TARGET_DEC"][:]
    else:
        d_dec[:] = tgdata["DEC"][:]

    if "PLATE_RA" in tgdata.dtype.names:
        print('Copying PLATE_RA from input file')
        d_pra[:] = tgdata["PLATE_RA"][:]
    else:
        print('Copying PLATE_RA from RA/TARGET_RA')
        d_pra[:] = d_ra[:]
    if "PLATE_DEC" in tgdata.dtype.names:
        print('Copying PLATE_DEC from input file')
        d_pdec[:] = tgdata["PLATE_DEC"][:]
    else:
        print('Copying PLATE_DEC from DEC/TARGET_DEC')
        d_pdec[:] = d_dec[:]
    if "PLATE_REF_EPOCH" in tgdata.dtype.names:
        print('Copying PLATE_REF_EPOCH from input file')
        d_pepoch[:] = tgdata["PLATE_REF_EPOCH"][:]
    else:
        print('No PLATE_REF_EPOCH')

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
            d_type[:] = desi_target_type(
                tgdata[typecol], sciencemask, stdmask, skymask, suppskymask,
                safemask, excludemask)

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

    if tagalong is not None:
        tagalong.add_data(d_targetid, tgdata, fake={'FA_TARGET': d_bits})

    # Append the data to our targets list.  This will print a
    # warning if there are duplicate target IDs.
    tgs.append(survey, d_targetid, d_nobs, d_prior,
               d_subprior, d_obscond, d_type)
    return


def load_target_table(tgs, tgdata, survey=None, typeforce=None, typecol=None,
                      sciencemask=None, stdmask=None, skymask=None,
                      suppskymask=None, safemask=None, excludemask=None,
                      tagalong=None):
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

    Returns:
        None

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
    if typecol == "FA_TYPE":
        if survey is None:
            msg = "When loading raw fiberassign tables, the survey must be \
                specified"
            log.error(msg)
            raise RuntimeError(msg)
        fsciencemask, fstdmask, fskymask, fsuppskymask, fsafemask, \
            fexcludemask = default_survey_target_masks(survey)
    else:
        fsurvey, fcol, fsciencemask, fstdmask, fskymask, fsuppskymask, \
            fsafemask, fexcludemask = default_target_masks(tgdata)
        if fcol is None:
            # File could not be identified.  In this case, the user must
            # completely specify the bitmask and column to use.
            if typeforce is None:
                if (typecol is None) or (sciencemask is None) \
                        or (stdmask is None) or (skymask is None) \
                        or (suppskymask is None) or (safemask is None) \
                        or (excludemask is None):
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

    log.debug("Target table using survey '{}', column {}:".format(
        survey, typecol))
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
    append_target_table(tgs, tgdata, survey, typeforce, typecol, sciencemask,
                        stdmask, skymask, suppskymask, safemask, excludemask,
                        tagalong=tagalong)
    return


def load_target_file(tgs, tfile, survey=None, typeforce=None, typecol=None,
                     sciencemask=None, stdmask=None, skymask=None,
                     suppskymask=None, safemask=None, excludemask=None,
                     rowbuffer=1000000, tagalong=None):
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

    Returns:
        (str): The survey type.

    """
    tm = Timer()
    tm.start()

    log = Logger.get()

    # Open file
    fits = fitsio.FITS(tfile, mode="r")

    # Total number of rows
    nrows = fits[1].get_nrows()
    log.info("Target file {} has {} rows.  Reading in chunks of {}"
             .format(tfile, nrows, rowbuffer))

    header = fits[1].read_header()
    if survey is None:
        if "FA_SURV" in header:
            survey = str(header["FA_SURV"]).rstrip()

    offset = 0
    n = rowbuffer
    while offset < nrows:
        if offset + n > nrows:
            n = nrows - offset
        data = fits[1].read(rows=np.arange(offset, offset+n, dtype=np.int64))
        log.debug("Target file {} read rows {} - {}"
                  .format(tfile, offset, offset+n-1))
        load_target_table(tgs, data, survey=survey,
                          typeforce=typeforce,
                          typecol=typecol,
                          sciencemask=sciencemask,
                          stdmask=stdmask,
                          skymask=skymask,
                          suppskymask=suppskymask,
                          safemask=safemask,
                          excludemask=excludemask,
                          tagalong=tagalong)
        offset += n

    tm.stop()
    tm.report("Read target file {}".format(tfile))

    return survey

def targets_in_tiles(hw, tgs, tiles, tagalong):
    '''
    Returns tile_targetids, tile_x, tile_y
    '''
    tile_targetids = {}
    tile_x = {}
    tile_y = {}

    #target_ra = tagalong.get('RA')
    #target_dec = tagalong.get('DEC')

    target_ids = tgs.ids()
    target_ra, target_dec = tagalong.get_for_ids(target_ids, ['RA', 'DEC'])

    #target_ra  = np.zeros(len(target_ids))
    #target_dec = np.zeros(len(target_ids))
    target_obscond = np.zeros(len(target_ids), np.int32)
    for i,tid in enumerate(target_ids):
        tg = tgs.get(tid)
        #target_ra[i] = tg.ra
        #target_dec[i] = tg.dec
        target_obscond[i] = tg.obscond
    #target_ra  = np.array([tgs.get(tid).ra  for tid in target_ids])
    #target_dec = np.array([tgs.get(tid).dec for tid in target_ids])
    kd = _radec2kd(target_ra, target_dec)

    for (tile_id, tile_ra, tile_dec, tile_obscond, tile_ha, tile_obstheta,
         tile_obstime) in zip(
            tiles.id, tiles.ra, tiles.dec, tiles.obscond, tiles.obshourang,
            tiles.obstheta, tiles.obstime):

        print('Tile', tile_id, 'at RA,Dec', tile_ra, tile_dec, 'obscond:', tile_obscond, 'HA', tile_ha, 'obstime', tile_obstime)

        inds = _kd_query_radec(kd, tile_ra, tile_dec, hw.focalplane_radius_deg)
        match = np.flatnonzero(target_obscond[inds] & tile_obscond)
        inds = inds[match]
        del match
        ras  = target_ra [inds]
        decs = target_dec[inds]
        tids = target_ids[inds]
        del inds
        print('Found', len(tids), 'targets near tile and matching obscond')

        fx,fy = radec2xy(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta,
                         tile_ha, ras, decs, False)

        tile_targetids[tile_id] = tids
        tile_x[tile_id] = fx
        tile_y[tile_id] = fy

    return tile_targetids, tile_x, tile_y


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
