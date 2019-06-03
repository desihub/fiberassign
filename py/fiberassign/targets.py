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

from desitarget.targetmask import desi_mask

from desitarget.cmx.cmx_targetmask import cmx_mask

from desitarget.sv1.sv1_targetmask import desi_mask as sv1_mask

from desitarget.targets import main_cmx_or_sv

from .utils import Logger, Timer

from ._internal import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                        TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                        Target, Targets, TargetTree, TargetsAvailable,
                        LocationsAvailable)


def str_to_target_type(input):
    if input == "science":
        return TARGET_TYPE_SCIENCE
    elif input == "sky":
        return TARGET_TYPE_SKY
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


def default_sv1_sciencemask():
    """Returns default mask of bits for science targets in SV1 survey.
    """
    sciencemask = 0
    sciencemask |= sv1_mask["LRG"].mask
    sciencemask |= sv1_mask["ELG"].mask
    sciencemask |= sv1_mask["QSO"].mask
    sciencemask |= sv1_mask["BGS_ANY"].mask
    sciencemask |= sv1_mask["MWS_ANY"].mask
    sciencemask |= sv1_mask["SECONDARY_ANY"].mask
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
    return sciencemask


def default_cmx_stdmask():
    """Returns default mask of bits for standards in CMX survey.
    """
    # Nothing in a CMX file is currently treated as a "standard".  The
    # objects are all things which should be assigned as science targets.
    stdmask = 0
    return stdmask


def default_cmx_skymask():
    """Returns default mask of bits for sky targets in CMX survey.
    """
    skymask = 0
    skymask |= cmx_mask["SKY"].mask
    return skymask


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
                     skymask, safemask, excludemask):
    """Determine fiber assign type from the data column.

    Args:
        desi_target (iterable):  Scalar or array-like integer values.
        sciencemask (int):  Integer value to bitwise-and when checking for
            science targets.
        stdmask (int):  Integer value to bitwise-and when checking for
            standards targets.
        skymask (int):  Integer value to bitwise-and when checking for
            sky targets.
        safemask (int):  Integer value to bitwise-and when checking for
            safe targets.
        excludemask (int):  Integer value to bitwise-and when checking for
            targets to exclude.

    Returns:
        (array):  The fiberassign target types.

    """
    if np.isscalar(desi_target):
        ttype = 0
        if desi_target & sciencemask != 0:
            ttype |= TARGET_TYPE_SCIENCE
        if desi_target & stdmask != 0:
            ttype |= TARGET_TYPE_STANDARD
        if desi_target & skymask != 0:
            ttype |= TARGET_TYPE_SKY
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
        ttype[desi_target & safemask != 0] |= TARGET_TYPE_SAFE
        ttype[desi_target & excludemask != 0] = 0

    return ttype


def default_survey_target_masks(survey):
    """Return the default masks for the survey.

    Args:
        survey (str): The survey name.

    Returns:
        (tuple): The science mask, standard mask, sky mask, safe mask,
            and exclude mask for the data.

    """
    sciencemask = None
    stdmask = None
    skymask = None
    safemask = None
    excludemask = None
    if survey == "main":
        sciencemask = default_main_sciencemask()
        stdmask = default_main_stdmask()
        skymask = default_main_skymask()
        safemask = default_main_safemask()
        excludemask = default_main_excludemask()
    elif survey == "cmx":
        sciencemask = default_cmx_sciencemask()
        stdmask = default_cmx_stdmask()
        skymask = default_cmx_skymask()
        safemask = default_cmx_safemask()
        excludemask = default_cmx_excludemask()
    elif survey == "sv1":
        sciencemask = default_sv1_sciencemask()
        stdmask = default_sv1_stdmask()
        skymask = default_sv1_skymask()
        safemask = default_sv1_safemask()
        excludemask = default_sv1_excludemask()
    return (sciencemask, stdmask, skymask, safemask, excludemask)


def default_target_masks(data):
    """Return the column name and default mask values for the data table.

    This identifies the type of target data and returns the defaults for
    the program type.

    Args:
        data (Table):  A Table or recarray.

    Returns:
        (tuple):  The survey, column name, science mask, standard mask,
            sky mask, safe mask, and exclude mask for the data.

    """
    col = None
    filecols, filemasks, filesurvey = main_cmx_or_sv(data)
    if filesurvey == "main":
        col = "DESI_TARGET"
    elif filesurvey == "cmx":
        col = filecols[0]
    elif filesurvey == "sv1":
        col = "SV1_DESI_TARGET"
    sciencemask, stdmask, skymask, safemask, excludemask = \
        default_survey_target_masks(filesurvey)
    return (filesurvey, col, sciencemask, stdmask, skymask, safemask,
            excludemask)


def append_target_table(tgs, tgdata, survey, typeforce, typecol, sciencemask,
                        stdmask, skymask, safemask, excludemask):
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
    d_bits = np.zeros(nrows, dtype=np.int64)
    d_type = np.zeros(nrows, dtype=np.uint8)
    d_nobs = np.zeros(nrows, dtype=np.int32)
    d_prior = np.zeros(nrows, dtype=np.int32)
    d_subprior = np.zeros(nrows, dtype=np.float64)
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
            d_type[:] = desi_target_type(
                tgdata[typecol], sciencemask, stdmask, skymask, safemask,
                excludemask)

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

    # Append the data to our targets list.  This will print a
    # warning if there are duplicate target IDs.
    tgs.append(survey, d_targetid, d_ra, d_dec, d_bits, d_nobs, d_prior,
               d_subprior, d_obscond, d_type)
    return


def load_target_table(tgs, tgdata, survey=None, typeforce=None, typecol=None,
                      sciencemask=None, stdmask=None, skymask=None,
                      safemask=None, excludemask=None):
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
    fsafemask = None
    fexcludemask = None
    if typecol == "FA_TYPE":
        if survey is None:
            msg = "When loading raw fiberassign tables, the survey must be \
                specified"
            log.error(msg)
            raise RuntimeError(msg)
        fsciencemask, fstdmask, fskymask, fsafemask, fexcludemask = \
            default_survey_target_masks(survey)
    else:
        fsurvey, fcol, fsciencemask, fstdmask, fskymask, fsafemask, \
            fexcludemask = default_target_masks(tgdata)
        if fcol is None:
            # File could not be identified.  In this case, the user must
            # completely specify the bitmask and column to use.
            if typeforce is None:
                if (typecol is None) or (sciencemask is None) \
                        or (stdmask is None) or (skymask is None) \
                        or (safemask is None) or (excludemask is None):
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
        log.debug("  safemask    {}".format(
            "|".join(sv1_mask.names(safemask))))
        log.debug("  excludemask {}".format(
            "|".join(sv1_mask.names(excludemask))))
    else:
        raise RuntimeError("unknown survey type, should never get here!")
    append_target_table(tgs, tgdata, survey, typeforce, typecol, sciencemask,
                        stdmask, skymask, safemask, excludemask)
    return


def load_target_file(tgs, tfile, survey=None, typeforce=None, typecol=None,
                     sciencemask=None, stdmask=None, skymask=None,
                     safemask=None, excludemask=None, rowbuffer=1000000):
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
                          stdmask=stdmask, skymask=skymask,
                          safemask=safemask,
                          excludemask=excludemask)
        offset += n

    tm.stop()
    tm.report("Read target file {}".format(tfile))

    return survey
