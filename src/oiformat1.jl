#
# oiformat1.jl --
#
# Define 1st revision of OI-FITS format.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015-2019: Éric Thiébaut.
#
#------------------------------------------------------------------------------


# OI-FITS FORMAT DESCRIPTION TABLES
#
# The format of the OI-FITS data block is described by a vector of strings
# like:
#
#   ["KEYWORD FORMAT DESCR",
#     ...,
#     ...,
#    "---------------------------",
#    "COLUMN  FORMAT DESCR",
#     ...,
#     ...,
#     ...]
#
# where:
#
#   KEYWORD = keyword for HDU header.
#   COLUMN = column name for table (TTYPE).
#   FORMAT = (may be prefixed with a ? to indicate optional field):
#       for keywords, a single letter indicating the type,
#       for columns, a letter indicating the type followed by list of
#                    dimensions in parenthesis (letter W for a dimension means
#                    NWAVE).
#   DESCR = description/comment; units, if any are indicated at the end
#           between square brackets.
#
# There may be any number of keyword definitions and any number of column
# definitions, the two parts are separated by a dash line like
# "--------------".
#

# OI_TARGET definition (1st revision):
add_def("OI_TARGET", 1,
        ["OI_REVN    I      revision number of the table definition",
         "---------------------------------------------------------",
         "TARGET_ID  I(1)   index number",
         "TARGET     A(16)  target name",
         "RAEP0      D(1)   RA at mean equinox [deg]",
         "DECEP0     D(1)   DEC at mean equinox [deg]",
         "EQUINOX    E(1)   equinox [yr]",
         "RA_ERR     D(1)   error in RA at mean equinox [deg]",
         "DEC_ERR    D(1)   error in DEC at mean equino [deg]",
         "SYSVEL     D(1)   systemic radial velocity [m/s]",
         "VELTYP     A(8)   reference for radial velocity",
         "VELDEF     A(8)   definition of radial velocity",
         "PMRA       D(1)   proper motion in RA [deg/yr]",
         "PMDEC      D(1)   proper motion in DEC [deg/yr]",
         "PMRA_ERR   D(1)   error of proper motion in RA [deg/yr]",
         "PMDEC_ERR  D(1)   error of proper motion in DEC [deg/yr]",
         "PARALLAX   E(1)   parallax [deg]",
         "PARA_ERR   E(1)   error in parallax [deg]",
         "SPECTYP    A(16)  spectral type"])

# OI_ARRAY definition (1st revision):
add_def("OI_ARRAY", 1,
        ["OI_REVN     I     revision number of the table definition",
         "ARRNAME     A     array name for cross-referencing",
         "FRAME       A     coordinate frame",
         "ARRAYX      D     array center X-coordinate [m]",
         "ARRAYY      D     array center Y-coordinate [m]",
         "ARRAYZ      D     array center Z-coordinate [m]",
         "------------------------------------------------------------",
         "TEL_NAME   A(16)  telescope name",
         "STA_NAME   A(16)  station name",
         "STA_INDEX  I(1)   station index",
         "DIAMETER   E(1)   element diameter [m]",
         "STAXYZ     D(3)   station coordinates relative to array center [m]"])

# OI_WAVELENGTH definition (1st revision):
add_def("OI_WAVELENGTH", 1,
        ["OI_REVN     I     revision number of the table definition",
         "INSNAME     A     name of detector for cross-referencing",
         "------------------------------------------------------------",
         "EFF_WAVE   E(1)   effective wavelength of channel [m]",
         "EFF_BAND   E(1)   effective bandpass of channel [m]"])

# OI_VIS definition (1st revision):
add_def("OI_VIS", 1,
        ["OI_REVN     I    revision number of the table definition",
         "DATE-OBS    A    UTC start date of observations",
         "ARRNAME    ?A    name of corresponding array",
         "INSNAME     A    name of corresponding detector",
         "------------------------------------------------------------",
         "TARGET_ID  I(1)  target number as index into OI_TARGET table",
         "TIME       D(1)  UTC time of observation [s]",
         "MJD        D(1)  modified Julian Day [day]",
         "INT_TIME   D(1)  integration time [s]",
         "VISAMP     D(W)  visibility amplitude",
         "VISAMPERR  D(W)  error in visibility amplitude",
         "VISPHI     D(W)  visibility phase [deg]",
         "VISPHIERR  D(W)  error in visibility phase [deg]",
         "UCOORD     D(1)  U coordinate of the data [m]",
         "VCOORD     D(1)  V coordinate of the data [m]",
         "STA_INDEX  I(2)  station numbers contributing to the data",
         "FLAG       L(W)  flag"])

# OI_VIS2 definition (1st revision):
add_def("OI_VIS2", 1,
        ["OI_REVN     I    revision number of the table definition",
         "DATE-OBS    A    UTC start date of observations",
         "ARRNAME    ?A    name of corresponding array",
         "INSNAME     A    name of corresponding detector",
         "------------------------------------------------------------",
         "TARGET_ID  I(1)  target number as index into OI_TARGET table",
         "TIME       D(1)  UTC time of observation [s]",
         "MJD        D(1)  modified Julian Day [day]",
         "INT_TIME   D(1)  integration time [s]",
         "VIS2DATA   D(W)  squared visibility",
         "VIS2ERR    D(W)  error in squared visibility",
         "UCOORD     D(1)  U coordinate of the data [m]",
         "VCOORD     D(1)  V coordinate of the data [m]",
         "STA_INDEX  I(2)  station numbers contributing to the data",
         "FLAG       L(W)  flag"])

# OI_T3 definition (1st revision):
add_def("OI_T3", 1,
        ["OI_REVN     I    revision number of the table definition",
         "DATE-OBS    A    UTC start date of observations",
         "ARRNAME    ?A    name of corresponding array",
         "INSNAME     A    name of corresponding detector",
         "------------------------------------------------------------",
         "TARGET_ID  I(1)  target number as index into OI_TARGET table",
         "TIME       D(1)  UTC time of observation [s]",
         "MJD        D(1)  modified Julian Day [day]",
         "INT_TIME   D(1)  integration time [s]",
         "T3AMP      D(W)  triple product amplitude",
         "T3AMPERR   D(W)  error in triple product amplitude",
         "T3PHI      D(W)  triple product phase [deg]",
         "T3PHIERR   D(W)  error in triple product phase [deg]",
         "U1COORD    D(1)  U coordinate of baseline AB of the triangle [m]",
         "V1COORD    D(1)  V coordinate of baseline AB of the triangle [m]",
         "U2COORD    D(1)  U coordinate of baseline BC of the triangle [m]",
         "V2COORD    D(1)  V coordinate of baseline BC of the triangle [m]",
         "STA_INDEX  I(3)  station numbers contributing to the data",
         "FLAG       L(W)  flag"])
