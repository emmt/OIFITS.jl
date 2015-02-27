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
# Copyright (C) 2015, Éric Thiébaut.
#
#------------------------------------------------------------------------------


# OI-FITS FORMAT DESCRIPTION TABLES
#
# The format of the OI-FITS data block is described by a vector of strings
# like:
#
#   ["KEYWORD FORMAT UNITS DESCR",
#     ...,
#     ...,
#    "---------------------------",
#    "COLUMN  FORMAT UNITS DESCR",
#     ...,
#     ...,
#     ...]
#
# where:
#
#   KEYWORD = keyword for HDU header.
#   COLUMN = column name for table (TTYPE).
#   FORMAT = nL where n is an integer and L a letter;
#       for keywords, 0 means optional and 1 means required;
#       for columns, a negative number means abs(n)*NWAVE.
#   UNITS = default units.
#   DESCR = description/comment.
#
# There may be any number of keyword definitions and any number of column
# definitions, the two parts are separated by a dash line like
# "--------------".
#

# OI_TARGET definition (1st revision):
add_def("OI_TARGET", 1,
        ["OI_REVN    1I -      revision number of the table definition",
         "------------------------------------------------------------",
         "TARGET_ID  1I -      index number",
         "TARGET    16A -      target name",
         "RAEP0      1D deg    RA at mean equinox",
         "DECEP0     1D deg    DEC at mean equinox",
         "EQUINOX    1E yr     equinox",
         "RA_ERR     1D deg    error in RA at mean equinox",
         "DEC_ERR    1D deg    error in DEC at mean equino",
         "SYSVEL     1D m/s    systemic radial velocity",
         "VELTYP     8A -      reference for radial velocity",
         "VELDEF     8A -      definition of radial velocity",
         "PMRA       1D deg/yr proper motion in RA",
         "PMDEC      1D deg/yr proper motion in DEC",
         "PMRA_ERR   1D deg/yr error of proper motion in RA",
         "PMDEC_ERR  1D deg/yr error of proper motion in DEC",
         "PARALLAX   1E deg    parallax",
         "PARA_ERR   1E deg    error in parallax",
         "SPECTYP   16A -      spectral type"])

# OI_ARRAY definition (1st revision):
add_def("OI_ARRAY", 1,
        ["OI_REVN    1I - revision number of the table definition",
         "ARRNAME    1A - array name for cross-referencing",
         "FRAME      1A - coordinate frame",
         "ARRAYX     1D m array center X-coordinate",
         "ARRAYY     1D m array center Y-coordinate",
         "ARRAYZ     1D m array center Z-coordinate",
         "------------------------------------------------------------",
         "TEL_NAME  16A - telescope name",
         "STA_NAME  16A - station name",
         "STA_INDEX  1I - station index",
         "DIAMETER   1E m element diameter",
         "STAXYZ     3D m station coordinates relative to array center"])

# OI_WAVELENGTH definition (1st revision):
add_def("OI_WAVELENGTH", 1,
        ["OI_REVN    1I - revision number of the table definition",
         "INSNAME    1A - name of detector for cross-referencing",
         "------------------------------------------------------------",
         "EFF_WAVE   1E m effective wavelength of channel",
         "EFF_BAND   1E m effective bandpass of channel"])

# OI_VIS definition (1st revision):
add_def("OI_VIS", 1,
        ["OI_REVN    1I -   revision number of the table definition",
         "DATE-OBS   1A -   UTC start date of observations",
         "ARRNAME    0A -   name of corresponding array",
         "INSNAME    1A -   name of corresponding detector",
         "------------------------------------------------------------",
         "TARGET_ID  1I -   target number as index into OI_TARGET table",
         "TIME       1D s   UTC time of observation",
         "MJD        1D day modified Julian Day",
         "INT_TIME   1D s   integration time",
         "VISAMP    -1D -   visibility amplitude",
         "VISAMPERR -1D -   error in visibility amplitude",
         "VISPHI    -1D deg visibility phase",
         "VISPHIERR -1D deg error in visibility phase",
         "UCOORD     1D m   U coordinate of the data",
         "VCOORD     1D m   V coordinate of the data",
         "STA_INDEX  2I -   station numbers contributing to the data",
         "FLAG      -1L -   flag"])

# OI_VIS2 definition (1st revision):
add_def("OI_VIS2", 1,
        ["OI_REVN    1I -   revision number of the table definition",
         "DATE-OBS   1A -   UTC start date of observations",
         "ARRNAME    0A -   name of corresponding array",
         "INSNAME    1A -   name of corresponding detector",
         "------------------------------------------------------------",
         "TARGET_ID  1I -   target number as index into OI_TARGET table",
         "TIME       1D s   UTC time of observation",
         "MJD        1D day modified Julian Day",
         "INT_TIME   1D s   integration time",
         "VIS2DATA  -1D -   squared visibility",
         "VIS2ERR   -1D -   error in squared visibility",
         "UCOORD     1D m   U coordinate of the data",
         "VCOORD     1D m   V coordinate of the data",
         "STA_INDEX  2I -   station numbers contributing to the data",
         "FLAG      -1L -   flag"])

# OI_T3 definition (1st revision):
add_def("OI_T3", 1,
        ["OI_REVN    1I -   revision number of the table definition",
         "DATE-OBS   1A -   UTC start date of observations",
         "ARRNAME    0A -   name of corresponding array",
         "INSNAME    1A -   name of corresponding detector",
         "------------------------------------------------------------",
         "TARGET_ID  1I -   target number as index into OI_TARGET table",
         "TIME       1D s   UTC time of observation",
         "MJD        1D day modified Julian Day",
         "INT_TIME   1D s   integration time",
         "T3AMP     -1D -   triple product amplitude",
         "T3AMPERR  -1D -   error in triple product amplitude",
         "T3PHI     -1D deg triple product phase",
         "T3PHIERR  -1D deg error in triple product phase",
         "U1COORD    1D m   U coordinate of baseline AB of the triangle",
         "V1COORD    1D m   V coordinate of baseline AB of the triangle",
         "U2COORD    1D m   U coordinate of baseline BC of the triangle",
         "V2COORD    1D m   V coordinate of baseline BC of the triangle",
         "STA_INDEX  3I -   station numbers contributing to the data",
         "FLAG      -1L -   flag"])

# OI_SPECTRUM definition (1st revision):
add_def("OI_SPECTRUM", 1,
        ["OI_REVN    1I -   revision number of the table definition",
         "DATE-OBS   1A -   UTC start date of observations",
         "INSNAME    1A -   name of corresponding detector",
        #"FOV        1D -   area of sky over which flux is integrated",
         "------------------------------------------------------------",
         "TARGET_ID  1I -   target number as index into OI_TARGET table",
         "MJD        1D day modified Julian Day",
         "INT_TIME   1D s   integration time",
         "FLUXDATA  -1D -   flux",
         "FLUXERR   -1D -   flux error"])

# Local Variables:
# mode: Julia
# tab-width: 8
# indent-tabs-mode: nil
# fill-column: 79
# coding: utf-8
# ispell-local-dictionary: "american"
# End:
