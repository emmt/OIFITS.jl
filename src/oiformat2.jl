#
# oiformat2.jl --
#
# Define 2nd revision of OI-FITS format.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015-2019: Éric Thiébaut, Jonathan Léger.
#
#------------------------------------------------------------------------------

# OI_TARGET definition (2nd revision):
add_def("OI_TARGET", 2,
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
         "SPECTYP    A(16)  spectral type",
         "CATEGORY  ?A(3)   CALibrator or SCIence target"])

# OI_ARRAY definition (2nd revision):
add_def("OI_ARRAY", 2,
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
         "STAXYZ     D(3)   station coordinates relative to array center [m]",
         "FOV        D(1)   photometric field of view [arcsec]",
         "FOVTYPE    A(6)   model for FOV: 'FWHM' or 'RADIUS'"])

# OI_WAVELENGTH definition (2nd revision):
add_def("OI_WAVELENGTH", 2,
        ["OI_REVN     I     revision number of the table definition",
         "INSNAME     A     name of detector for cross-referencing",
         "------------------------------------------------------------",
         "EFF_WAVE   E(1)   effective wavelength of channel [m]",
         "EFF_BAND   E(1)   effective bandpass of channel [m]"])

# OI_VIS2 definition (2nd revision):
add_def("OI_VIS2", 2,
        ["OI_REVN            I    revision number of the table definition",
         "DATE-OBS           A    UTC start date of observations",
         "ARRNAME            A    name of corresponding array",
         "INSNAME            A    name of corresponding detector",
         "CORRNAME          ?A    name of corresponding correlation table",
         "------------------------------------------------------------",
         "TARGET_ID          I(1)  target number as index into OI_TARGET table",
         "TIME               D(1)  UTC time of observation [s]",
         "MJD                D(1)  modified Julian Day [day]",
         "INT_TIME           D(1)  integration time [s]",
         "VIS2DATA           D(W)  squared visibility",
         "VIS2ERR            D(W)  error in squared visibility",
         "CORRINDX_VIS2DATA ?J(1)  index into correlation matrix for 1st VIS2DATA element",
         "UCOORD             D(1)  U coordinate of the data [m]",
         "VCOORD             D(1)  V coordinate of the data [m]",
         "STA_INDEX          I(2)  station numbers contributing to the data",
         "FLAG               L(W)  flag"])

# FIXME: In OI-FITS Rev. 2 only MJD and DATE-OBS must be used to express
#        time. The TIME column is retained only for backwards compatibility
#        and must contain zeros.

# OI_VIS definition (2nd revision):
add_def("OI_VIS", 2,
        ["OI_REVN          I       revision number of the table definition",
         "DATE-OBS         A       UTC start date of observations",
         "ARRNAME          A       name of corresponding array",
         "INSNAME          A       name of corresponding detector",
         "CORRNAME        ?A       name of corresponding correlation table",
         "AMPTYP          ?A       'ABSOLUTE', 'DIFFERENTIAL' or 'CORRELATED FLUX'",
         "PHITYP          ?A       'ABSOLUTE', or 'DIFFERENTIAL'",
         "AMPORDER        ?I       polynomial fit order for differential amplitudes",
         "PHIORDER        ?I       polynomial fit order for differential phases",
         "------------------------------------------------------------------------",
         "TARGET_ID        I(1)    target number as index into OI_TARGET table",
         "TIME             D(1)    UTC time of observation [s]",
         "MJD              D(1)    modified Julian Day [day]",
         "INT_TIME         D(1)    integration time [s]",
         "VISAMP           D(W)    visibility amplitude",
         "VISAMPERR        D(W)    error in visibility amplitude",
         "CORRINDX_VISAMP ?J(1)    index into correlation matrix for 1st VISAMP element",
         "VISPHI           D(W)    visibility phase [deg]",
         "VISPHIERR        D(W)    error in visibility phase [deg]",
         "CORRINDX_VISPHI ?J(1)    index into correlation matrix for 1st VISPHI element",
         "VISREFMAP       ?L(W,W)  true where spectral channels were taken as reference for differential visibility computation",
         "RVIS            ?D(W)    real part of complex coherent flux",
         "RVISERR         ?D(W)    error on RVIS",
         "CORRINDX_RVIS   ?J(1)    index into correlation matrix for 1st RVIS element",
         "IVIS            ?D(W)    imaginary part of complex coherent flux",
         "IVISERR         ?D(W)    error on IVIS",
         "CORRINDX_IVIS   ?J(1)    index into correlation matrix for 1st IVIS element",
         "UCOORD           D(1)    U coordinate of the data [m]",
         "VCOORD           D(1)    V coordinate of the data [m]",
         "STA_INDEX        I(2)    station numbers contributing to the data",
         "FLAG             L(W)    flag"])

# OI_T3 definition (2nd revision):
add_def("OI_T3", 2,
        ["OI_REVN         I     revision number of the table definition",
         "DATE-OBS        A     UTC start date of observations",
         "ARRNAME         A     name of corresponding array",
         "INSNAME         A     name of corresponding detector",
         "CORRNAME       ?A     name of corresponding correlation table",
         "------------------------------------------------------------",
         "TARGET_ID       I(1)  target number as index into OI_TARGET table",
         "TIME            D(1)  UTC time of observation [s]",
         "MJD             D(1)  modified Julian Day [day]",
         "INT_TIME        D(1)  integration time [s]",
         "T3AMP           D(W)  triple product amplitude",
         "T3AMPERR        D(W)  error in triple product amplitude",
         "CORRINDX_T3AMP ?J(1)  index into correlation matrix for 1st T3AMP element",
         "T3PHI           D(W)  triple product phase [deg]",
         "T3PHIERR        D(W)  error in triple product phase [deg]",
         "CORRINDX_T3PHI ?J(1)  index into correlation matrix for 1st T3PHI element",
         "U1COORD         D(1)  U coordinate of baseline AB of the triangle [m]",
         "V1COORD         D(1)  V coordinate of baseline AB of the triangle [m]",
         "U2COORD         D(1)  U coordinate of baseline BC of the triangle [m]",
         "V2COORD         D(1)  V coordinate of baseline BC of the triangle [m]",
         "STA_INDEX       I(3)  station numbers contributing to the data",
         "FLAG            L(W)  flag"])

# OI_SPECTRUM definition (1st revision):
add_def("OI_SPECTRUM", 1,
        ["OI_REVN            I     revision number of the table definition",
         "DATE-OBS           A     UTC start date of observations",
         "INSNAME            A     name of corresponding detector",
         "ARRNAME           ?A     name of corresponding array",
         "CORRNAME          ?A     name of corresponding correlation table",
         "FOV                D     area of sky over which flux is integrated [arcsec]",
         "FOVTYPE            A     model for FOV: 'FWHM' or 'RADIUS'",
         "CALSTAT            A     'C': spectrum is calibrated, 'U': uncalibrated",
         "------------------------------------------------------------",
         "TARGET_ID          I(1)  target number as index into OI_TARGET table",
         "MJD                D(1)  modified Julian Day [day]",
         "INT_TIME           D(1)  integration time [s]",
         "FLUXDATA           D(W)  flux",
         "FLUXERR            D(W)  flux error",
         "CORRINDX_FLUXDATA ?J(1)  index into correlation matrix for 1st FLUXDATA element",
         "STA_INDEX          I(1)  station number contributing to the data"])

# OI_CORR definition (1st revision):
add_def("OI_CORR", 1,
        ["OI_REVN            I     revision number of the table definition",
         "CORRNAME           A     name of correlation data set",
         "NDATA              I     number of correlated data",
         "------------------------------------------------------------",
         "IINDX              J(1)  1st index of correlation matrix element",
         "JINDX              J(1)  2nd index of correlation matrix element",
         "CORR               D(1)  matrix element"])

# OI_INSPOL definition (1st revision):
add_def("OI_INSPOL", 1,
        ["OI_REVN    I     revision number of the table definition",
         "DATE-OBS   A     UTC start date of observations",
         "NPOL       I     number of polarization types in this table",
         "ARRNAME    A     identifies corresponding OI_ARRAY",
         "ORIENT     A     orientation of the Jones Matrix: 'NORTH' (for on-sky orientation), or 'LABORATORY'",
         "MODEL      A     describe the way the Jones matrix is estimated",
         "---------------------------------------------------------------------------------------",
         "TARGET_ID  I(1)  target number as index into OI_TARGET table",
         "INSNAME    A(16) INSNAME of this polarization",
         "MJD_OBS    D(1)  modified Julian day, start of time lapse",
         "MJD_END    D(1)  modified Julian day, end of time lapse",
         "JXX        C(W)  complex Jones Matrix component along X axis",
         "JYY        C(W)  complex Jones Matrix component along Y axis",
         "JXY        C(W)  complex Jones Matrix component between Y and X axis",
         "JYX        C(W)  complex Jones Matrix component between Y and X axis",
         "STA_INDEX ?I(1)  station number for the above matrices"])
