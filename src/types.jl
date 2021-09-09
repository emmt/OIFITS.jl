#
# types.jl -
#
# Definitions of OI-FITS data types.
#
#------------------------------------------------------------------------------

"""

`OIFITS.Undef` is the type of `undef`.

"""
const Undef = typeof(undef)

"""

`OIDataBlock{T}` is the abstract super-type of any OI data-block.  Parameter
`T` is the floating-point type of the stored data.

"""
abstract type OIDataBlock{T<:AbstractFloat} end

# The data fields are specified as separate arrays, instead of arrays of
# structures.
#
# - Pros: (i) some fields may be absent (optional or not yet defined), (ii)
#   easier to write simple code, (iii) some fields are multi-dimensional, (iv)
#   it is possible to check for defined fields in Julia.
#
# - Cons: It is not certain that the sizes of the fields are compatible, but
#   code to ensure that may be provided.
#
# The sections called "Header Part" and "Data Part" correspond to the HDU FITS
# Table.  The section called "Dependencies" consists in additional fields
# managed by the package to keep track of the data-set structure owning the
# data-block and to links to other data-blocks such as the instrument and the
# telescope array for complex visibility measurements.
#
# There is (to my knowledge) no official means to "undefine" a field in a
# structure so the OIFITS API does not provide means to pop a data-block out of
# a data-set.  One has to rebuild a data-set.  Hopefully this is easy (with the
# API), fast and memory efficient.
#
mutable struct OITarget{T<:AbstractFloat} <: OIDataBlock{T}
    # Header Part
    revn::Int               # revision number of the table definition
    # Data Part
    target_id::Vector{Int}  # index number
    target::Vector{String}  # target name
    raep0::Vector{T}        # RA at mean equinox [deg]
    decep0::Vector{T}       # DEC at mean equinox [deg]
    equinox::Vector{T}      # equinox [yr]
    ra_err::Vector{T}       # error in RA at mean equinox [deg]
    dec_err::Vector{T}      # error in DEC at mean equino [deg]
    sysvel::Vector{T}       # systemic radial velocity [m/s]
    veltyp::Vector{String}  # reference for radial velocity
    veldef::Vector{String}  # definition of radial velocity
    pmra::Vector{T}         # proper motion in RA [deg/yr]
    pmdec::Vector{T}        # proper motion in DEC [deg/yr]
    pmra_err::Vector{T}     # error of proper motion in RA [deg/yr]
    pmdec_err::Vector{T}    # error of proper motion in DEC [deg/yr]
    parallax::Vector{T}     # parallax [deg]
    para_err::Vector{T}     # error in parallax [deg]
    spectyp::Vector{String} # spectral type
    category::Vector{String}# "CAL"ibrator or "SCI"ence target
    OITarget{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        return db
    end
end

mutable struct OIArray{T<:AbstractFloat} <: OIDataBlock{T}
    # Header Part
    revn::Int               # revision number of the table definition
    arrname::String         # array name for cross-referencing
    frame::String           # coordinate frame
    arrayx::Cdouble         # array center X-coordinate [m]
    arrayy::Cdouble         # array center Y-coordinate [m]
    arrayz::Cdouble         # array center Z-coordinate [m]
    # Data Part
    tel_name::Vector{String}# telescope name
    sta_name::Vector{String}# station name
    sta_index::Vector{Int}  # station index
    diameter::Vector{T}     # element diameter [m]
    staxyz::Matrix{T}       # station coordinates relative to array center [m]
    fov::Vector{T}          # photometric field of view [arcsec]
    fovtype::Vector{String} # model for FOV: "FWHM" or "RADIUS"
    OIArray{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        db.arrayx = NaN
        db.arrayy = NaN
        db.arrayz = NaN
        return db
    end
end

mutable struct OIWavelength{T<:AbstractFloat} <: OIDataBlock{T}
    # Header Part
    revn::Int               # revision number of the table definition
    insname::String         # name of detector for cross-referencing
    # Data Part
    eff_wave::Vector{T}     # effective wavelength of channel [m]
    eff_band::Vector{T}     # effective bandpass of channel [m]
    OIWavelength{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        return db
    end
end

mutable struct OICorr{T<:AbstractFloat} <: OIDataBlock{T}
    # Header Part
    revn::Int               # revision number of the table definition
    corrname::String        # name of correlation data-block
    ndata::Int              # number of correlated data
    # Data Part
    iindx::Vector{Int}      # 1st index of correlation matrix element
    jindx::Vector{Int}      # 2nd index of correlation matrix element
    corr::Vector{T}         # matrix element
    OICorr{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        db.ndata = 0
        return db
    end
end

# FIXME: In OI-FITS Rev. 2 only MJD and DATE-OBS must be used to express
#        time. The TIME column is retained only for backwards compatibility
#        and must contain zeros.

mutable struct OIVis{T<:AbstractFloat} <: OIDataBlock{T}
    # Dependencies
    array::OIArray{T}       # related telescope array
    instr::OIWavelength{T}  # related instrument wavelengths
    correl::OICorr{T}       # related correlation data
    # Header Part
    revn::Int               # revision number of the table definition
    date_obs::String        # UTC start date of observations
    arrname::String         # name of corresponding array
    insname::String         # name of corresponding detector
    corrname::String        # name of corresponding correlation table
    amptyp::String          # "ABSOLUTE", "DIFFERENTIAL" or "CORRELATED FLUX"
    phityp::String          # "ABSOLUTE", or "DIFFERENTIAL"
    amporder::Int           # polynomial fit order for differential amplitudes
    phiorder::Int           # polynomial fit order for differential phases
    # Data Part
    target_id::Vector{Int}  # target number as index into OI_TARGET table
    time::Vector{T}         # UTC time of observation [s]
    mjd::Vector{T}          # modified Julian Day [day]
    int_time::Vector{T}     # integration time [s]
    visamp::Matrix{T}       # visibility amplitude
    visamperr::Matrix{T}    # error in visibility amplitude
    corrindx_visamp::Vector{Int} # index into correlation matrix for 1st VISAMP element
    visphi::Matrix{T}       # visibility phase [deg]
    visphierr::Matrix{T}    # error in visibility phase [deg]
    corrindx_visphi::Vector{Int} # index into correlation matrix for 1st VISPHI element
    visrefmap::Array{Bool,3}# true where spectral channels were taken as reference for differential visibility computation
    rvis::Matrix{T}         # real part of complex coherent flux
    rviserr::Matrix{T}      # error on RVIS
    corrindx_rvis::Vector{Int} # index into correlation matrix for 1st RVIS element
    ivis::Matrix{T}         # imaginary part of complex coherent flux
    iviserr::Matrix{T}      # error on IVIS
    corrindx_ivis::Vector{Int} # index into correlation matrix for 1st IVIS element
    ucoord::Vector{T}       # U coordinate of the data [m]
    vcoord::Vector{T}       # V coordinate of the data [m]
    sta_index::Matrix{Int}  # station numbers contributing to the data
    flag::Matrix{Bool}      # flags
    OIVis{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        db.amporder = 0
        db.phiorder = 0
       return db
    end
end

mutable struct OIVis2{T} <: OIDataBlock{T}
    # Dependencies
    array::OIArray{T}       # related telescope array
    instr::OIWavelength{T}  # related instrument wavelengths
    correl::OICorr{T}       # related correlation data
    # Header Part
    revn::Int               # revision number of the table definition
    date_obs::String        # UTC start date of observations
    arrname::String         # name of corresponding array
    insname::String         # name of corresponding detector
    corrname::String        # name of corresponding correlation table
    # Data Part
    target_id::Vector{Int}  # target number as index into OI_TARGET table
    time::Vector{T}         # UTC time of observation [s]
    mjd::Vector{T}          # modified Julian Day [day]
    int_time::Vector{T}     # integration time [s]
    vis2data::Matrix{T}     # squared visibility
    vis2err::Matrix{T}      # error in squared visibility
    corrindx_vis2data::Vector{Int} # index into correlation matrix for 1st VIS2DATA element
    ucoord::Vector{T}       # U coordinate of the data [m]
    vcoord::Vector{T}       # V coordinate of the data [m]
    sta_index::Matrix{Int}  # station numbers contributing to the data
    flag::Matrix{Bool}      # flags
    OIVis2{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        return db
    end
end

mutable struct OIT3{T<:AbstractFloat} <: OIDataBlock{T}
    # Dependencies
    array::OIArray{T}       # related telescope array
    instr::OIWavelength{T}  # related instrument wavelengths
    correl::OICorr{T}       # related correlation data
    # Header Part
    revn::Int               # revision number of the table definition
    date_obs::String        # UTC start date of observations
    arrname::String         # name of corresponding array
    insname::String         # name of corresponding detector
    corrname::String        # name of corresponding correlation table
    # Data Part
    target_id::Vector{Int}  # target number as index into OI_TARGET table
    time::Vector{T}         # UTC time of observation [s]
    mjd::Vector{T}          # modified Julian Day [day]
    int_time::Vector{T}     # integration time [s]
    t3amp::Matrix{T}        # triple product amplitude
    t3amperr::Matrix{T}     # error in triple product amplitude
    corrindx_t3amp::Vector{Int} # index into correlation matrix for 1st T3AMP element
    t3phi::Matrix{T}        # triple product phase [deg]
    t3phierr::Matrix{T}     # error in triple product phase [deg]
    corrindx_t3phi::Vector{Int} # index into correlation matrix for 1st T3PHI element
    u1coord::Vector{T}      # U coordinate of baseline AB of the triangle [m]
    v1coord::Vector{T}      # V coordinate of baseline AB of the triangle [m]
    u2coord::Vector{T}      # U coordinate of baseline BC of the triangle [m]
    v2coord::Vector{T}      # V coordinate of baseline BC of the triangle [m]
    sta_index::Matrix{Int}  # station numbers contributing to the data
    flag::Matrix{Bool}      # flags
    OIT3{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        return db
    end
end

mutable struct OIFlux{T<:AbstractFloat} <: OIDataBlock{T}
    # Dependencies
    array::OIArray{T}       # related telescope array
    instr::OIWavelength{T}  # related instrument wavelengths
    correl::OICorr{T}       # related correlation data
    # Header Part
    revn::Int               # revision number of the table definition
    date_obs::String        # UTC start date of observations
    insname::String         # name of corresponding detector
    arrname::String         # name of corresponding array
    corrname::String        # name of corresponding correlation table
    fov::Cdouble            # area of sky over which flux is integrated [arcsec]
    fovtype::String         # model for FOV: "FWHM" or "RADIUS"
    calstat::String         # "C": spectrum is calibrated, "U": uncalibrated
    # Data Part
    target_id::Vector{Int}  # target number as index into OI_TARGET table
    mjd::Vector{T}          # modified Julian Day [day]
    int_time::Vector{T}     # integration time [s]
    fluxdata::Matrix{T}     # flux
    fluxerr::Matrix{T}      # flux error
    corrindx_fluxdata::Vector{Int}# index into correlation matrix for 1st FLUXDATA element
    sta_index::Vector{Int}  # station number contributing to the data
    flag::Matrix{Bool}      # flags
    OIFlux{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        db.fov = NaN
        return db
    end
end

mutable struct OIInsPol{T<:AbstractFloat} <: OIDataBlock{T}
    # Dependencies
    array::OIArray{T}       # related telescope array
    instr::OIWavelength{T}  # related instrument wavelengths
    # Header Part
    revn::Int               # revision number of the table definition
    date_obs::String        # UTC start date of observations
    npol::Int               # number of polarization types in this table
    arrname::String         # identifies corresponding OI_ARRAY
    orient::String          # orientation of the Jones Matrix: "NORTH" (for on-sky orientation), or "LABORATORY"
    model::String           # describe the way the Jones matrix is estimated
    # Data Part
    target_id::Vector{Int}  # target number as index into OI_TARGET table
    insname::Vector{String} # INSNAME of this polarization
    mjd_obs::Vector{T}      # modified Julian day, start of time lapse
    mjd_end::Vector{T}      # modified Julian day, end of time lapse
    jxx::Matrix{Complex{T}} # complex Jones Matrix component along X axis
    jyy::Matrix{Complex{T}} # complex Jones Matrix component along Y axis
    jxy::Matrix{Complex{T}} # complex Jones Matrix component between Y and X axis
    jyx::Matrix{Complex{T}} # complex Jones Matrix component between Y and X axis
    sta_index::Vector{Int}  # station number for the above matrices
    OIInsPol{T}() where {T<:AbstractFloat} = begin
        db = new{T}()
        db.revn = 0
        db.npol = 0
        return db
    end
end

# Union of datablocks with dependencies.
const DataBlocksWithDependencies = Union{OIVis,OIVis2,OIT3,OIFlux,OIInsPol}

"""

`OIData` stores the contents of an OI-FITS file.  All data-blocks containing
measurements (OI_VIS, OI_VIS2, OI_T3, OI_FLUX and OI_INSPOL) are stored into a
vector and thus indexed by an integer.  Named data-blocks (OI_ARRAY,
OI_WAVELENGTH and OI_CORR) are indexed by their names (converted to upper case
letters, with leading and trailing spaces stripped, multiple spaces replaced by
a single ordinary space).

    # Loop over OIVis instances in data-set.
    for db in data.vis
        ...
    end
    data.target           # yields an OITarget instance
    data.instr[insname]   # yields an OIWavelength instance
    data.array[arrname]   # yields an OIArray instance
    data.correl[corrname] # yields an OICorr instance

"""
mutable struct OIData{T<:AbstractFloat}
    target::OITarget{T}           # OI_TARGET data-block (unique)
    array::Vector{OIArray{T}}     # OI_ARRAY data-blocks
    instr::Vector{OIWavelength{T}}# OI_WAVELENGTH data-blocks
    correl::Vector{OICorr{T}}     # OI_CORR data-blocks
    vis::Vector{OIVis{T}}         # OI_VIS data-blocks
    vis2::Vector{OIVis2{T}}       # OI_VIS2 data-blocks
    t3::Vector{OIT3{T}}           # OI_T3 data-blocks
    flux::Vector{OIFlux{T}}       # OI_FLUX data-blocks
    inspol::Vector{OIInsPol{T}}   # OI_INSPOL data-blocks

    # The inner constructor creates an empty structure.  Outer constructors are
    # provided to populate this structure with data-blocks.
    OIData{T}() where {T<:AbstractFloat} = begin
        dat = new{T}()
        dat.array  = OIArray{T}[]
        dat.instr  = OIWavelength{T}[]
        dat.correl = OICorr{T}[]
        dat.vis    = OIVis{T}[]
        dat.vis2   = OIVis2{T}[]
        dat.t3     = OIT3{T}[]
        dat.flux   = OIFlux{T}[]
        dat.inspol = OIInsPol{T}[]
        return dat
    end
end
