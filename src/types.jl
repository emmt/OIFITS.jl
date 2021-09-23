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

# Private singleton type used to indicate unspecified arguments/keywords.
struct Unspecified; end
const unspecified = Unspecified()

# Exception thrown when a keyword is not found in a FITS header.
struct MissingKeyword <: Exception
    key::String
    ext::String
end

# Exception thrown when a column is not found in a FITS table.
struct MissingColumn <: Exception
    col::String
    ext::String
end

# An empty string.
const empty_string = ""

"""

`OIDataBlock` is the abstract super-type of any OI-FITS data-block.

"""
abstract type OIDataBlock end

# Except for `OITarget`, the data fields are specified as separate arrays,
# instead of arrays of structures.
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
# a data-set.  One has to build another data-set.  Hopefully this is easy (with
# the API), fast and memory efficient.

struct OITargetEntry
    target_id ::Int     # index number
    target    ::String  # target name
    raep0     ::Float64 # RA at mean equinox [deg]
    decep0    ::Float64 # DEC at mean equinox [deg]
    equinox   ::Float64 # equinox [yr]
    ra_err    ::Float64 # error in RA at mean equinox [deg]
    dec_err   ::Float64 # error in DEC at mean equino [deg]
    sysvel    ::Float64 # systemic radial velocity [m/s]
    veltyp    ::String  # reference for radial velocity
    veldef    ::String  # definition of radial velocity
    pmra      ::Float64 # proper motion in RA [deg/yr]
    pmdec     ::Float64 # proper motion in DEC [deg/yr]
    pmra_err  ::Float64 # error of proper motion in RA [deg/yr]
    pmdec_err ::Float64 # error of proper motion in DEC [deg/yr]
    parallax  ::Float64 # parallax [deg]
    para_err  ::Float64 # error in parallax [deg]
    spectyp   ::String  # spectral type
    category  ::String  # "CAL"ibrator or "SCI"ence target
end

mutable struct OITarget <: OIDataBlock
    # Header Part
    revn::Int
    # Data Part
    rows::Vector{OITargetEntry}
    OITarget(::Undef) = new(0, OITargetEntry[])
    OITarget(; revn::Integer, rows::AbstractVector{OITargetEntry}) =
        new(revn, rows)
end

mutable struct OIArray <: OIDataBlock
    # Header Part
    revn::Int                # revision number of the table definition
    arrname::String          # array name for cross-referencing
    frame::String            # coordinate frame
    arrayx::Float64          # array center X-coordinate [m]
    arrayy::Float64          # array center Y-coordinate [m]
    arrayz::Float64          # array center Z-coordinate [m]
    # Data Part
    tel_name::Vector{String} # telescope name
    sta_name::Vector{String} # station name
    sta_index::Vector{Int}   # station index
    diameter::Vector{Float64}# element diameter [m]
    staxyz::Matrix{Float64}  # station coordinates relative to array center [m]
    fov::Vector{Float64}     # photometric field of view [arcsec]
    fovtype::Vector{String}  # model for FOV: "FWHM" or "RADIUS"
    OIArray(::Undef) = begin
        db = new()
        db.revn = 0
        db.arrayx = NaN
        db.arrayy = NaN
        db.arrayz = NaN
        return db
    end
end

mutable struct OIWavelength <: OIDataBlock
    # Header Part
    revn::Int                # revision number of the table definition
    insname::String          # name of detector for cross-referencing
    # Data Part
    eff_wave::Vector{Float64}# effective wavelength of channel [m]
    eff_band::Vector{Float64}# effective bandpass of channel [m]
    OIWavelength(::Undef) = begin
        db = new()
        db.revn = 0
        return db
    end
end

mutable struct OICorr <: OIDataBlock
    # Header Part
    revn::Int               # revision number of the table definition
    corrname::String        # name of correlation data-block
    ndata::Int              # number of correlated data
    # Data Part
    iindx::Vector{Int}      # 1st index of correlation matrix element
    jindx::Vector{Int}      # 2nd index of correlation matrix element
    corr::Vector{Float64}   # matrix element
    OICorr(::Undef) = begin
        db = new()
        db.revn = 0
        db.ndata = 0
        return db
    end
end

# FIXME: In OI-FITS Rev. 2 only MJD and DATE-OBS must be used to express
#        time. The TIME column is retained only for backwards compatibility
#        and must contain zeros.

mutable struct OIVis <: OIDataBlock
    # Dependencies
    array::OIArray              # related telescope array
    instr::OIWavelength         # related instrument wavelengths
    correl::OICorr              # related correlation data
    # Header Part
    revn::Int                   # revision number of the table definition
    date_obs::String            # UTC start date of observations
    arrname::String             # name of corresponding array
    insname::String             # name of corresponding detector
    corrname::String            # name of corresponding correlation table
    amptyp::String              # "ABSOLUTE", "DIFFERENTIAL" or "CORRELATED FLUX"
    phityp::String              # "ABSOLUTE", or "DIFFERENTIAL"
    amporder::Int               # polynomial fit order for differential amplitudes
    phiorder::Int               # polynomial fit order for differential phases
    # Data Part
    target_id::Vector{Int}      # target number as index into OI_TARGET table
    time::Vector{Float64}       # UTC time of observation [s]
    mjd::Vector{Float64}        # modified Julian Day [day]
    int_time::Vector{Float64}   # integration time [s]
    visamp::Matrix{Float64}     # visibility amplitude
    visamperr::Matrix{Float64}  # error in visibility amplitude
    corrindx_visamp::Vector{Int}# index into correlation matrix for 1st VISAMP element
    visphi::Matrix{Float64}     # visibility phase [deg]
    visphierr::Matrix{Float64}  # error in visibility phase [deg]
    corrindx_visphi::Vector{Int}# index into correlation matrix for 1st VISPHI element
    visrefmap::Array{Bool,3}    # true where spectral channels were taken as reference for differential visibility computation
    rvis::Matrix{Float64}       # real part of complex coherent flux
    rviserr::Matrix{Float64}    # error on RVIS
    corrindx_rvis::Vector{Int}  # index into correlation matrix for 1st RVIS element
    ivis::Matrix{Float64}       # imaginary part of complex coherent flux
    iviserr::Matrix{Float64}    # error on IVIS
    corrindx_ivis::Vector{Int}  # index into correlation matrix for 1st IVIS element
    ucoord::Vector{Float64}     # U coordinate of the data [m]
    vcoord::Vector{Float64}     # V coordinate of the data [m]
    sta_index::Matrix{Int}      # station numbers contributing to the data
    flag::Matrix{Bool}          # flags
    OIVis(::Undef) = begin
        db = new()
        db.revn = 0
        db.amporder = 0
        db.phiorder = 0
       return db
    end
end

mutable struct OIVis2 <: OIDataBlock
    # Dependencies
    array::OIArray                # related telescope array
    instr::OIWavelength           # related instrument wavelengths
    correl::OICorr                # related correlation data
    # Header Part
    revn::Int                     # revision number of the table definition
    date_obs::String              # UTC start date of observations
    arrname::String               # name of corresponding array
    insname::String               # name of corresponding detector
    corrname::String              # name of corresponding correlation table
    # Data Part
    target_id::Vector{Int}        # target number as index into OI_TARGET table
    time::Vector{Float64}         # UTC time of observation [s]
    mjd::Vector{Float64}          # modified Julian Day [day]
    int_time::Vector{Float64}     # integration time [s]
    vis2data::Matrix{Float64}     # squared visibility
    vis2err::Matrix{Float64}      # error in squared visibility
    corrindx_vis2data::Vector{Int}# index into correlation matrix for 1st VIS2DATA element
    ucoord::Vector{Float64}       # U coordinate of the data [m]
    vcoord::Vector{Float64}       # V coordinate of the data [m]
    sta_index::Matrix{Int}        # station numbers contributing to the data
    flag::Matrix{Bool}            # flags
    OIVis2(::Undef) = begin
        db = new()
        db.revn = 0
        return db
    end
end

mutable struct OIT3 <: OIDataBlock
    # Dependencies
    array::OIArray             # related telescope array
    instr::OIWavelength        # related instrument wavelengths
    correl::OICorr             # related correlation data
    # Header Part
    revn::Int                  # revision number of the table definition
    date_obs::String           # UTC start date of observations
    arrname::String            # name of corresponding array
    insname::String            # name of corresponding detector
    corrname::String           # name of corresponding correlation table
    # Data Part
    target_id::Vector{Int}     # target number as index into OI_TARGET table
    time::Vector{Float64}      # UTC time of observation [s]
    mjd::Vector{Float64}       # modified Julian Day [day]
    int_time::Vector{Float64}  # integration time [s]
    t3amp::Matrix{Float64}     # triple product amplitude
    t3amperr::Matrix{Float64}  # error in triple product amplitude
    corrindx_t3amp::Vector{Int}# index into correlation matrix for 1st T3AMP element
    t3phi::Matrix{Float64}     # triple product phase [deg]
    t3phierr::Matrix{Float64}  # error in triple product phase [deg]
    corrindx_t3phi::Vector{Int}# index into correlation matrix for 1st T3PHI element
    u1coord::Vector{Float64}   # U coordinate of baseline AB of the triangle [m]
    v1coord::Vector{Float64}   # V coordinate of baseline AB of the triangle [m]
    u2coord::Vector{Float64}   # U coordinate of baseline BC of the triangle [m]
    v2coord::Vector{Float64}   # V coordinate of baseline BC of the triangle [m]
    sta_index::Matrix{Int}     # station numbers contributing to the data
    flag::Matrix{Bool}         # flags
    OIT3(::Undef) = begin
        db = new()
        db.revn = 0
        return db
    end
end

mutable struct OIFlux <: OIDataBlock
    # Dependencies
    array::OIArray                # related telescope array
    instr::OIWavelength           # related instrument wavelengths
    correl::OICorr                # related correlation data
    # Header Part
    revn::Int                     # revision number of the table definition
    date_obs::String              # UTC start date of observations
    insname::String               # name of corresponding detector
    arrname::String               # name of corresponding array
    corrname::String              # name of corresponding correlation table
    fov::Cdouble                  # area of sky over which flux is integrated [arcsec]
    fovtype::String               # model for FOV: "FWHM" or "RADIUS"
    calstat::String               # "C": spectrum is calibrated, "U": uncalibrated
    # Data Part
    target_id::Vector{Int}        # target number as index into OI_TARGET table
    mjd::Vector{Float64}          # modified Julian Day [day]
    int_time::Vector{Float64}     # integration time [s]
    fluxdata::Matrix{Float64}     # flux
    fluxerr::Matrix{Float64}      # flux error
    corrindx_fluxdata::Vector{Int}# index into correlation matrix for 1st FLUXDATA element
    sta_index::Vector{Int}        # station number contributing to the data
    flag::Matrix{Bool}            # flags
    OIFlux(::Undef) = begin
        db = new()
        db.revn = 0
        db.fov = NaN
        return db
    end
end

mutable struct OIInsPol <: OIDataBlock
    # Dependencies
    array::OIArray               # related telescope array
    instr::OIWavelength          # related instrument wavelengths
    # Header Part
    revn::Int                    # revision number of the table definition
    date_obs::String             # UTC start date of observations
    npol::Int                    # number of polarization types in this table
    arrname::String              # identifies corresponding OI_ARRAY
    orient::String               # orientation of the Jones Matrix: "NORTH" (for on-sky orientation), or "LABORATORY"
    model::String                # describe the way the Jones matrix is estimated
    # Data Part
    target_id::Vector{Int}       # target number as index into OI_TARGET table
    insname::Vector{String}      # INSNAME of this polarization
    mjd_obs::Vector{Float64}     # modified Julian day, start of time lapse
    mjd_end::Vector{Float64}     # modified Julian day, end of time lapse
    jxx::Matrix{Complex{Float64}}# complex Jones Matrix component along X axis
    jyy::Matrix{Complex{Float64}}# complex Jones Matrix component along Y axis
    jxy::Matrix{Complex{Float64}}# complex Jones Matrix component between Y and X axis
    jyx::Matrix{Complex{Float64}}# complex Jones Matrix component between Y and X axis
    sta_index::Vector{Int}       # station number for the above matrices
    OIInsPol(::Undef) = begin
        db = new()
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
mutable struct OIData
    target::OITarget           # OI_TARGET data-block (unique)
    array::Vector{OIArray}     # OI_ARRAY data-blocks
    instr::Vector{OIWavelength}# OI_WAVELENGTH data-blocks
    correl::Vector{OICorr}     # OI_CORR data-blocks
    vis::Vector{OIVis}         # OI_VIS data-blocks
    vis2::Vector{OIVis2}       # OI_VIS2 data-blocks
    t3::Vector{OIT3}           # OI_T3 data-blocks
    flux::Vector{OIFlux}       # OI_FLUX data-blocks
    inspol::Vector{OIInsPol}   # OI_INSPOL data-blocks

    # The inner constructor creates an empty structure.  Outer constructors are
    # provided to populate this structure with data-blocks.
    OIData(::Undef) = begin
        dat = new()
        dat.array  = OIArray[]
        dat.instr  = OIWavelength[]
        dat.correl = OICorr[]
        dat.vis    = OIVis[]
        dat.vis2   = OIVis2[]
        dat.t3     = OIT3[]
        dat.flux   = OIFlux[]
        dat.inspol = OIInsPol[]
        return dat
    end
end
