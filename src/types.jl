#
# types.jl -
#
# Definitions of OI-FITS data types.
#
#-------------------------------------------------------------------------------------------

"""
    OIFITS.Undef

Type of `undef`, an alias for `UndefInitializer`.

"""
const Undef = typeof(undef)

"""
    OIFITS.unspecified

Singleton object (of type `OIFITS.Unspecified`) used to indicate unspecified
arguments/keywords in `OIFITS` package.

"""
struct Unspecified; end
const unspecified = Unspecified()
@doc @doc(Unspecified) unspecified

"""
    OIFITS.MissingColumn(col, ext)

Return an exception to be thrown when column `col` is not found in the table of a FITS
Table extension whose name is `ext`.

"""
struct MissingColumn <: Exception
    col::String
    ext::String
end

"""
    OIDataBlock

Abstract super-type of any OI-FITS data-block.

"""
abstract type OIDataBlock end

# Except for `OI_TARGET`, the data fields are specified as separate arrays, instead of
# arrays of structures.
#
# - Pros: (i) some fields may be absent (optional or not yet defined), (ii) easier to write
#   simple code, (iii) some fields are multi-dimensional, (iv) it is possible to check for
#   defined fields in Julia.
#
# - Cons: It is not certain that the sizes of the fields are compatible, but code to ensure
#   that may be provided.
#
# The sections called "Header Part" and "Data Part" correspond to the HDU FITS Table. The
# section called "Dependencies" consists in additional fields managed by the package to keep
# track of the data-set structure owning the data-block and to links to other data-blocks
# such as the instrument and the telescope array for complex visibility measurements.
#
# There is (to my knowledge) no official means to "undefine" a field in a structure so the
# OIFITS API does not provide means to pop a data-block out of a data-set. One has to build
# another data-set. Hopefully this is easy (with the API), fast, and memory efficient.

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

mutable struct OI_TARGET <: OIDataBlock
    # Header Part
    revn::Int
    # Data Part
    list::Vector{OITargetEntry}
    function OI_TARGET(list::AbstractVector{OITargetEntry}=OITargetEntry[];
                       revn::Integer=0)
        new(revn, list)
    end
end

mutable struct OI_ARRAY <: OIDataBlock
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
    OI_ARRAY(::Undef) = begin
        db = new()
        db.revn = 0
        db.arrayx = NaN
        db.arrayy = NaN
        db.arrayz = NaN
        return db
    end
end

mutable struct OI_WAVELENGTH <: OIDataBlock
    # Header Part
    revn::Int                # revision number of the table definition
    insname::String          # name of detector for cross-referencing
    # Data Part
    eff_wave::Vector{Float64}# effective wavelength of channel [m]
    eff_band::Vector{Float64}# effective bandpass of channel [m]
    OI_WAVELENGTH(::Undef) = begin
        db = new()
        db.revn = 0
        return db
    end
end

mutable struct OI_CORR <: OIDataBlock
    # Header Part
    revn::Int               # revision number of the table definition
    corrname::String        # name of correlation data-block
    ndata::Int              # number of correlated data
    # Data Part
    iindx::Vector{Int}      # 1st index of correlation matrix element
    jindx::Vector{Int}      # 2nd index of correlation matrix element
    corr::Vector{Float64}   # matrix element
    OI_CORR(::Undef) = begin
        db = new()
        db.revn = 0
        db.ndata = 0
        return db
    end
end

# FIXME: In OI-FITS Rev. 2 only MJD and DATE-OBS must be used to express time. The TIME
#        column is retained only for backwards compatibility and must contain zeros.

mutable struct OI_VIS <: OIDataBlock
    # Dependencies
    array::OI_ARRAY             # related telescope array
    instr::OI_WAVELENGTH        # related instrument wavelengths
    correl::OI_CORR             # related correlation data
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
    OI_VIS(::Undef) = begin
        db = new()
        db.revn = 0
        db.amporder = 0
        db.phiorder = 0
       return db
    end
end

mutable struct OI_VIS2 <: OIDataBlock
    # Dependencies
    array::OI_ARRAY               # related telescope array
    instr::OI_WAVELENGTH          # related instrument wavelengths
    correl::OI_CORR               # related correlation data
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
    OI_VIS2(::Undef) = begin
        db = new()
        db.revn = 0
        return db
    end
end

mutable struct OI_T3 <: OIDataBlock
    # Dependencies
    array::OI_ARRAY            # related telescope array
    instr::OI_WAVELENGTH       # related instrument wavelengths
    correl::OI_CORR            # related correlation data
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
    OI_T3(::Undef) = begin
        db = new()
        db.revn = 0
        return db
    end
end

mutable struct OI_FLUX <: OIDataBlock
    # Dependencies
    array::OI_ARRAY               # related telescope array
    instr::OI_WAVELENGTH          # related instrument wavelengths
    correl::OI_CORR               # related correlation data
    # Header Part
    revn::Int                     # revision number of the table definition
    date_obs::String              # UTC start date of observations
    insname::String               # name of corresponding detector
    arrname::String               # name of corresponding array
    corrname::String              # name of corresponding correlation table
    fov::Float64                  # area of sky over which flux is integrated [arcsec]
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
    OI_FLUX(::Undef) = begin
        db = new()
        db.revn = 0
        db.fov = NaN
        return db
    end
end

mutable struct OI_INSPOL <: OIDataBlock
    # Dependencies
    array::OI_ARRAY              # related telescope array
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
    OI_INSPOL(::Undef) = begin
        db = new()
        db.revn = 0
        db.npol = 0
        return db
    end
end

# Union of types of data-blocks that have a name (and are dependencies of
# others).
const NamedDataBlock = Union{OI_ARRAY,OI_WAVELENGTH,OI_CORR}

"""
    OIDDataSet

Type of objects storing the content of an OI-FITS file. All data-blocks containing
measurements (`OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, and `OI_INSPOL`) are stored into a
vector and thus indexed by an integer. Named data-blocks (`OI_ARRAY`, `OI_WAVELENGTH`, and
`OI_CORR`) are indexed by their names (converted to upper case letters, with leading and
trailing spaces stripped, multiple spaces replaced by a single ordinary space).

Reading an OI-FITS file is as simple as one of:

    data = OIDataSet(filename)
    data = read(OIDataSet, filename)

with `filename` the name of the file. Keyword `hack_revn` can be used to force the revision
number (FITS keyword `OI-REVN`) of all OI-FITS data-blocks; `hack_revn` may be set to an
integer, the revision number to assume for all data-blocks, or to a function that takes 2
arguments, the type and actual revision number of the current data-block, and that returns
the revision number to assume. For example, to force revision number 1 for all `OI_VIS`
data-blocks and left others unchanged:

    data = OIDataSet(filename; hack_revn = (T, revn) -> (T === OI_VIS ? 1 : revn))

Looping over `OI_VIS` data-blocks in the data-set can be done as follows:

    for db in data.vis
        ...
    end

and similarly for fields `vis2`, `t3`, `flux`, and `inspol` to access `OI_VIS2`, `OI_T3`,
`OI_FLUX`, and `OI_INSPOL` data-blocks.

The data-set has a number of other public properties:

    data.target           # the OI_TARGET data-block
    data.instr[insname]   # the OI_WAVELENGTH data-block named `insname`
    data.array[arrname]   # the OI_ARRAY data-block named `arrname`
    data.correl[corrname] # the OI_CORR data-block named `corrname`

"""
struct OIDataSet
    # Named data-blocks and their dictionaries to map names to indices. The target
    # identifiers `target_id` are rewritten to match Julia indexing in the vector of target
    # entries.
    target::OI_TARGET            # list of targets data-block
    target_dict::Dict{String,Int}# to map target name to index in `target` list
    target_id_map::Dict{Int,Int} # to re-write target identifiers
    array::Vector{OI_ARRAY}      # list of OI_ARRAY data-blocks
    array_dict::Dict{String,Int} # to map array names to index in `array` list
    instr::Vector{OI_WAVELENGTH} # list of OI_WAVELENGTH data-blocks
    instr_dict::Dict{String,Int} # to map instrument names to index in `instr` list
    correl::Vector{OI_CORR}      # list of OI_CORR data-blocks
    correl_dict::Dict{String,Int}# to map correlation names to index in `correl` list

    # Measurement data-blocks.
    vis::Vector{OI_VIS}          # list of OI_VIS data-blocks
    vis2::Vector{OI_VIS2}        # list of OI_VIS2 data-blocks
    t3::Vector{OI_T3}            # list of OI_T3 data-blocks
    flux::Vector{OI_FLUX}        # list of OI_FLUX data-blocks
    inspol::Vector{OI_INSPOL}    # list of OI_INSPOL data-blocks

    # The inner constructor creates an empty structure. Outer constructors are
    # provided to populate this structure with data-blocks.
    OIDataSet() = new(
        OI_TARGET(),     Dict{String,Int}(), Dict{Int,Int}(),
        OI_ARRAY[],      Dict{String,Int}(),
        OI_WAVELENGTH[], Dict{String,Int}(),
        OI_CORR[],       Dict{String,Int}(),
        OI_VIS[],
        OI_VIS2[],
        OI_T3[],
        OI_FLUX[],
        OI_INSPOL[])
end
