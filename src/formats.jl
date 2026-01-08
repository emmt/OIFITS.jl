#
# formats.jl -
#
# Definitions of OI-FITS tables.
#
#-------------------------------------------------------------------------------------------

"""
    OIFITS.Formats

This module provides definitions of OI-FITS tables.

"""
module Formats

"""
    def = OIFITS.Formats.FieldDefinition(...)

Return the definition of a keyword or column of a OI-FITS data-block. The following
properties are available:

- `def.name` specifies the name of the keyword or column in the FITS file.

- `def.symb` specifies the symbolic field name in the corresponding Julia structure.

- `def.type` specifies the data type. It is `:A` for strings, `:C` for complexes, `:D` or
  `:E` for reals (floating-points), `:I` or `:J` for integers, and `:L` for booleans.

- `def.rank` is `0` for keywords, `n > 0` for `n` elements per row, `n < 0` for cell of `-n`
  dimensions of length equal to the number of spectral channels.

- `def.optional` specifies whether the keyword / column is optional according to the
  standard.

- `def.descr` is a string describing the keyword / column.

These instances are built by the macros [`OIFITS.Formats.@header`](@ref) and
[`OIFITS.Formats.@column`](@ref) for a FITS keyword and column respectively.

Base method `ndims(def)` yields the dimensionality of `def`, that is `0` for a FITS keyword,
`1` for a column of scalar cells and `n ≥ 2` for a column of wavelength-wise
multi-dimensional cells.

"""
struct FieldDefinition
    name::String   # Keyword/column name in FITS file.
    symb::Symbol   # Keyword/column symbolic name in structures.
    type::Symbol   # Data type.
    rank::Int      # 0 for keywords, n > 0 for n elements per row, n < 0 for
                   # cell of -n dimensions of length equal to the number of
                   # spectral channels.
    optional::Bool # Optional field?
    descr::String  # Description.
    units::String  # Units.
end

Base.ndims(s::FieldDefinition) =
    s.rank == 0   ? 0          : # keyword
    s.type === :A ? 1          : # column of strings
    s.rank == 1   ? 1          : # column of values
    s.rank < 0    ? 1 - s.rank : # column of wavelength-wise multi-dim. cells
    2                            # column of multiple values

const FORMATS = Dict{Tuple{Symbol,Int},Vector{FieldDefinition}}()

"""
    OIFITS.get_format(ext, revn; throw_errors=false) -> A

Return a vector of [`OIFITS.Formats.Field`](@ref) instances corresponding to the FITS
keywords and columns of an OI-FITS data-block `ext` with revision number `revn`. If `ext`
and `revn` do not correspond to any known definition, `nothing` is returned if
`throw_errors` is `false` (the default), otherwise an exception is thrown.

For example:

    for def in OIFITS.get_format(:OI_VIS2, 2)
        println(
            (def.rank == 0 ? "Keyword" : "Column "), " \"",
            def.name, "\" => `", def.symb, "` -- ", def.descr)
    end

"""
function get_format(ext::Symbol, rev::Integer; throw_errors::Bool=false)
    spec = get(FORMATS, (ext,Int(rev)), nothing)
    if spec === nothing && throw_errors
        error("unknown revision ", rev, " of OI-FITS extension \"", ext, "\"")
    end
    return spec
end

"""
    @header key type descr

Expand as an instance of [`OIFITS.Formats.FieldDefinition`](@ref) for an OI-FITS header
keyword `key` (a string) with value type and description specified by `type` and `descr`.
The syntax closely follows that of the tables of the papers specifying the OI-FITS format.
Arguments `key` and `descr` are strings, `type` is one of the letter `I` (for integer), `D`
(for floating-point), or `A` (for string). Optional keywords have their `type` suffixed by a
single quote (as the `ARRNAME` keyword in the example below).

# Examples

```julia
@header "OI_REVN"   I   "revision number of the table definition"
@header "DATE-OBS"  A   "UTC start date of observations"
@header "ARRNAME"   A'  "name of corresponding array"
@header "INSNAME"   A   "name of corresponding detector"
@header "ARRAYZ"    D   "array center Z-coordinate [m]"
```

# References

- G. Duvert, J. S. Young & Ch. A. Hummel: "*OIFITS 2: the 2nd version of the data exchange
  standard for optical interferometry*", Astron. & Astrophys., **597**, A8 (2017).

- T. A. Pauls, J. S. Young, W. D. Cotton & J. D. Monnier: "*A Data Exchange Standard for
  Optical (Visible/IR) Interferometry*", PASP **117**, 1255-1262 (2005).

"""
macro header(key, type, str)
    opt, sym, rank = parse_type(type)
    isa(key, String) || error("expecting a string for keyword name, got ", key)
    rank == 0 || error("invalid dimensions in `", type, "` for header keyword")
    descr, units = parse_description(str)
    quote
        FieldDefinition($key, $(QuoteNode(to_field_name(key))),
                        $(QuoteNode(sym)), $rank, $opt, $descr, $units)
    end
end

"""
    @column key type descr

Expand as an instance of [`OIFITS.Formats.FieldDefinition`](@ref) for an OI-FITS column
`key` (a string or symbol) with cell type and dimensions specified by `type` and description
given by `descr`. The syntax closely follows that of the tables of the papers specifying the
OI-FITS format. Arguments `key` and `descr` are strings, `type` is one of the letters `L`
(for logical), `I`, `J`, `K` (for 16-,32-,64-bit integers), `E`, `D` (for 32-,64-bit
floating-point), or `A` (for string), followed by the cell dimensions in parenthesis. A
wavelength-wise is indicated by the letter `W`. Optional columns have their `type` suffixed
by single quote.

# Examples

```julia
@column "VISPHI"           D(W)     "visibility phase [deg]"
@column "VISPHIERR"        D(W)     "error in visibility phase [deg]"
@column "CORRINDX_VISPHI"  J(1)'    "index into correlation matrix for 1st VISPHI element"
@column "VISREFMAP"        L(W,W)'  "true where spectral channels were taken as reference for differential visibility computation"
@column "TARGET_ID"        I(1)     "target number as index into OI_TARGET table"
@column "TIME"             D(1)     "UTC time of observation [s]"
```

# References

- G. Duvert, J. S. Young & Ch. A. Hummel: "*OIFITS 2: the 2nd version of the data exchange
  standard for optical interferometry*", Astron. & Astrophys., **597**, A8 (2017).

- T. A. Pauls, J. S. Young, W. D. Cotton & J. D. Monnier: "*A Data Exchange Standard for
  Optical (Visible/IR) Interferometry*", PASP **117**, 1255-1262 (2005).

"""
macro column(key, type, str)
    opt, sym, rank = parse_type(type)
    isa(key, String) || error("expecting a string for column name, got ", key)
    rank != 0 || error("invalid dimensions in `", type, "` for column")
    descr, units = parse_description(str)
    quote
        FieldDefinition($key, $(QuoteNode(to_field_name(key))),
                        $(QuoteNode(sym)), $rank, $opt, $descr, $units)
    end
end

to_field_name(str::AbstractString) =
    Symbol(lowercase(replace(replace(str, r"^OI_" => ""), '-' => '_')))

"""
    parse_description(str) -> (descr, units)

Return a 2-tuple of strings extracted from string `str` which shall have the form:

    str = "\$descr [\$units]"

with any leading/trailing spaces and where the units part is optional (returned as an empty
string if not specified).

"""
function parse_description(str::String)
    txt = strip(str)
    if endswith(txt, ']')
        i = findlast('[', txt)
        if i !== nothing
            units = String(txt[nextind(txt,i):prevind(txt,end)])
            i_first = firstindex(txt)
            if i > i_first
                return (String(strip(txt[i_first:prevind(txt,i)])), units)
            else
                return ("", units)
            end
        end
    end
    return (String(txt), "")
end

parse_description(x::Any) = error("description must be a string")

"""
     parse_type(x) -> (opt, sym, rank)

Parse type `x` in OI-FITS format specification and yields a 3-tuple `(opt,type,dims)` where
`opt` is a boolean indicating whether the field is optional, `sym` is a single letter symbol
corresponding to the type of the field and `rank` specifies the dimensions.

Argument `x` is a symbol or an expression.

""" parse_type

const APOSTROPHE = Symbol(Char(39))

# 1st pass is to figure out the optional flag (an apostrophe).
parse_type(x::Symbol) = parse_type(false, x)
parse_type(x::Expr) = begin
    if x.head === APOSTROPHE
        parse_type(true, x.args[1])
    else
        parse_type(false, x)
    end
end

# 2nd pass is to parse type letter and dimensions.
parse_type(opt::Bool, x::Symbol) = (opt, parse_type_letter(x), 0)
parse_type(opt::Bool, x::Expr) = begin
    if x.head === :call && length(x.args) ≥ 1
        symb = parse_type_letter(x.args[1])
        rank = parse_type_dims(x.args[2:end]...)
        return (opt, symb, rank)
    else
        invalid_type(opt, x)
    end
end
parse_type(opt::Bool, x::Any) = invalid_type(opt, x)

parse_type_letter(x::Symbol) = (x ∈ (:A, :C, :D, :E, :I, :J, :L) ? x :
                                error("unknown type letter/symbol: ", x))

parse_type_dims(dims...) = invalid_type_dims(dims)
parse_type_dims(dim::Integer) =
    (dim ≥ 1 ? Int(dim) : invalid_type_dims((dim,)))
parse_type_dims(dim::Symbol) =
    (dim === :W ? -1 : invalid_type_dims((dim,)))
parse_type_dims(dim1::Symbol, dim2::Symbol) =
    (dim1 === :W && dim2 === :W ? -2 : invalid_type_dims((dim1,dim2)))

@noinline invalid_type_dims(dims) = error("invalid type dimensions: ", dims)

@noinline invalid_type(opt::Bool, x) =
    error("invalid type specification: ", x, (opt ? APOSTROPHE : nothing))


function define(db::AbstractString,
                numbs::Tuple{Vararg{Integer}},
                specs::AbstractVector{FieldDefinition})
    for revn in numbs
        define(Symbol(db), revn, specs)
    end
end

define(db::AbstractString, revn::Integer, specs::AbstractVector{FieldDefinition}) =
    define(Symbol(db), revn, specs)

function define(db::Symbol, revn::Integer, specs::AbstractVector{FieldDefinition})
    FORMATS[(db,Int(revn))] = specs
end

# OI_TARGET definition (1st revision):
define("OI_TARGET", 1, [
    @header "OI_REVN"    I      "revision number of the table definition"
    #---------------------------------------------------------
    @column "TARGET_ID"  I(1)   "index number"
    @column "TARGET"     A(16)  "target name"
    @column "RAEP0"      D(1)   "RA at mean equinox [deg]"
    @column "DECEP0"     D(1)   "DEC at mean equinox [deg]"
    @column "EQUINOX"    E(1)   "equinox [yr]"
    @column "RA_ERR"     D(1)   "error in RA at mean equinox [deg]"
    @column "DEC_ERR"    D(1)   "error in DEC at mean equino [deg]"
    @column "SYSVEL"     D(1)   "systemic radial velocity [m/s]"
    @column "VELTYP"     A(8)   "reference for radial velocity"
    @column "VELDEF"     A(8)   "definition of radial velocity"
    @column "PMRA"       D(1)   "proper motion in RA [deg/yr]"
    @column "PMDEC"      D(1)   "proper motion in DEC [deg/yr]"
    @column "PMRA_ERR"   D(1)   "error of proper motion in RA [deg/yr]"
    @column "PMDEC_ERR"  D(1)   "error of proper motion in DEC [deg/yr]"
    @column "PARALLAX"   E(1)   "parallax [deg]"
    @column "PARA_ERR"   E(1)   "error in parallax [deg]"
    @column "SPECTYP"    A(16)  "spectral type"
])

# OI_TARGET definition (2nd revision):
define("OI_TARGET", 2, [
    @header "OI_REVN"    I      "revision number of the table definition"
    #---------------------------------------------------------
    @column "TARGET_ID"  I(1)   "index number"
    @column "TARGET"     A(16)  "target name"
    @column "RAEP0"      D(1)   "RA at mean equinox [deg]"
    @column "DECEP0"     D(1)   "DEC at mean equinox [deg]"
    @column "EQUINOX"    E(1)   "equinox [yr]"
    @column "RA_ERR"     D(1)   "error in RA at mean equinox [deg]"
    @column "DEC_ERR"    D(1)   "error in DEC at mean equino [deg]"
    @column "SYSVEL"     D(1)   "systemic radial velocity [m/s]"
    @column "VELTYP"     A(8)   "reference for radial velocity"
    @column "VELDEF"     A(8)   "definition of radial velocity"
    @column "PMRA"       D(1)   "proper motion in RA [deg/yr]"
    @column "PMDEC"      D(1)   "proper motion in DEC [deg/yr]"
    @column "PMRA_ERR"   D(1)   "error of proper motion in RA [deg/yr]"
    @column "PMDEC_ERR"  D(1)   "error of proper motion in DEC [deg/yr]"
    @column "PARALLAX"   E(1)   "parallax [deg]"
    @column "PARA_ERR"   E(1)   "error in parallax [deg]"
    @column "SPECTYP"    A(16)  "spectral type"
    @column "CATEGORY"   A(3)'  "CALibrator or SCIence target"
])

# OI_ARRAY definition (1st revision):
define("OI_ARRAY", 1, [
    @header "OI_REVN"    I      "revision number of the table definition"
    @header "ARRNAME"    A      "array name for cross-referencing"
    @header "FRAME"      A      "coordinate frame"
    @header "ARRAYX"     D      "array center X-coordinate [m]"
    @header "ARRAYY"     D      "array center Y-coordinate [m]"
    @header "ARRAYZ"     D      "array center Z-coordinate [m]"
    #------------------------------------------------------------
    @column "TEL_NAME"   A(16)  "telescope name"
    @column "STA_NAME"   A(16)  "station name"
    @column "STA_INDEX"  I(1)   "station index"
    @column "DIAMETER"   E(1)   "element diameter [m]"
    @column "STAXYZ"     D(3)   "station coordinates relative to array center [m]"
])

# OI_ARRAY definition (2nd revision):
define("OI_ARRAY", 2, [
    @header "OI_REVN"    I      "revision number of the table definition"
    @header "ARRNAME"    A      "array name for cross-referencing"
    @header "FRAME"      A      "coordinate frame"
    @header "ARRAYX"     D      "array center X-coordinate [m]"
    @header "ARRAYY"     D      "array center Y-coordinate [m]"
    @header "ARRAYZ"     D      "array center Z-coordinate [m]"
    #------------------------------------------------------------
    @column "TEL_NAME"   A(16)  "telescope name"
    @column "STA_NAME"   A(16)  "station name"
    @column "STA_INDEX"  I(1)   "station index"
    @column "DIAMETER"   E(1)   "element diameter [m]"
    @column "STAXYZ"     D(3)   "station coordinates relative to array center [m]"
    @column "FOV"        D(1)   "photometric field of view [arcsec]"
    @column "FOVTYPE"    A(6)   "model for FOV: 'FWHM' or 'RADIUS'"
])

# OI_WAVELENGTH definition (1st and 2nd revisions):
define("OI_WAVELENGTH", (1, 2), [
    @header "OI_REVN"    I    "revision number of the table definition"
    @header "INSNAME"    A    "name of detector for cross-referencing"
    #------------------------------------------------------------
    @column "EFF_WAVE"  E(1)  "effective wavelength of channel [m]"
    @column "EFF_BAND"  E(1)  "effective bandpass of channel [m]"
])

# OI_VIS definition (1st revision):
define("OI_VIS", 1, [
    @header "OI_REVN"    I     "revision number of the table definition"
    @header "DATE-OBS"   A     "UTC start date of observations"
    @header "ARRNAME"    A'    "name of corresponding array"
    @header "INSNAME"    A     "name of corresponding detector"
    #------------------------------------------------------------
    @column "TARGET_ID"  I(1)  "target number as index into OI_TARGET table"
    @column "TIME"       D(1)  "UTC time of observation [s]"
    @column "MJD"        D(1)  "modified Julian Day [day]"
    @column "INT_TIME"   D(1)  "integration time [s]"
    @column "VISAMP"     D(W)  "visibility amplitude"
    @column "VISAMPERR"  D(W)  "error in visibility amplitude"
    @column "VISPHI"     D(W)  "visibility phase [deg]"
    @column "VISPHIERR"  D(W)  "error in visibility phase [deg]"
    @column "UCOORD"     D(1)  "U coordinate of the data [m]"
    @column "VCOORD"     D(1)  "V coordinate of the data [m]"
    @column "STA_INDEX"  I(2)  "station numbers contributing to the data"
    @column "FLAG"       L(W)  "flag"
])

# FIXME: In OI-FITS Rev. 2, only MJD and DATE-OBS must be used to express time. The TIME
#        column is retained only for backwards compatibility and must contain zeros.

# OI_VIS definition (2nd revision):
define("OI_VIS", 2, [
    @header "OI_REVN"          I        "revision number of the table definition"
    @header "DATE-OBS"         A        "UTC start date of observations"
    @header "ARRNAME"          A        "name of corresponding array"
    @header "INSNAME"          A        "name of corresponding detector"
    @header "CORRNAME"         A'       "name of corresponding correlation table"
    @header "AMPTYP"           A'       "'ABSOLUTE', 'DIFFERENTIAL' or 'CORRELATED FLUX'"
    @header "PHITYP"           A'       "'ABSOLUTE', or 'DIFFERENTIAL'"
    @header "AMPORDER"         I'       "polynomial fit order for differential amplitudes"
    @header "PHIORDER"         I'       "polynomial fit order for differential phases"
    #------------------------------------------------------------------------
    @column "TARGET_ID"        I(1)     "target number as index into OI_TARGET table"
    @column "TIME"             D(1)     "UTC time of observation [s]"
    @column "MJD"              D(1)     "modified Julian Day [day]"
    @column "INT_TIME"         D(1)     "integration time [s]"
    @column "VISAMP"           D(W)     "visibility amplitude"
    @column "VISAMPERR"        D(W)     "error in visibility amplitude"
    @column "CORRINDX_VISAMP"  J(1)'    "index into correlation matrix for 1st VISAMP element"
    @column "VISPHI"           D(W)     "visibility phase [deg]"
    @column "VISPHIERR"        D(W)     "error in visibility phase [deg]"
    @column "CORRINDX_VISPHI"  J(1)'    "index into correlation matrix for 1st VISPHI element"
    @column "VISREFMAP"        L(W,W)'  "true where spectral channels were taken as reference for differential visibility computation"
    @column "RVIS"             D(W)'    "real part of complex coherent flux"
    @column "RVISERR"          D(W)'    "error on RVIS"
    @column "CORRINDX_RVIS"    J(1)'    "index into correlation matrix for 1st RVIS element"
    @column "IVIS"             D(W)'    "imaginary part of complex coherent flux"
    @column "IVISERR"          D(W)'    "error on IVIS"
    @column "CORRINDX_IVIS"    J(1)'    "index into correlation matrix for 1st IVIS element"
    @column "UCOORD"           D(1)     "U coordinate of the data [m]"
    @column "VCOORD"           D(1)     "V coordinate of the data [m]"
    @column "STA_INDEX"        I(2)     "station numbers contributing to the data"
    @column "FLAG"             L(W)     "flag"
])

# OI_VIS2 definition (1st revision):
define("OI_VIS2", 1, [
    @header "OI_REVN"    I     "revision number of the table definition"
    @header "DATE-OBS"   A     "UTC start date of observations"
    @header "ARRNAME"    A'    "name of corresponding array"
    @header "INSNAME"    A     "name of corresponding detector"
    #------------------------------------------------------------
    @column "TARGET_ID"  I(1)  "target number as index into OI_TARGET table"
    @column "TIME"       D(1)  "UTC time of observation [s]"
    @column "MJD"        D(1)  "modified Julian Day [day]"
    @column "INT_TIME"   D(1)  "integration time [s]"
    @column "VIS2DATA"   D(W)  "squared visibility"
    @column "VIS2ERR"    D(W)  "error in squared visibility"
    @column "UCOORD"     D(1)  "U coordinate of the data [m]"
    @column "VCOORD"     D(1)  "V coordinate of the data [m]"
    @column "STA_INDEX"  I(2)  "station numbers contributing to the data"
    @column "FLAG"       L(W)  "flag"
])

# OI_VIS2 definition (2nd revision):
define("OI_VIS2", 2, [
    @header "OI_REVN"            I      "revision number of the table definition"
    @header "DATE-OBS"           A      "UTC start date of observations"
    @header "ARRNAME"            A      "name of corresponding array"
    @header "INSNAME"            A      "name of corresponding detector"
    @header "CORRNAME"           A'     "name of corresponding correlation table"
    #------------------------------------------------------------
    @column "TARGET_ID"          I(1)   "target number as index into OI_TARGET table"
    @column "TIME"               D(1)   "UTC time of observation [s]"
    @column "MJD"                D(1)   "modified Julian Day [day]"
    @column "INT_TIME"           D(1)   "integration time [s]"
    @column "VIS2DATA"           D(W)   "squared visibility"
    @column "VIS2ERR"            D(W)   "error in squared visibility"
    @column "CORRINDX_VIS2DATA"  J(1)'  "index into correlation matrix for 1st VIS2DATA element"
    @column "UCOORD"             D(1)   "U coordinate of the data [m]"
    @column "VCOORD"             D(1)   "V coordinate of the data [m]"
    @column "STA_INDEX"          I(2)   "station numbers contributing to the data"
    @column "FLAG"               L(W)   "flag"
])

# OI_T3 definition (1st revision):
define("OI_T3", 1, [
    @header "OI_REVN"    I     "revision number of the table definition"
    @header "DATE-OBS"   A     "UTC start date of observations"
    @header "ARRNAME"    A'    "name of corresponding array"
    @header "INSNAME"    A     "name of corresponding detector"
    #------------------------------------------------------------
    @column "TARGET_ID"  I(1)  "target number as index into OI_TARGET table"
    @column "TIME"       D(1)  "UTC time of observation [s]"
    @column "MJD"        D(1)  "modified Julian Day [day]"
    @column "INT_TIME"   D(1)  "integration time [s]"
    @column "T3AMP"      D(W)  "triple product amplitude"
    @column "T3AMPERR"   D(W)  "error in triple product amplitude"
    @column "T3PHI"      D(W)  "triple product phase [deg]"
    @column "T3PHIERR"   D(W)  "error in triple product phase [deg]"
    @column "U1COORD"    D(1)  "U coordinate of baseline AB of the triangle [m]"
    @column "V1COORD"    D(1)  "V coordinate of baseline AB of the triangle [m]"
    @column "U2COORD"    D(1)  "U coordinate of baseline BC of the triangle [m]"
    @column "V2COORD"    D(1)  "V coordinate of baseline BC of the triangle [m]"
    @column "STA_INDEX"  I(3)  "station numbers contributing to the data"
    @column "FLAG"       L(W)  "flag"
])

# OI_T3 definition (2nd revision):
define("OI_T3", 2, [
    @header "OI_REVN"         I      "revision number of the table definition"
    @header "DATE-OBS"        A      "UTC start date of observations"
    @header "ARRNAME"         A      "name of corresponding array"
    @header "INSNAME"         A      "name of corresponding detector"
    @header "CORRNAME"        A'     "name of corresponding correlation table"
    #------------------------------------------------------------
    @column "TARGET_ID"       I(1)   "target number as index into OI_TARGET table"
    @column "TIME"            D(1)   "UTC time of observation [s]"
    @column "MJD"             D(1)   "modified Julian Day [day]"
    @column "INT_TIME"        D(1)   "integration time [s]"
    @column "T3AMP"           D(W)   "triple product amplitude"
    @column "T3AMPERR"        D(W)   "error in triple product amplitude"
    @column "CORRINDX_T3AMP"  J(1)'  "index into correlation matrix for 1st T3AMP element"
    @column "T3PHI"           D(W)   "triple product phase [deg]"
    @column "T3PHIERR"        D(W)   "error in triple product phase [deg]"
    @column "CORRINDX_T3PHI"  J(1)'  "index into correlation matrix for 1st T3PHI element"
    @column "U1COORD"         D(1)   "U coordinate of baseline AB of the triangle [m]"
    @column "V1COORD"         D(1)   "V coordinate of baseline AB of the triangle [m]"
    @column "U2COORD"         D(1)   "U coordinate of baseline BC of the triangle [m]"
    @column "V2COORD"         D(1)   "V coordinate of baseline BC of the triangle [m]"
    @column "STA_INDEX"       I(3)   "station numbers contributing to the data"
    @column "FLAG"            L(W)   "flag"
])

# OI_FLUX definition (1st revision):
define("OI_FLUX", 1, [
    @header "OI_REVN"            I      "revision number of the table definition"
    @header "DATE-OBS"           A      "UTC start date of observations"
    @header "INSNAME"            A      "name of corresponding detector"
    @header "ARRNAME"            A'     "name of corresponding array"
    @header "CORRNAME"           A'     "name of corresponding correlation table"
    @header "FOV"                D'     "area of sky over which flux is integrated [arcsec]"
    @header "FOVTYPE"            A'     "model for FOV: 'FWHM' or 'RADIUS'"
    @header "CALSTAT"            A      "'C': spectrum is calibrated, 'U': uncalibrated"
    #------------------------------------------------------------
    @column "TARGET_ID"          I(1)   "target number as index into OI_TARGET table"
    @column "MJD"                D(1)   "modified Julian Day [day]"
    @column "INT_TIME"           D(1)   "integration time [s]"
    @column "FLUXDATA"           D(W)   "flux"
    @column "FLUXERR"            D(W)   "flux error"
    @column "CORRINDX_FLUXDATA"  J(1)'  "index into correlation matrix for 1st FLUXDATA element"
    @column "STA_INDEX"          I(1)'  "station number contributing to the data"
    @column "FLAG"               L(W)   "flag"
])

# OI_CORR definition (1st revision):
define("OI_CORR", 1, [
    @header "OI_REVN"   I     "revision number of the table definition"
    @header "CORRNAME"  A     "name of correlation data set"
    @header "NDATA"     I     "number of correlated data"
    #------------------------------------------------------------
    @column "IINDX"     J(1)  "1st index of correlation matrix element"
    @column "JINDX"     J(1)  "2nd index of correlation matrix element"
    @column "CORR"      D(1)  "matrix element"
])

# OI_INSPOL definition (1st revision):
define("OI_INSPOL", 1, [
    @header "OI_REVN"    I      "revision number of the table definition"
    @header "DATE-OBS"   A      "UTC start date of observations"
    @header "NPOL"       I      "number of polarization types in this table"
    @header "ARRNAME"    A      "identifies corresponding OI_ARRAY"
    @header "ORIENT"     A      "orientation of the Jones Matrix: 'NORTH' (for on-sky orientation), or 'LABORATORY'"
    @header "MODEL"      A      "describe the way the Jones matrix is estimated"
    #---------------------------------------------------------------------------------------
    @column "TARGET_ID"  I(1)   "target number as index into OI_TARGET table"
    @column "INSNAME"    A(16)  "INSNAME of this polarization"
    @column "MJD_OBS"    D(1)   "modified Julian day, start of time lapse"
    @column "MJD_END"    D(1)   "modified Julian day, end of time lapse"
    @column "JXX"        C(W)   "complex Jones Matrix component along X axis"
    @column "JYY"        C(W)   "complex Jones Matrix component along Y axis"
    @column "JXY"        C(W)   "complex Jones Matrix component between Y and X axis"
    @column "JYX"        C(W)   "complex Jones Matrix component between Y and X axis"
    @column "STA_INDEX"  I(1)'  "station number for the above matrices"
])

end # module
