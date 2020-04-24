#
# utils.jl --
#
# Implement utility methods for dealing with FITS files and for converting
# values.
#
#------------------------------------------------------------------------------

"""
   OIFITS.get_hdu_type(arg)

yields the FITS Header Data Unit (HDU) type of `arg` as one of the following
symbols: `:binary_table` for a binary table, `:image_hdu` for an image,
`:ascii_table` for an ASCII table, or `:unknown` otherwise.  Argument `arg` can
be a FITS header, a FITS HDU or the value of the XTENSION keyword.

The returned symbol should match the result of the low level method
`FITSIO.Libcfitsio.fits_get_hdu_type`.

"""
get_hdu_type(xtension::AbstractString) =
    (xtension == "BINTABLE"   ? :binary_table :
     xtension == "IMAGE"      ? :image_hdu    :
     xtension == "TABLE"      ? :ascii_table  : :unknown)
get_hdu_type(::TableHDU)      = :binary_table
get_hdu_type(::ImageHDU)      = :image_hdu
get_hdu_type(::ASCIITableHDU) = :ascii_table
get_hdu_type(::HDU)           = :unknown
function get_hdu_type(hdr::FITSHeader)
    val = get(hdr, "XTENSION", nothing)
    if isa(val, AbstractString)
        return get_hdu_type(fix_name(val))
    elseif get(hdr, "SIMPLE", false) == true
        return :image_hdu
    else
        return :unknown
    end
end

"""
   OIFITS.get_key_index(hdr, key[, def])

yields the index of the keyword `key` in FITS header `hdr`.  If `key` is not
found in `hdr`, the default value `def` is returned if it is specified,
otherwise a `KeyError` exception is thrown.

"""
get_key_index(hdr::FITSHeader, key::String) = get(hdr.map, key)
get_key_index(hdr::FITSHeader, key::String, def) = get(hdr.map, key, def)
# FIXME: This method is a hack, it should be part of FITSIO.

"""
    OIFITS.get_value(hdr, key[, def])

yields the value of keyword `key` in FITS header `hdr`.  If `key` is not found
in `hdr`, the default value `def` is returned if it is specified, otherwise an
error is thrown.

"""
get_value(hdr::FITSHeader, key::String) =
    ((i = get_key_index(hdr, key, 0)) > 0 ? hdr[i] : missing_fits_keyword(key))
get_value(hdr::FITSHeader, key::String, def) =
    ((i = get_key_index(hdr, key, 0)) > 0 ? hdr[i] : def)

@noinline missing_fits_keyword(key) =
    error("missing FITS keyword \"", key, "\"")

"""
    OIFITS.get_comment(hdr, key[, def])

yields the comment of keyword `key` in FITS header `hdr`.  If `key` is not
found in `hdr`, the default value `def` is returned if it is specified,
otherwise an error is thrown.

"""
get_comment(hdr::FITSHeader, key::String) =
    ((i = get_key_index(hdr, key, 0)) > 0 ? FITSIO.get_comment(hdr, i) :
     missing_fits_keyword(key))
get_comment(hdr::FITSHeader, key::String, def) =
    ((i = get_key_index(hdr, key, 0)) > 0 ? FITSIO.get_comment(hdr, i) : def)

"""
    OIFITS.get_integer(hdr, key[, def])

yields the value of keyword `key` in FITS header `hdr` as an integer.  If `key`
is not found in `hdr`, the default value `def` is returned if it is specified,
otherwise an error is thrown.  If `key` is found, it is checked that it can be
converted into an `Int` which is returned.  The possible types of the result
are thus `Int` (if `key` is found) or `typeof(def)` (if `key` not found and
`def` specified).

"""
function get_integer(hdr::FITSHeader, key::String, def)
    i = get_key_index(hdr, key, 0)
    i == 0 && return def
    val = hdr[i]
    isa(val, Integer) || noninteger_fits_keyword(key)
    to_integer(val)
end

function get_integer(hdr::FITSHeader, key::String)
    i = get_key_index(hdr, key, 0)
    i == 0 && missing_fits_keyword(key)
    val = hdr[i]
    isa(val, Integer) || noninteger_fits_keyword(key)
    to_integer(val)
end

@noinline noninteger_fits_keyword(key) =
    error("non-integer value for FITS keyword \"", key, "\"")

"""
    OIFITS.get_float(hdr, key[, def])

yields the value of keyword `key` in FITS header `hdr` as a floating-point.  If
`key` is not found in `hdr`, the default value `def` is returned if it is
specified, otherwise an error is thrown.  If `key` is found, it is checked that
it can be converted into a `Float64` which is returned.  The possible types of
the result are thus `Float64` (if `key` is found) or `typeof(def)` (if `key`
not found and `def` specified).

"""
function get_float(hdr::FITSHeader, key::String, def)
    i = get_key_index(hdr, key, 0)
    i == 0 && return def
    val = hdr[i]
    isa(val, Real) || nonfloat_fits_keyword(key)
    to_float(val)
end

function get_float(hdr::FITSHeader, key::String)
    i = get_key_index(hdr, key, 0)
    i == 0 && missing_fits_keyword(key)
    val = hdr[i]
    isa(val, Real) || nonfloat_fits_keyword(key)
    to_float(val)
end

@noinline nonfloat_fits_keyword(key) =
    error("non-floating point value for FITS keyword \"", key, "\"")

"""
    OIFITS.get_logical(hdr, key[, def])

yields the value of keyword `key` in FITS header `hdr` as a boolean.  If `key`
is not found in `hdr`, the default value `def` is returned if it is specified,
otherwise an error is thrown.  If `key` is found, it is checked that it can be
converted into a `Bool` which is returned.  The possible types of the result
are thus `Bool` (if `key` is found) or `typeof(def)` (if `key` not found and
`def` specified).

""" get_logical

function get_logical(hdr::FITSHeader, key::String, def)
    i = get_key_index(hdr, key, 0)
    i == 0 && return def
    val = hdr[i]
    isa(val, Bool) || nonlogical_fits_keyword(key)
    to_logical(val)
end

function get_logical(hdr::FITSHeader, key::String)
    i = get_key_index(hdr, key, 0)
    i == 0 && missing_fits_keyword(key)
    val = hdr[i]
    isa(val, Bool) || nonlogical_fits_keyword(key)
    to_logical(val)
end

@noinline nonlogical_fits_keyword(key) =
    error("non-logical value for FITS keyword \"", key, "\"")

"""
    OIFITS.get_string(hdr, key[, def]; fix=false)

yields the value of keyword `key` in FITS header `hdr` as a string.  If `key`
is not found in `hdr`, the default value `def` is returned if it is specified,
otherwise an error is thrown.  If `key` is found, it is checked that it can be
converted into a `String` which is returned.  The possible types of the result
are thus `String` (if `key` is found) or `typeof(def)` (if `key` not found and
`def` specified).

If keyword `fix` is true, trailing spaces are removed and characters converted
to upper case letters in the returned value.  This is to follow FITS
conventions that letter case and trailing spaces are insignificant when
comparing names.

"""
function get_string(hdr::FITSHeader, key::String, def; fix::Bool = false)
    i = get_key_index(hdr, key, 0)
    i == 0 && return def
    val = hdr[i]
    isa(val, AbstractString) || nonstring_fits_keyword(key)
    fix ? fix_name(val) : to_string(val)
end

function get_string(hdr::FITSHeader, key::String; fix::Bool = false)
    i = get_key_index(hdr, key, 0)
    i == 0 && missing_fits_keyword(key)
    val = hdr[i]
    isa(val, AbstractString) || nonstring_fits_keyword(key)
    fix ? fix_name(val) : to_string(val)
end

@noinline nonstring_fits_keyword(key) =
    error("non-string value for FITS keyword \"", key, "\"")

"""
   OIFITS.fix_name(str)

converts string `str` to uppercase letters and removes trailing spaces.  This
is useful to compare names according to FITS conventions that letter case and
trailing spaces are insignificant.

"""
fix_name(name::AbstractString) = to_string(uppercase(rstrip(name)))

"""
    OIFITS.to_logical(arg)

converts argument `arg` to boolean(s) of type `Bool`.

"""
to_logical(x::Bool) = x
to_logical(x::Integer) = Bool(x)
to_logical(x::Array{Bool}) = x
to_logical(x::AbstractArray{<:Integer,N}) where {N} = convert(Array{Bool,N}, x)
to_logical(x::Tuple{Vararg{Bool}}) = x
to_logical(x::Tuple{Vararg{Integer}}) = map(to_logical, x)

"""
    OIFITS.to_integer(arg)

converts argument `arg` to integer(s) of type `Int`.

"""
to_integer(x::Int) = x
to_integer(x::Integer) = Int(x)
to_integer(x::Array{Int}) = x
to_integer(x::AbstractArray{<:Integer,N}) where {N} = convert(Array{Int,N}, x)
to_integer(x::Tuple{Vararg{Int}}) = x
to_integer(x::Tuple{Vararg{Integer}}) = map(to_integer, x)

"""
    OIFITS.to_float(arg)

converts argument `arg` to float(s) of type `Float64`.

"""
to_float(x::Float64) = x
to_float(x::Real) = Float64(x)
to_float(x::Array{Float64}) = x
to_float(x::AbstractArray{<:Real,N}) where {N} = convert(Array{Float64,N}, x)
to_float(x::Tuple{Vararg{Float64}}) = x
to_float(x::Tuple{Vararg{Real}}) = map(to_float, x)

"""
    OIFITS.to_complex(arg)

converts argument `arg` to complex(es) of type `Complex{Float64}`.

"""
to_complex(x::Complex{Float64}) = x
to_complex(x::Union{Real,Complex{<:Real}}) = convert(Complex{Float64}, x)
to_complex(x::Array{Complex{Float64}}) = x
to_complex(x::AbstractArray{<:Union{Real,Complex{<:Real}},N}) where {N} =
    convert(Array{Complex{Float64},N}, x)
to_complex(x::Tuple{Vararg{Complex{Float64}}}) = x
to_complex(x::Tuple{Vararg{Union{Real,Complex{<:Real}}}}) = map(to_complex, x)

"""
    OIFITS.to_string(arg)

converts argument `arg` to string(s) of type `String`.  Argument can be a
scalar, an array or a tuple.

""" to_string
# `string(x::SubString)` calls `String(x)` so the two are as
# fast. `string(x::Symbol)` is 2-3 tiles slower than `String(x::Symbol)`.
# Hence we call `String`, not `string` in all other cases.
to_string(x::String) = x
to_string(x::AbstractString) = String(x)
to_string(x::Symbol) = String(x)
to_string(x::Array{String}) = x
to_string(x::AbstractArray{<:Union{<:AbstractString,Symbol},N}) where {N} =
    convert(Array{String,N}, x)
to_string(x::Tuple{Vararg{String}}) = x
to_string(x::Tuple{Vararg{Union{AbstractString,Symbol}}}) = map(to_string, x)

"""
    OIFITS.is_logical(arg)

yields whether argument `arg` is a logical value, or an array or a tuple of
logical values.

"""
is_logical(::Any) = false
is_logical(::Bool) = true
is_logical(::Array{<:Bool}) = true
is_logical(::Tuple{Vararg{Bool}}) = true

"""
    OIFITS.is_string(arg)

yields whether argument `arg` is a string value, or an array or a tuple of
string values.

"""
is_string(::Any) = false
is_string(::AbstractString) = true
is_string(::Array{<:AbstractString}) = true
is_string(::Tuple{Vararg{AbstractString}}) = true

"""
    OIFITS.is_integer(arg)

yields whether argument `arg` is an integer value, or an array or a tuple of
integers.

"""
is_integer(::Any) = false
is_integer(::Integer) = true
is_integer(::Array{<:Integer}) = true
is_integer(::Tuple{Vararg{Integer}}) = true

"""
    OIFITS.is_float(arg)

yields whether argument `arg` is a real value, or an array or a tuple of real
values.  The name means that `to_float(arg)` can be applied to `arg`.

"""
is_float(::Any) = false
is_float(::Real) = true
is_float(::Array{<:Real}) = true
is_float(::Tuple{Vararg{Real}}) = true

"""
    OIFITS.is_complex(arg)

yields whether argument `arg` is a complex value, or an array or a tuple of
complex values.  The name means that `to_complex(arg)` can be applied to `arg`.

"""
is_complex(::Any) = false
is_complex(::Complex) = true
is_complex(::Array{<:Complex{<:Real}}) = true
is_complex(::Tuple{Vararg{Complex{<:Real}}}) = true

"""
    OIFITS.to_fieldname(s)

converts OI-FITS column name `s` into a `Symbol` that can be used as a field
name in OI-FITS structures.  Letters are converted to lowercase and
non-alphanumeric characters replaced by underscores.

"""
to_fieldname(sym::Symbol) = sym
function to_fieldname(name::AbstractString)
    key = lowercase(name)
    key == "oi_revn" ? :revn : Symbol(replace(key, r"[^a-z0-9_]" => '_'))
end

"""
    OIFITS.warn(args...)

prints a warning message to the standard error output.

"""
warn(args...) = message(stderr, args...; head="WARNING", color=:yellow)

"""
    OIFITS.inform(args...)

prints an informational message to the standard output.

"""
inform(args...) = message(stdout, args...; head="INFO", color=:blue)

"""
    OIFITS.message(io, args...; head=nothing, color=:normal)

prints a message to `io` made of `args...` and using specified `color`.  If
`head` is not `nothing`, it is printed first in bold and followed by a space.

"""
@noinline function message(io::IO, args...; head = nothing, color::Symbol = :normal)
    head === nothing || printstyled(io, head, " "; bold=true, color=color)
    printstyled(io, args...; bold=false, color=color)
end

#
# The functions `type_to_eqcoltype`, `eqcoltype_to_type`, `bitpix_to_type` and
# `type_to_bitpix` deal with conversion between CFITSIO type code or BITPIX
# value and actual Julia data types.  They can be used as follows (assuming `T`
# is a Julia data type, while `eqcoltype` and `bitpix` are integers):
#
#     type_to_eqcoltype(T) -----------> eqcoltype (e.g., TBYTE, TFLOAT, etc.)
#     eqcoltype_to_type(eqcoltype) ---> T
#
#     type_to_bitpix(T) --------------> bitpix (e.g., BYTE_IMG, FLOAT_IMG, etc.)
#     bitpix_to_type(bitpix) ---------> T
#
# The following table gives the correspondances between CFITSIO "types",
# the BITPIX keyword and Julia types.
#
#     -------------------------------------------------
#     BITPIX   CFISTIO         Julia     Comments
#     -------------------------------------------------
#              int             Cint
#              long            Clong
#              LONGLONG        Int64     64-bit integer
#     -------------------------------------------------
#        8     BYTE_IMG        UInt8
#       16     SHORT_IMG       Int16
#       32     LONG_IMG        Int32
#       64     LONGLONG_IMG    Int64
#      -32     FLOAT_IMG       Float32
#      -64     DOUBLE_IMG      Float64
#     -------------------------------------------------
#              TBIT
#              TBYTE           Cuchar = UInt8
#              TSBYTE          Cchar = Int8
#              TLOGICAL        Bool
#              TSHORT          Cshort
#              TUSHORT         Cushort
#              TINT            Cint
#              TUINT           Cuint
#              TLONG           Clong
#              TLONGLONG       Int64
#              TULONG          Culong
#              TFLOAT          Cfloat
#              TDOUBLE         Cdouble
#              TCOMPLEX        Complex{Cfloat}
#              TDBLCOMPLEX     Complex{Cdouble}
#     -------------------------------------------------

"""
    OIFITS.type_to_bitpix(T) -> val

yields the FITS BITPIX value `val` which corresponds to Julia type `T`.

"""
type_to_bitpix(::Type{UInt8})   =   8
type_to_bitpix(::Type{Int16})   =  16
type_to_bitpix(::Type{Int32})   =  32
type_to_bitpix(::Type{Int64})   =  64
type_to_bitpix(::Type{Float32}) = -32
type_to_bitpix(::Type{Float64}) = -64

"""
    OIFITS.bitpix_to_type(val) -> T

yields the Julia type `T` which corresponds to the BITPIX value `val`.

"""
bitpix_to_type(val::Integer) =
    (val ==   8 ? UInt8   :
     val ==  16 ? Int16   :
     val ==  32 ? Int32   :
     val ==  64 ? Int64   :
     val == -32 ? Float32 :
     val == -64 ? Float64 : bad_bitpix(val))

@noinline bad_bitpix(val) =
    error(string("bad BITPIX value ", val))

"""
    OIFITS.eqcoltype_to_type(val) -> T

yields the Julia type `T` which best corresponds to the FITSIO equivalent
column type `val` or `Nothing` if there is no correspondence.

"""
eqcoltype_to_type(val::Integer) =
    (#   ==   1 ? ...              : # TBIT
     val ==  11 ? UInt8            : # TBYTE
     val ==  12 ? Int8             : # TSBYTE
     val ==  14 ? Bool             : # TLOGICAL
     val ==  16 ? String           : # TSTRING
     val ==  20 ? Cushort          : # TUSHORT      (Uint16?)
     val ==  21 ? Cshort           : # TSHORT       (Int16?)
     val ==  30 ? Cuint            : # TUINT        (Uint32?)
     val ==  31 ? Cint             : # TINT         (Int32?)
     val ==  40 ? Culong           : # TULONG
     val ==  41 ? Clong            : # TLONG
     #   ==  42 ? ...              : # TINT32BIT
     val ==  42 ? Cfloat           : # TFLOAT      (Float32?)
     val ==  81 ? Int64            : # TLONGLONG
     val ==  82 ? Cdouble          : # TDOUBLE     (Float64?)
     val ==  83 ? Complex{Cfloat}  : # TCOMPLEX    (ComplexF32?)
     val == 163 ? Complex{Cdouble} : # TDBLCOMPLEX (ComplexF64?)
     Nothing)

"""
    OIFITS.type_to_eqcoltype(T) -> val

yields the FITSIO equivalent column type `val` which best corresponds to the
Julia type `T`.

"""
type_to_eqcoltype(::Type{UInt8})            =  11
type_to_eqcoltype(::Type{Int8})             =  12
type_to_eqcoltype(::Type{Bool})             =  14
type_to_eqcoltype(::Type{String})           =  16
type_to_eqcoltype(::Type{Cushort})          =  20
type_to_eqcoltype(::Type{Cshort})           =  21
type_to_eqcoltype(::Type{Cuint})            =  30
type_to_eqcoltype(::Type{Cint})             =  31
@static if Culong !== Cuint
    type_to_eqcoltype(::Type{Culong})       =  40
end
@static if Clong !== Cint && Clong !== Int64
    type_to_eqcoltype(::Type{Clong})        =  41
end
type_to_eqcoltype(::Type{Cfloat})           =  42
type_to_eqcoltype(::Type{Int64})            =  81
type_to_eqcoltype(::Type{Cdouble})          =  82
type_to_eqcoltype(::Type{Complex{Cfloat} }) =  83
type_to_eqcoltype(::Type{Complex{Cdouble}}) = 163
