#
# oifile.jl --
#
# Implement reading/writing of OI-FITS data from/to FITS files.
#
#------------------------------------------------------------------------------

using FITSIO
using FITSIO.Libcfitsio

# Read a column from an OI-FITS table.
function read_column(ff::FITSFile, colnum::Integer, multiplier::Integer)
    # Minimal check.
    fits_assert_open(ff)

    # Get the type and the dimensions of the data stored in the column.
    (typecode, repcnt, width) = fits_get_eqcoltype(ff, colnum)
    dims = fits_read_tdim(ff, colnum)
    nrows = fits_get_num_rows(ff)

    # Allocate the array and read the column contents.
    T = coltype_to_type(typecode) # FIXME: improve type stability
    if T <: AbstractString
        # Column contains an array of strings.  Strip the leading dimension
        # which is the maximum length of each strings.  On return trailing
        # spaces are removed (they are insignificant according to the FITS
        # norm).
        if length(dims) == 1
            dims = nrows
        else
            dims[1:end-1] = dims[2:end]
            dims[end] = nrows
        end
        data = Array{T}(undef, dims...)
        fits_read_col(ff, colnum, 1, 1, data)
        return map(rstrip, data)
    elseif T === Nothing
        error("unsupported column data")
    else
        # Column contains numerical data.
        if length(dims) == 1 && dims[1] == multiplier == 1
            # Result will be a simple vector.
            resize!(dims, 1)
            dims[1] = nrows
        else
            # Result will be a multi-dimensional array.
            push!(dims, nrows)
        end
        data = Array{T}(undef, dims...)
        fits_read_col(ff, colnum, 1, 1, data)
        return data
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
    fix ? fixname(val) : to_string(val)
end

function get_string(hdr::FITSHeader, key::String; fix::Bool = false)
    i = get_key_index(hdr, key, 0)
    i == 0 && missing_fits_keyword(key)
    val = hdr[i]
    isa(val, AbstractString) || nonstring_fits_keyword(key)
    fix ? fixname(val) : to_string(val)
end

@noinline nonstring_fits_keyword(key) =
    error("non-string value for FITS keyword \"", key, "\"")


function get_file_handle(hdu::HDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    return hdu.fitsfile
end

"""
    OIFITS.read_datablock(hdu; quiet=false)

reads the OI-FITS data in FITS Header Data Units `hdu` and returns a 3-tuple
`(extname, revn, dict)` or `nothing` if `hdu` does not contain OI-FITS data.
The entries of the 3-tuple are the name of the FITS extension, the OI-FITS
revision number and a dictionary of the OI-FITS keywords and columns.  This
3-tuple can be directly provided to [`OIFITS.build_datablock`](@ref).

If keyword `quiet` is `true`, no warning messages are printed.

""" read_datablock

# OI-FITS data-blocks are stored as FITS binary tables, hence returns nothing
# for any other HDU type.
read_datablock(hdu::HDU; kwds...) = nothing

function read_datablock(hdu::TableHDU; quiet::Bool=false)
    # Read the header of the binary table and check extension name.
    local extname::String
    hdr = read_header(hdu)
    let val = get_value(hdr, "EXTNAME", nothing)
        if ! isa(val, AbstractString)
            quiet || warn(val === nothing
                          ? "missing keyword \"EXTNAME\""
                          : "EXTNAME value is not a string")
            return nothing
        end
        extname = fixname(val)
    end
    startswith(extname, "OI_") || return nothing
    if ! haskey(_DATABLOCKS, extname)
        quiet || warn("unknown OI-FITS extension \"", extname, "\"")
        return nothing
    end

    # Get revision number.
    local revn::Int
    let val = get_value(hdr, "OI_REVN", nothing)
        if ! isa(val, Integer)
            quiet || warn(val === nothing
                          ? "missing keyword \"OI_REVN\""
                          : "\"OI_REVN\" value is not an integer")
            return nothing
        end
        revn = val
    end
    if revn < 1
        quiet || warn("invalid \"OI_REVN\" value: ", revn)
        return nothing
    end

    # Get format definition.
    local defn::OIDataBlockDef
    let val = get_def(extname, revn, nothing)
        if val === nothing
            quiet || warn("unknown OI-FITS extension \"", extname,
                          "\" revision ", revn)
            return nothing
        end
        defn = val
    end

    # So far so good, make a dictionary of the columns of the table.
    columns = Dict{String,Int}()
    for k in 1:get_integer(hdr, "TFIELDS", 0)
        ttype = get_string(hdr, "TTYPE$k")
        columns[fixname(ttype)] = k
    end

    # Read columns contents as a dictionary.
    ff = get_file_handle(hdu)
    nerrs = 0
    dict = Dict{Symbol,Any}(:revn => revn)
    for field in defn.fields
        spec = defn.spec[field]
        name = spec.name
        if spec.iskeyword
            let val = get_value(hdr, name, nothing)
                if val === nothing
                    if spec.isoptional == false
                        warn("missing keyword \"", name,
                             "\" in OI-FITS extension ", extname)
                        nerrs += 1
                    end
                else
                    dict[field] = val
                end
            end
        else
            colnum = get(columns, name, 0)
            if colnum < 1
                if spec.isoptional == false
                    warn("missing column \"", name,
                         "\" in OI-FITS extension ", extname)
                    nerrs += 1
                end
            else
                dict[field] = read_column(ff, colnum, spec.multiplier)
            end
        end
    end
    nerrs > 0 && error("bad OI-FITS extension ", extname)
    return (extname, revn, dict)
end

"""
    OIFITS.load(f)

reads all OI-FITS data-blocks from FITS file `f` and return an instance of
`OIMaster`.  Argument `f` can be a file name of a FITS handle.  If keyword
`quiet` is `true`, no warning messages are printed.

""" load

load(filename::AbstractString; kwds...) =
    load(FITS(filename, "r"); kwds...)

function load(f::FITS; quiet::Bool=true)
    # Read all contents, skipping first HDU.
    master = new_master()
    for hdu in 2:length(f)
        let dat = read_datablock(f[hdu], quiet=quiet)
            if dat === nothing
                quiet || println("skipping HDU ", hdu, " (no OI-FITS data)")
            else
                (extname, revn, dict) = dat
                quiet || println("reading OI-FITS ", extname, " in HDU ", hdu)
                _push!(master, build_datablock(extname, revn, dict))
            end
        end
    end
    _update_links!(master)
end
