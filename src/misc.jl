#
# misc.jl --
#
# Implement reading/writing of FITS data from/to FITS files.
#
#------------------------------------------------------------------------------

using FITSIO
using FITSIO.Libcfitsio

import FITSIO: TableHDU, ASCIITableHDU
import Base: read

# FIXME: Extend `get(hdr,key,def)` for FITSHeader if not yet done.
if length(methods(get, (FITSHeader, String, Any))) < 1
    Base.get(hdr::FITSHeader, key::String, def) =
        ((i = get(hdr.map, key, -1)) > 0 ? hdr.values[i] : def)
end

"""
   OIFITS.get_hdu_type(arg)

yields the FITS Header Data Unit (HDU) type of `arg` as one of the folloing
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
        return get_hdu_type(fixname(val))
    elseif get(hdr, "SIMPLE", false) == true
        return :image_hdu
    else
        return :unknown
    end
end

# Convert low-level handle into high level HDU type.
function make_hdu(ff::FITSFile)
    fits_assert_open(ff)
    hdutype = fits_get_hdu_type(ff)
    n = fits_get_hdu_num(ff)
    hdutype == :image_hdu    ? ImageHDU(ff, n) :
    hdutype == :binary_table ? TableHDU(ff, n) :
    hdutype == :ascii_table  ? ASCIITableHDU(ff, n) :
    error("current FITS HDU is not a table")
end

# Low level version.
read_table(ff::FITSFile) = read_table(make_hdu(ff))

# FIXME: This can't work.
function read_table(hdu::Union{TableHDU,ASCIITableHDU})
    hdr = read_header(hdu)
    data = Dict{String,Any}()
    ncols = get_integer(hdr, "TFIELDS", 0)
    for k in 1:ncols
        name = uppercase(strip(get_string(hdr, "TTYPE$k", "")))
        if haskey(data, name)
            @warn "duplicate column name: \"$name\""
            continue
        end
        data[name] = read_column(hdu.fitsfile, k)
        units = strip(get_string(hdr, "TUNIT$k", ""))
        if length(units) > 0
            data[name*".units"] = units
        end
    end
    return data
end

# Read the entire table from disk. (High level version.)
read(hdu::Union{TableHDU,ASCIITableHDU}) = read_table(hdu)
