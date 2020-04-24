#
# misc.jl --
#
# Implement reading/writing of FITS data from/to FITS files.
#
#------------------------------------------------------------------------------

import Base: read # FIXME:

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
# FIXME: type-piracy!
read(hdu::Union{TableHDU,ASCIITableHDU}) = read_table(hdu)
