#
# misc.jl --
#
# Implement reading/writing of FITS data from/to FITS files.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015-2019: Éric Thiébaut.
#
#------------------------------------------------------------------------------

using FITSIO
using FITSIO.Libcfitsio

import FITSIO: TableHDU, ASCIITableHDU
import Base: read

const _EXTENSION = Dict("IMAGE" => :image_hdu,
                        "TABLE" => :ascii_table,
                        "BINTABLE" => :binary_table)

# Guess HDU type (no warranty to work for primary HDU nor for incomplete
# header).
function get_hdutype(hdr::FITSHeader)
    if haskey(hdr, "XTENSION")
        return get(_EXTENSION, uppercase(rstrip(hdr["XTENSION"])), :unknown)
    elseif haskey(hdr, "SIMPLE") && hdr["SIMPLE"] == true
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

function read_table(hdu::Union{TableHDU,ASCIITableHDU})
    hdr = read_header(hdu)
    data = Dict{Name,Any}()
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
