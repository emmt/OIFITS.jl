#
# oifile.jl --
#
# Implement reading/writing of OI-FITS data from/to FITS files.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015, Éric Thiébaut.
#
#------------------------------------------------------------------------------

using FITSIO

include("fix-fitsio.jl")

# Read a column from a table.
function oifits_read_column(ff::FITSFile, colnum::Integer)
    # Make sure FIT file is still open.
    fits_assert_open(ff)

    # Get the type and the dimensions of the data stored in the column.
    (typecode, repcnt, width) = fits_get_eqcoltype(ff, colnum)
    dims = fits_read_tdim(ff, colnum)
    nrows = fits_get_num_rows(ff)

    # Allocate the array and read the column contents.
    T = fits_datatype(typecode)
    if T <: String
        # Column contains an array of strings.  Strip the leading dimension
        # which is the maximum length of each strings.  On return trailing
        # spaces are removed (they are insignificant according to the FITS
        # norm).
        T = ASCIIString
        if length(dims) == 1
            dims = nrows
        else
            dims[1:end-1] = dims[2:end]
            dims[end] = nrows
        end
        data = Array(T, dims...)
        fits_read_col(ff, T, colnum, 1, 1, data)
        return map(rstrip, data)
    elseif T == Nothing
        error("unsupported column data")
    else
        # Column contains numerical data.
        if length(dims) == 1 && dims[1] == 1
            # Result will be a simple vector.
            dims = nrows
        else
            # Result will be a multi-dimensional array.
            push!(dims, nrows)
        end
        data = Array(T, dims...)
        fits_read_col(ff, T, colnum, 1, 1, data)
        return data
    end
end

# `OIHeader` is the type of the result returned by `oifits_read_header()`.
# For commentary cards, the value and the comment are the same thing.
immutable OIHeader
    hdutype::Symbol                  # HDU type, one of :image_hdu,
                                     # :ascii_table, :binary_table
    contents::Dict{ASCIIString,Any}  # dictionary of (value, comment) tuples
    columns::Dict{ASCIIString,Int}   # dictionary of column indexes
end

# Make `OIHeader` indexable although read-only.
getindex(hdr::OIHeader, key) = get(hdr.contents, key, (nothing,))
function setindex!(hdr::OIHeader, key::ASCIIString, value, comment::ASCIIString)
    setindex!(hdr.contents, key, (value, comment))
end
function setindex!(hdr::OIHeader, key::ASCIIString, value::(Any,ASCIIString))
    setindex!(hdr.contents, key, value)
end
haskey(hdr::OIHeader, key) = haskey(hdr.contents, key)
keys(hdr::OIHeader) = keys(hdr.contents)

# Make `OIHeader` a valid iterator.
start(hdr::OIHeader) = start(hdr.contents)
done(hdr::OIHeader, state) = done(hdr.contents, state)
next(hdr::OIHeader, state) = next(hdr.contents, state)

const _COMMENT = Set(["HISTORY", "COMMENT"])

const _EXTENSION = ["IMAGE" => :image_hdu,
                    "TABLE" => :ascii_table,
                    "BINTABLE" => :binary_table]

# The function `oifits_read_header()` reads the header of the current HDU
# and returns it as a dictionary whose keys are the header keywords (into
# upper case letters and leading/trailing white spaces removed) and values
# are vectors of strings for commentary cards (e.g., "HISTORY", "COMMENT")
# and 2-element tuples for other cards.  These tuples have the form
# (`value`, `comment`) where `value` is the card value converted to a
# suitable type (`Bool`, `Clong`, `Cdouble` or `ASCIIString` respectively
# for FITS logical, integer, real or string cards) and `comment` is the
# comment part.  According to FITS standard, trailing spaces are
# insignificant in strings, thus string values and comments have their
# trailing spaces removed, if any.
function oifits_read_header(ff::FITSFile)
    fits_assert_open(ff)
    (nkeys,) = fits_get_hdrspace(ff)
    contents = Dict{ASCIIString,Any}()
    for k in 1:nkeys
        (keyword, value, comment) = fits_read_keyn(ff, k)
        keyword = uppercase(strip(keyword))
        comment = rstrip(comment)
        if keyword ∈ _COMMENT
            if haskey(contents, keyword)
                push!(contents[keyword][1], comment)
            else
                # Create a new entry, the tuple consists in two references to
                # the same things.
                entry = [comment]
                contents[keyword] = (entry, entry)
            end
        elseif length(keyword) >= 1
            # Guess the type of the valueaccording to its literal
            # representation.
            if haskey(contents, keyword)
                warn("duplicate FITS keyword: $keyword")
                continue
            end
            value = strip(value)
            len = length(value)
            val = nothing
            if len >= 1
                c = value[1]
                if c == 'T' && len == 1
                    val = true
                elseif c == 'F' && len == 1
                    val = false
                elseif c == '\'' && len >= 2 && value[end] == '\''
                    # This is a character string; according to FITS
                    # standard, trailing spaces are insignificant, thus
                    # we remove them.
                    val = rstrip(value[2:end-1])
                else
                    try
                        val = parseint(Clong, value)
                    catch
                        try
                            val = parsefloat(Cdouble, value)
                        catch
                        end
                    end
                end
            end
            if val == nothing
                warn("unexpected FITS value: $keyword = $value")
            else
                contents[keyword] = (val, comment)
            end
        end
    end
    if fits_get_hdu_num(ff) == 1
        hdutype = _EXTENSION["IMAGE"]
    else
        xtension, = get(contents, "XTENSION", (nothing,))
        hdutype = get(_EXTENSION, xtension, :unknown)
    end
    columns = Dict{ASCIIString,Int}()
    hdr = OIHeader(hdutype, contents, columns)
    if hdutype == :binary_table || hdutype == :ascii_table
        ncols = oifits_get_integer(hdr, "TFIELDS") # fits_get_num_cols(ff)
        for k in 1:ncols
            ttype = oifits_get_string(hdr, "TTYPE$k")
            columns[uppercase(strip(ttype))] = k
        end
    end
    return hdr
end

for (fn, T, chk) in ((:oifits_get_integer, Integer, :is_integer),
                     (:oifits_get_real,    Real,    :is_real),
                     (:oifits_get_logical, Bool,    :is_logical),
                     (:oifits_get_string,  String,  :is_string))
    @eval begin
        function $fn(hdr::OIHeader, key::String, def::$T)
            haskey(hdr.contents, key) || return def
            val = hdr.contents[key][1]
            $chk(val) || error("bad type for FITS keyword $key")
            return val
        end
        function $fn(hdr::OIHeader, key::String)
            haskey(hdr.contents, key) || error("missing FITS keyword $key")
            val = hdr.contents[key][1]
            $chk(val) || error("bad type for FITS keyword $key")
            return val
        end
    end
end
function oifits_get_comment(hdr::OIHeader, key::String, def::String)
    haskey(hdr.contents, key) ? hdr.contents[key][2] : def
end
function oifits_get_comment(hdr::OIHeader, key::String)
    haskey(hdr.contents, key) || error("missing FITS keyword $key")
    hdr.contents[key][2]
end

function oifits_read_datablock(ff::FITSFile, hdr::OIHeader)
    # Column contains an array of strings.  Strip the leading dimension
    # which is the maximum length of each strings.  On return trailing
    # spaces are removed (they are insignificant according to the FITS
    # norm).

    ncols = fits_get_num_cols(ff)

end

# Local Variables:
# mode: Julia
# tab-width: 8
# indent-tabs-mode: nil
# fill-column: 79
# coding: utf-8
# ispell-local-dictionary: "american"
# End:
