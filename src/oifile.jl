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
    # Make sure FIT file is open.
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

# Make `OIHeader` indexable.
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

# Get the type of the HDU.
oifits_get_hdutype(hdr::OIHeader) = hdr.hdutype

# Get the column number, return -1 if keyword not present.
function oifits_get_colnum(hdr::OIHeader, colname::String)
    get(hdr.columns, fixname(colname), -1)
end

# Get the type of the data-block.
function oifits_get_dbtype(hdr::OIHeader)
    if hdr.hdutype == :binary_table
        extname = fixname(oifits_get_string(hdr, "EXTNAME", ""))
        if beginswith(extname, "OI_")
            return symbol(replace(extname, r"[^A-Z0-9_]", '_'))
        end
    end
    :unknown
end

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
            columns[fixname(ttype)] = k
        end
    end
    return hdr
end

function get_part(hdr::OIHeader, key::String, idx::Integer, def)
    haskey(hdr.contents, key) ? hdr.contents[key][idx] : def
end

function get_part(hdr::OIHeader, key::String, idx::Integer)
    haskey(hdr.contents, key) || error("missing FITS keyword $key")
    hdr.contents[key][idx]
end

oifits_get_value(hdr::OIHeader, key::String) = get_part(hdr, key, 1)
oifits_get_value(hdr::OIHeader, key::String, def) = get_part(hdr, key, 1, def)
oifits_get_comment(hdr::OIHeader, key::String) = get_part(hdr, key, 2)
oifits_get_comment(hdr::OIHeader, key::String,
                   def::String) = get_part(hdr, key, 2, def)

for (fn, T, S) in ((:oifits_get_integer, Integer, Int),
                   (:oifits_get_real,    Real,    Cdouble),
                   (:oifits_get_logical, Bool,    Bool),
                   (:oifits_get_string,  String,  ASCIIString))
    @eval begin
        function $fn(hdr::OIHeader, key::String, def::$T)
            val = get_part(hdr, key, 1, def)
            isa(val, $T) || error("bad type for FITS keyword $key")
            return typeof(val) != $S ? convert($S, val) : val
        end
        function $fn(hdr::OIHeader, key::String)
            val = get_part(hdr, key, 1)
            isa(val, $T) || error("bad type for FITS keyword $key")
            return typeof(val) != $S ? convert($S, val) : val
        end
    end
end

# Returns invalid result if not a valid OI-FITS data-block.
# Unless quiet is true, print warn message.
function check_datablock(hdr::OIHeader; quiet::Bool=false)
    # Values returned in case of error.
    dbname = ""
    dbrevn = -1
    dbdefn = nothing

    # Use a while loop to break out whenever an error occurs.
    while hdr.hdutype == :binary_table
        # Get extension name.
        extname = oifits_get_value(hdr, "EXTNAME", nothing)
        if ! isa(extname, String)
            quiet || warn(extname == nothing ? "missing keyword EXTNAME"
                                             : "EXTNAME value is not a string")
            break
        end
        extname = fixname(extname)
        beginswith(extname, "OI_") || break
        dbname = extname
        if ! haskey(_DATABLOCKS, dbname)
            quiet || warn("unknown OI-FITS data-block \"$extname\"")
            break
        end

        # Get revision number.
        revn = oifits_get_value(hdr, "OI_REVN", nothing)
        if ! isa(revn, Integer)
            quiet || warn(revn == nothing ? "missing keyword OI_REVN"
                                          : "OI_REVN value is not an integer")
            break
        end
        dbrevn = revn
        if dbrevn <= 0
            quiet || warn("invalid OI_REVN value ($dbrevn)")
            break
        end
        if dbrevn > length(_FORMATS)
            quiet || warn("unsupported OI_REVN value ($dbrevn)")
            break
        end
        if ! haskey(_FORMATS[dbrevn], dbname)
            quiet || warn("unknown OI-FITS data-block \"$extname\"")
        end
        dbdefn = _FORMATS[dbrevn][dbname]
        break
    end
    return (dbname, dbrevn, dbdefn)
end

function oifits_read_datablock(ff::FITSFile; quiet::Bool=false)
    oifits_read_datablock(ff, oifits_read_header(ff), quiet=quiet)
end

function oifits_read_datablock(ff::FITSFile, hdr::OIHeader; quiet::Bool=false)
    (dbtype, revn, defn) = check_datablock(hdr, quiet=quiet)
    defn == nothing && return nothing
    nerrs = 0
    data = Dict{Symbol,Any}([:revn => revn])
    for field in defn.fields
        spec = defn.spec[field]
        name = spec.name
        if spec.keyword
            value = oifits_get_value(hdr, name, nothing)
            if value == nothing
                warn("missing keyword \"$name\" in OI-FITS $dbtype data-block")
                ++nerrs
            else
                data[field] = value
            end
        else
            colnum = oifits_get_colnum(hdr, name)
            if colnum < 1
                warn("missing column \"$name\" in OI-FITS $dbtype data-block")
                ++nerrs
            else
                data[field] = oifits_read_column(ff, colnum)
            end
        end
    end
    nerrs > 0 && error("bad OI-FITS $dbtype data-block")
    return build_datablock(dbtype, revn, data)
end

function oifits_load(filename::String; quiet::Bool=false, update::Bool=true)
    return oifits_load(fits_open_file(filename), quiet=quiet, update=update)
end

function oifits_load(ff::FITSFile; quiet::Bool=false, update::Bool=true)
    master = oifits_new_master()

    # Read all contents, skipping first HDU.
    for hdu in 2:fits_get_num_hdus(ff)
        fits_movabs_hdu(ff, hdu)
        db = oifits_read_datablock(ff, quiet=quiet)
        if db == nothing
            quiet || println("skipping HDU $hdu (no OI-FITS data)")
            continue
        end
        dbname = _EXTNAMES[typeof(db)]
        quiet || println("reading OI-FITS $dbname in HDU $hdu")
        oifits_attach!(master, db)
    end
    update && oifits_update!(master)
    return master
end

# Local Variables:
# mode: Julia
# tab-width: 8
# indent-tabs-mode: nil
# fill-column: 79
# coding: utf-8
# ispell-local-dictionary: "american"
# End:
