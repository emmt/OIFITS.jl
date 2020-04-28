#
# files.jl --
#
# Implement reading/writing of OI-FITS data from/to FITS files.
#
#------------------------------------------------------------------------------

# Read a column from an OI-FITS table.
function read_column(ff::FITSFile, colnum::Integer, multiplier::Integer)
    # Minimal check.
    fits_assert_open(ff)

    # Get the type and the dimensions of the data stored in the column.
    (eqcoltype, repcnt, width) = fits_get_eqcoltype(ff, colnum)
    dims = fits_read_tdim(ff, colnum)
    nrows = fits_get_num_rows(ff)

    # Allocate the array and read the column contents.
    T = eqcoltype_to_type(eqcoltype) # FIXME: improve type stability
    if T <: String
        # Column contains an array of strings.  Strip the leading dimension
        # which is the maximum length of each strings.  On return trailing
        # spaces are removed (they are insignificant according to the FITS
        # norm) and the result converted to a regular string.
        if length(dims) == 1
            dims = nrows
        else
            dims[1:end-1] = dims[2:end]
            dims[end] = nrows
        end
        data = Array{T}(undef, dims...)
        fits_read_col(ff, colnum, 1, 1, data)
        return map(str -> String(rstrip(c -> c == ' ', str)), data)
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

function get_file_handle(hdu::HDU)
    fh = hdu.fitsfile
    fits_assert_open(fh)
    fits_movabs_hdu(fh, hdu.ext)
    return fh
end

"""
    OIFITS.read_table(hdu, sel=everything)

reads the FITS table in FITS Header Data Units `hdu` and returns a dictionary
whose keys are the column names and whose values are the column data.  By
default all columns are returned but a selector `sel` may be specified to
select a subset of columns, it must be a callable object such that `sel(key)`
yields `true` if column named `key` (in uppercase letters and trailing spaces
removed) is to be returned.  Units of returned columns, if any, are available
under the name `"\$(key).units"`.

"""
function read_table(hdu::Union{TableHDU,ASCIITableHDU},
                    sel = everything)
    fh = get_file_handle(hdu)
    hdr = read_header(hdu)
    data = Dict{String,Any}()
    ncols = get_integer(hdr, "TFIELDS", 0)
    for k in 1:ncols
        key = fix_name(get_string(hdr, "TTYPE$k", ""))
        if sel(key)
            if haskey(data, key)
                warn("Duplicate column name: \"", key, "\"")
                continue
            end
            data[key] = read_column(fh, k)
            units = strip(get_string(hdr, "TUNIT$k", ""))
            if length(units) > 0
                data[key*".units"] = units
            end
        end
    end
    return data
end

everything(args...) = true

"""
    OIFITS.read_datablock([T = Float64,] hdu; quiet=false)

reads the OI-FITS data in FITS Header Data Units `hdu` and returns an instance
of one of the `OIDataBlock` sub-types or `nothing` if `hdu` does not contain
OI-FITS data.  Optional argument `T` is the floating-point type to use for the
returned object.  If keyword `quiet` is `true`, no warning messages are
printed.

If keyword `quiet` is `true`, no warning messages are printed.

"""
read_datablock(hdu::HDU; kwds...) = read_datablock(Float64, hdu; kwds...),

# OI-FITS data-blocks are stored as FITS binary tables, hence returns nothing
# for any other HDU type.
read_datablock(::Type{<:AbstractFloat}, hdu::HDU; kwds...) = nothing

function read_datablock(::Type{T},
                        hdu::TableHDU;
                        quiet::Bool=false) where {T<:AbstractFloat}
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
        extname = fix_name(val)
    end
    startswith(extname, "OI_") || return nothing
    if ! (extname in EXTNAMES)
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
    local defn::Builder.DataBlockDefinition
    let val = Builder.get_definition(extname, revn, nothing)
        if val === nothing
            quiet || warn("unknown OI-FITS extension \"", extname,
                          "\" revision ", revn)
            return nothing
        end
        defn = val
    end

    # So far so good, make a dictionary of the column numbers of the table
    # indexed by the column names.
    cols = Dict{String,Int}()
    for k in 1:get_integer(hdr, "TFIELDS", 0)
        ttype = get_string(hdr, "TTYPE$k")
        cols[fix_name(ttype)] = k
    end

    # Read the columns into a new data-block instance.
    _read_datablock(get_datablock_type(T, extname), defn, hdu, hdr, cols)
end

function _read_datablock(::Type{T},
                         defn::Builder.DataBlockDefinition,
                         hdu::TableHDU,
                         hdr::FITSHeader,
                         cols::Dict{String,Int}) where {T<:OIDataBlock}
    obj = T()
    obj.revn = defn.revn
    nerrs = 0
    fh = get_file_handle(hdu)
    for field in defn.fields
        spec = defn.spec[field]
        name = spec.name
        if spec.iskeyword
            let val = get_value(hdr, name, nothing)
                if val === nothing
                    if spec.isoptional == false
                        warn("missing keyword \"", name,
                             "\" in OI-FITS extension ", defn.extname)
                        nerrs += 1
                    end
                else
                    _setproperty!(obj, spec.symb, val)
                end
            end
        else
            col = get(cols, name, 0)
            if col < 1
                if spec.isoptional == false
                    warn("missing column \"", name,
                         "\" in OI-FITS extension ", defn.extname)
                    nerrs += 1
                end
            else
                let val = read_column(fh, col, spec.multiplier)
                    _setproperty!(obj, spec.symb, val)
                end
            end
        end
    end
    nerrs > 0 && error("bad OI-FITS extension ", defn.extname)
    Builder.check_contents(obj)
    return obj
end

"""
    OIFITS.load([T = Float64,] src) -> master

reads all OI-FITS data-blocks from FITS file `src` and return an instance of
`OIMaster`.  Argument `src` can be a file name of a FITS handle.  Optional
argument `T` is the floating-point type to use for the returned object.  If
keyword `quiet` is `true`, no warning messages are printed.

"""
load(src::Union{AbstractString,FITS}; kwds...) = load(Float64, src; kwds...)

function load(::Type{T}, src::Union{AbstractString,FITS};
              kwds...) where {T<:AbstractFloat}
    load!(OIMaster{T}(), src; kwds...)
end

"""
    OIFITS.load!(dst, src) -> dst

stores all OI-FITS data-blocks from FITS file `src` in `OIMaster` instance
`dst` and return `dst`.  Argument `src` can be a file name of a FITS handle.
If keyword `quiet` is `true`, no warning messages are printed.

"""
load!(dst::OIMaster, src::AbstractString; kwds...) =
    load!(dst, FITS(src, "r"); kwds...)

function load!(dst::OIMaster{T}, src::FITS; quiet::Bool=true) where {T}
    # Read all contents, skipping first HDU.
    for hdu in 2:length(src)
        let db = read_datablock(T, src[hdu], quiet=quiet)
            if db === nothing
                quiet || println("skipping HDU ", hdu, " (no OI-FITS data)")
            else
                quiet || println("reading OI-FITS ", db.extname, " in HDU ", hdu)
                Builder._push!(dst, db)
            end
        end
    end
    Builder._update_links!(dst)
end
