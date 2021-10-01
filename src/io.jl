#
# io.jl -
#
# Read and write OI-FITS files.
#
#------------------------------------------------------------------------------

# Union of possible FITS keyword types.
const KeywordTypes = Union{Bool,Int,Cdouble,String}

"""
    OIFITS.column_type(sym)

yields the exact element type in a column of a FITS table whose type is
specified by the letter `sym`.

"""
column_type(sym::Symbol) = (
    sym === :A ? String :
    sym === :I ? Int16 :
    sym === :J ? Int32 :
    sym === :K ? Int64 :
    sym === :L ? Bool :
    sym === :E ? Float32 :
    sym === :D ? Float64 :
    sym === :C ? Complex{Float32} :
    sym === :M ? Complex{Float64} :
    error("unknown FITS column type `", sym, "`"))


# Exceptions.

MissingColumn(col::AbstractString, src::Union{TableHDU,OIDataBlock}) =
    MissingColumn(col, extname(src))

MissingKeyword(key::AbstractString, src::Union{TableHDU,OIDataBlock}) =
    MissingKeyword(key, extname(src))

show(io::IO, ::MIME"text/plain", e::MissingKeyword) =
    print(io, "keyword \"", e.key, "\" not found in FITS extension ", e.ext)

show(io::IO, ::MIME"text/plain", e::MissingColumn) =
    print(io, "column \"", e.col, "\" not found in FITS extension ", e.ext)

#------------------------------------------------------------------------------
# WRITING OF OI-FITS FILES

function write(filename::AbstractString, ds::OIDataSet;
               overwrite::Bool=false, kwds...)
    if !overwrite && isfile(filename)
        error("file \"", filename, "\" already exists, ",
              "use keyword `overwrite=true` to overwrite")
    end
    FITS(filename, "w") do f
        write(f, ds; kwds...)
    end
    return nothing
end

# Writing an OI-FITS file is easy: just write all datablocks.  To make reading
# easier, dependencies are written first.
function write(f::FITS, ds::OIDataSet; quiet::Bool=false)
    write(f, ds.target)
    for db in ds.array
        write(f, db)
    end
    for db in ds.instr
        write(f, db)
    end
    for db in ds.correl
        write(f, db)
    end
    for db in ds.vis
        write(f, db)
    end
    for db in ds.vis2
        write(f, db)
    end
    for db in ds.t3
        write(f, db)
    end
    for db in ds.flux
        write(f, db)
    end
    for db in ds.inspol
        write(f, db)
    end
end

function write(f::FITS, db::OIDataBlock)
    hdr_names = String[]
    hdr_values = Any[]
    hdr_comments = String[]
    col_names = String[]
    col_data = Any[]
    col_units = Dict{String,String}()
    for spec in get_format(db; throw_errors=true)
        if !isdefined(db, spec.symb)
            spec.optional || error(
                "mandatory field `", spec.symb, "` is not defined in ",
                db.extname, " data-block")
            continue
        end
        if ndims(spec) == 0
            # Header keyword.
            push!(hdr_names, spec.name)
            push!(hdr_values, getfield(db, spec.symb))
            if spec.units == ""
                push!(hdr_comments, spec.descr)
            else
                push!(hdr_comments, "["*spec.units*"] "*spec.descr)
            end
        else
            # Column keyword.
            push!(col_names, spec.name)
            push!(col_data, getfield(db, spec.symb))
            col_units[spec.name] = spec.units
        end
    end
    write(f, col_names, col_data;
          hdutype = TableHDU,
          name = db.extname,
          header = FITSHeader(hdr_names, hdr_values, hdr_comments),
          units = col_units)
end

function write(f::FITS, db::OI_TARGET)
    hdr_names = String[]
    hdr_values = Any[]
    hdr_comments = String[]
    col_names = String[]
    col_data = Any[]
    col_units = Dict{String,String}()
    for spec in get_format(db)
        if ndims(spec) == 0
            # Header keyword.
            push!(hdr_names, spec.name)
            push!(hdr_values, getfield(db, spec.symb))
            if spec.units == empty_string
                push!(hdr_comments, spec.descr)
            else
                push!(hdr_comments, "["*spec.units*"] "*spec.descr)
            end
        else
            # Column keyword.
            push!(col_names, spec.name)
            push!(col_data, get_column(column_type(spec.type),
                                       db, spec.symb))
            col_units[spec.name] = spec.units
        end
    end
    write(f, col_names, col_data;
          hdutype = TableHDU,
          name = db.extname,
          header = FITSHeader(hdr_names, hdr_values, hdr_comments),
          units = col_units)
end

#------------------------------------------------------------------------------
# READING OF OI-FITS FILES

# Yields whether an exception was due to a missing FITS keyword.
missing_keyword(ex::CFITSIO.CFITSIOError) = ex.errcode == 202
missing_keyword(ex::Exception) = false

# Yields whether an exception was due to a missing FITS column.
missing_column(ex::CFITSIO.CFITSIOError) = ex.errcode == 219
missing_column(ex::Exception) = false

"""
    OIFITS.read_keyword(T, hdu, key [, def]) -> val

yields the value of keyword `key` in FITS header of `hdu` converted to type
`T`.  If the keyword is not part of the header, then the default value `def` is
returned if specified, otherwise a `MissingKeyword` exception is thrown.  This
method provides some type-stability.

"""
function read_keyword(T::Type{<:KeywordTypes}, hdu::HDU, key::String,
                      def = unspecified)
    try
        val, com = FITSIO.read_key(hdu, key)
        return _convert_keyword(T, key, val)
    catch ex
        if missing_keyword(ex)
            def === unspecified && throw(MissingKeyword(key, hdu))
            return def
        end
        rethrow()
    end
end

for (T, S) in ((Bool, Bool),
               (Int, Integer),
               (Cdouble, AbstractFloat),
               (String, AbstractString))
    @eval _convert_keyword(::Type{T}, key::String, val::$S) where {T<:$T} =
        convert(T, val)
end

@noinline _convert_keyword(::Type{T}, key::String, val::S) where {T,S} =
    error("FITS keyword \"", key, "\" of type ", S,
          " cannot be converted to type ", T)

"""
    OIFITS.read_column([T,] hdu, col [, def]) -> val

yields the contents of column `col` in FITS table `hdu` and converted to array
type `T`.  If the column is not part of the table, then the default value `def`
is returned if specified, otherwise a `MissingColumn` exception is thrown.
This method provides some type-stability and add missing leading dimensions as
needed.

"""
function read_column(T::Type{<:Array}, hdu::TableHDU, col::String,
                     def = unspecified)
    try
        val = read(hdu, col; case_sensitive=false)
        return _convert_column(T, col, val)
    catch ex
        if missing_column(ex)
            def === unspecified && throw(MissingColumn(key, hdu))
            return def
        end
        rethrow()
    end
end

# Convert column data, provide missing dimensions.
function _convert_column(::Type{T}, col::String, val::T) where {T}
    return val
end

for T in (Bool, Integer, AbstractFloat, Complex{<:AbstractFloat}, AbstractString)
    @eval begin
        function _convert_column(::Type{Array{T,N}}, col::String,
                                 val::AbstractArray{<:$T,N}) where {T<:$T,N}
            return convert(Array{T,N}, val)
        end
        function _convert_column(::Type{Array{T,2}}, col::String,
                                 val::AbstractArray{<:$T,1}) where {T<:$T}
            # FIXME: add trailing dimension(s) instead?
            return copyto!(Array{T,2}(undef, 1, length(val)), val)
        end
        function _convert_column(::Type{Array{T,3}}, col::String,
                                 val::AbstractArray{<:$T,1}) where {T<:$T}
            # FIXME: add trailing dimension(s) instead?
            return copyto!(Array{T,3}(undef, 1, 1, length(val)), val)
        end
    end

end
@noinline _convert_column(::Type{T}, col::String, val::S) where {T,S} =
    error("FITS table column \"", col, "\" of type ", S,
          " cannot be converted to type ", T)

# To read into an existing data-set, first read the OI-FITS file (to make sure
# it is a consistent data-set), then merge.
function read!(ds::OIDataSet, args::Union{AbstractString,FITS}...; kwds...)
    for arg in args
        merge!(ds, read(OIDataSet, arg; kwds...))
    end
    return ds
end

OIDataSet(args::Union{AbstractString,FITS}...; kwds...) =
    read(OIDataSet, args...; kwds...)

function read(::Type{OIDataSet}, arg::Union{AbstractString,FITS},
              args::Union{AbstractString,FITS}...; kwds...)
    read!(read(OIDataSet, arg; kwds...),  args...; kwds...)
end

read(::Type{OIDataSet}, filename::AbstractString; kwds...) =
    FITS(filename, "r") do f
        read(OIDataSet, f; kwds...)
    end

# Thanks to the implemented methods, reading an OI-FITS file is not too
# difficult.  The dependencies must however be read first and the OI-FITS file
# must be a consistent data-set in itself.
function read(::Type{OIDataSet}, f::FITS; kwds...)
    # Starting with an empty data-set, first read all dependencies, then all
    # data.
    ds = OIDataSet()
    for pass in 1:2
        for i in 2:length(f)
            _read!(ds, f[i], pass; kwds...)
        end
    end
    return ds
end

# skip non-table HDUs
_read!(::OIDataSet, ::HDU, ::Integer; kwds...) = nothing

function _read!(ds::OIDataSet, hdu::TableHDU, pass::Integer; kwds...)
    extn = read_keyword(String, hdu, "EXTNAME", nothing)
    if extn !== nothing
        if pass == 1
            # Read dependencies.
            if extn == "OI_TARGET"
                isempty(ds.target.list) || error(
                    "only one OI_TARGET data-block is allowed in an OI-FITS file")
                push!(ds, _read(OI_TARGET, hdu; kwds...))
            elseif extn == "OI_ARRAY"
                push!(ds, _read(OI_ARRAY, hdu; kwds...))
            elseif extn == "OI_WAVELENGTH"
                push!(ds, _read(OI_WAVELENGTH, hdu; kwds...))
            elseif extn == "OI_CORR"
                push!(ds, _read(OI_CORR, hdu; kwds...))
            end
        elseif pass == 2
            # Read data.
            if extn == "OI_VIS"
                push!(ds, _read(OI_VIS, hdu; kwds...))
            elseif extn == "OI_VIS2"
                push!(ds, _read(OI_VIS2, hdu; kwds...))
            elseif extn == "OI_T3"
                push!(ds, _read(OI_T3, hdu; kwds...))
            elseif extn == "OI_FLUX"
                push!(ds, _read(OI_FLUX, hdu; kwds...))
            elseif extn == "OI_INSPOL"
                push!(ds, _read(OI_INSPOL, hdu; kwds...))
            end
        end
    end
end

function _read(T::Type{<:OIDataBlock}, hdu::TableHDU; hack_revn = undef)
    db = T(undef)
    db.revn = _read_revn(T, hdu, hack_revn)
    for spec in get_format(db)
        if spec.rank == 0
            if spec.symb !== :revn
                _read_keyword!(db, Val(spec.symb), hdu, spec)
            end
        else
            _read_column!(db, Val(spec.symb), hdu, spec)
        end
    end
    return db
end

function _read(T::Type{<:OI_TARGET}, hdu::TableHDU; hack_revn = undef)
    # Read keywords.
    revn = _read_revn(T, hdu, hack_revn)
    nrows = read_keyword(Int, hdu, "NAXIS2")

    # Read columns.
    target_id = read_column(Vector{Int16  }, hdu, "TARGET_ID")
    target    = read_column(Vector{String }, hdu, "TARGET")
    raep0     = read_column(Vector{Float64}, hdu, "RAEP0")
    decep0    = read_column(Vector{Float64}, hdu, "DECEP0")
    equinox   = read_column(Vector{Float32}, hdu, "EQUINOX")
    ra_err    = read_column(Vector{Float64}, hdu, "RA_ERR")
    dec_err   = read_column(Vector{Float64}, hdu, "DEC_ERR")
    sysvel    = read_column(Vector{Float64}, hdu, "SYSVEL")
    veltyp    = read_column(Vector{String }, hdu, "VELTYP")
    veldef    = read_column(Vector{String }, hdu, "VELDEF")
    pmra      = read_column(Vector{Float64}, hdu, "PMRA")
    pmdec     = read_column(Vector{Float64}, hdu, "PMDEC")
    pmra_err  = read_column(Vector{Float64}, hdu, "PMRA_ERR")
    pmdec_err = read_column(Vector{Float64}, hdu, "PMDEC_ERR")
    parallax  = read_column(Vector{Float32}, hdu, "PARALLAX")
    para_err  = read_column(Vector{Float32}, hdu, "PARA_ERR")
    spectyp   = read_column(Vector{String }, hdu, "SPECTYP")
    if revn ≥ 2
        category = read_column(Vector{String}, hdu, "CATEGORY")
    end
    list = Vector{OITargetEntry}(undef, nrows)
    for i in 1:nrows
        list[i] = OITargetEntry(
            target_id[i],
            target[i],
            raep0[i],
            decep0[i],
            equinox[i],
            ra_err[i],
            dec_err[i],
            sysvel[i],
            veltyp[i],
            veldef[i],
            pmra[i],
            pmdec[i],
            pmra_err[i],
            pmdec_err[i],
            parallax[i],
            para_err[i],
            spectyp[i],
            (revn ≥ 2 ? category[i] : empty_string),
        )
    end
    return OI_TARGET(list; revn=revn)
end

# Methods to allow for hacking the revision number.
_read_revn(T::Type{<:OIDataBlock}, hdu::TableHDU, hack::Integer) = hack
_read_revn(T::Type{<:OIDataBlock}, hdu::TableHDU, ::typeof(undef)) =
   read_keyword(Int, hdu, "OI_REVN")
_read_revn(T::Type{<:OIDataBlock}, hdu::TableHDU, hack) =
    if applicable(hack, hdu)
        hack(hdu)
    else
        hack(T, read_keyword(Int, hdu, "OI_REVN"))
    end

function _read_keyword!(db::T, ::Val{S}, hdu::TableHDU,
                        spec::FieldDefinition) where {T<:OIDataBlock,S}
    try
        val = read_keyword(fieldtype(T, S), hdu, spec.name)
        setfield!(db, S, val)
    catch ex
        (spec.optional && isa(ex, MissingKeyword)) || rethrow()
    end
    nothing
end

function _read_column!(db::T, ::Val{S}, hdu::TableHDU,
                       spec::FieldDefinition) where {T<:OIDataBlock,S}
    try
        val = read_column(fieldtype(T, S), hdu, spec.name)
        setfield!(db, S, val)
    catch ex
        (spec.optional && isa(ex, MissingColumn)) || rethrow()
    end
    nothing
end
