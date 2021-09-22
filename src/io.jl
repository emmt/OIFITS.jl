#
# io.jl -
#
# Read and write OI-FITS files.
#
#------------------------------------------------------------------------------

# Union of types from which OIData can be read.
const ReadInputs = Union{AbstractString,FITS}

# Union of possible FITS keyword types.
const KeywordTypes = Union{Bool,Int,Cdouble,String}

"""
    OIFITS.get_format(db,  revn=db.revn; throwerrors=false)
    OIFITS.get_format(T,   revn;         throwerrors=false)
    OIFITS.get_format(ext, revn;         throwerrors=false)

all yield OI-FITS definitions for data-block `db`, of data-block type `T`, or
of OI-FITS extension `ext` (a string or a symbol).

"""
get_format(::Type{T}, revn::Integer; kwds...) where {T<:OIDataBlock} =
    get_format(extname(Symbol, T), revn; kwds...)
get_format(db::OIDataBlock, revn::Integer=db.revn; kwds...) =
    get_format(extname(Symbol, db), revn; kwds...)
get_format(extname::Union{Symbol,AbstractString}, revn::Integer; kwds...) =
    get_format(Symbol(extname), revn; kwds...)

MissingColumn(col::AbstractString, src::Union{TableHDU,OIDataBlock}) =
    MissingColumn(col, extname(src))

MissingKeyword(key::AbstractString, src::Union{TableHDU,OIDataBlock}) =
    MissingKeyword(key, extname(src))

show(io::IO, ::MIME"text/plain", e::MissingKeyword) =
    print(io, "keyword \"", e.key, "\" not found in FITS extension ", e.ext)

show(io::IO, ::MIME"text/plain", e::MissingColumn) =
    print(io, "column \"", e.col, "\" not found in FITS extension ", e.ext)

extname(hdu::TableHDU) = read_keyword(String, hdu, "EXTNAME", "")

#------------------------------------------------------------------------------
# WRITING OF OI-FITS FILES

function write(filename::AbstractString, data::OIData;
               overwrite::Bool=false, kwds...)
    if !overwrite && isfile(filename)
        error("file \"", filename, "\" already exists, ",
              "use keyword `overwrite=true` to overwrite")
    end
    FITS(filename, "w") do f
        write(f, data; kwds...)
    end
    return nothing
end

function write(f::FITS, data::OIData; quiet::Bool=false)
    if isdefined(data, :target)
        write(f, data.target)
    elseif !quiet
        @warn "Missing OI_TARGET extension, OI-FITS file will be incomplete"
    end
    if isdefined(data, :array)
        for db in data.array
            write(f, db)
        end
    elseif !quiet
        @warn "No OI_ARRAY extensions, OI-FITS file will be incomplete"
    end
    if isdefined(data, :instr)
        for db in data.instr
            write(f, db)
        end
    elseif !quiet
        @warn "No OI_WAVELENGTH extensions, OI-FITS file will be incomplete"
    end
    if isdefined(data, :correl)
        for db in data.correl
            write(f, db)
        end
    end
    if isdefined(data, :vis)
        for db in data.vis
            write(f, db)
        end
    end
    if isdefined(data, :vis2)
        for db in data.vis2
            write(f, db)
        end
    end
    if isdefined(data, :t3)
        for db in data.t3
            write(f, db)
        end
    end
    if isdefined(data, :flux)
        for db in data.flux
            write(f, db)
        end
    end
    if isdefined(data, :inspol)
        for db in data.inspol
            write(f, db)
        end
    end
end

function write(f::FITS, db::OIDataBlock)
    hdr_names = String[]
    hdr_values = Any[]
    hdr_comments = String[]
    col_names = String[]
    col_data = Any[]
    col_units = Dict{String,String}()
    for entry in get_format(db)
        if !isdefined(db, entry.symb)
            entry.optional || error(
                "mandatory field `", entry.symb, "` is not defined in ",
                db.extname, " data-block")
            continue
        end
        if ndims(entry) == 0
            # Header keyword.
            push!(hdr_names, entry.name)
            push!(hdr_values, getfield(db, entry.symb))
            if entry.units == ""
                push!(hdr_comments, entry.descr)
            else
                push!(hdr_comments, "["*entry.units*"] "*entry.descr)
            end
        else
            # Column keyword.
            push!(col_names, entry.name)
            push!(col_data, getfield(db, entry.symb))
            col_units[entry.name] = entry.units
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

OIData(arg::ReadInputs; kwds...) = OIData{Float64}(arg; kwds...)
OIData{T}(arg::ReadInputs; kwds...) where {T<:AbstractFloat} =
    read(OIData{T}, arg; kwds...)

read(::Type{OIData}, args...; kwds...) =
    read(OIData{Float64}, args...; kwds...)

read(::Type{OIData{T}}, filename::AbstractString) where {T<:AbstractFloat} =
    read(OIData{T}, FITS(filename))

read(::Type{OIData{T}}, f::FITS) where {T<:AbstractFloat} =
    read!(OIData{T}(), f)

function read!(dest::Union{OIData{T},Vector{OIDataBlock{T}}},
               f::FITS) where {T<:AbstractFloat}
    for i in 2:length(f)
        hdu = f[i]
        extn = read_keyword(String, hdu, "EXTNAME", nothing)
        if extn !== nothing
            if extn == "OI_TARGET"
                push!(dest, _read(OITarget{T}, hdu))
            elseif extn == "OI_ARRAY"
                push!(dest, _read(OIArray{T}, hdu))
            elseif extn == "OI_WAVELENGTH"
                push!(dest, _read(OIWavelength{T}, hdu))
            elseif extn == "OI_CORR"
                push!(dest, _read(OICorr{T}, hdu))
            elseif extn == "OI_VIS"
                push!(dest, _read(OIVis{T}, hdu))
            elseif extn == "OI_VIS2"
                push!(dest, _read(OIVis2{T}, hdu))
            elseif extn == "OI_T3"
                push!(dest, _read(OIT3{T}, hdu))
            elseif extn == "OI_FLUX"
                push!(dest, _read(OIFlux{T}, hdu))
            elseif extn == "OI_INSPOL"
                push!(dest, _read(OIInsPol{T}, hdu))
            end
        end
    end
    return dest
end

_read(T::Type{<:OIDataBlock}, hdu::TableHDU) = _read!(T(), hdu)

function _read!(db::OIDataBlock, hdu::TableHDU)
    revn = read_keyword(Int, hdu, "OI_REVN")
    nrows = read_keyword(Int, hdu, "NAXIS2")
    nwaves = -1
    for spec in get_format(db, revn)
        if spec.rank == 0
            _read_keyword!(db, Val(spec.symb), hdu, spec)
        else
            _read_column!(db, Val(spec.symb), hdu, spec)
        end
    end
    return db
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
