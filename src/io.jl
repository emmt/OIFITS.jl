#
# io.jl -
#
# Read and write OI-FITS files.
#
#-------------------------------------------------------------------------------------------

# Union of possible FITS keyword types.
const KeywordTypes = Union{Bool,Int,Cdouble,String}

# Exceptions.

MissingColumn(col::AbstractString, src::Union{FitsTableHDU,OIDataBlock}) =
    MissingColumn(col, extname(src))

show(io::IO, ::MIME"text/plain", e::MissingColumn) =
    print(io, "column \"", e.col, "\" not found in FITS extension ", e.ext)

#---------------------------------------------------------------- Writing of OI-FITS files -

function write(filename::AbstractString, ds::OIDataSet;
               overwrite::Bool=false, kwds...)
    FitsFile(filename, overwrite ? "w!" : "w") do f
        write(f, ds; kwds...)
    end
    return nothing
end

# Writing an OI-FITS file is easy: just write all datablocks. To make reading easier,
# dependencies are written first.
function write(f::FitsFile, ds::OIDataSet; quiet::Bool=false)
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

function push_keyword!(header::Vector{<:Pair{String,<:Any}}, spec::FieldDefinition, value)
    spec.optional && isa(value, AbstractFloat) && isnan(value) && return nothing
    push_keyword!(header, spec.name, value, spec.units, spec.descr)
    return nothing
end

function push_keyword!(header::Vector{<:Pair{String,<:Any}},
                       name::String, value, units::String, descr::String)
    push!(header, name => (value, (units == "" ? descr : "["*units*"] "*descr)))
    return nothing
end

push_column!(data::Vector{<:Pair{String,<:Any}}, spec::FieldDefinition, value) =
    push_column!(data, spec.name, value, spec.units)

function push_column!(data::Vector{<:Pair{String,<:Any}},
                      name::String, value, units::String)
    if units == ""
        push!(data, name => value)
    else
        push!(data, name => (value, units))
    end
    return nothing
end

function write(f::FitsFile, db::OIDataBlock)
    header = Pair{String,Any}[]
    data = Pair{String,Any}[]
    check(db) # recheck!
    for spec in get_format(db; throw_errors=true)
        # Skip undefined fields.
        isdefined(db, spec.symb) || continue
        if ndims(spec) == 0
            # Header keyword.
            push_keyword!(header, spec, getfield(db, spec.symb))
        else
            # Column keyword.
            if db isa OI_TARGET
                push_column!(data, spec, get_column(column_type(spec.type), db, spec.symb))
            else
                push_column!(data, spec, getfield(db, spec.symb))
            end
        end
    end
    write(f, header, data)
end

#---------------------------------------------------------------- Reading of OI-FITS files -

# Yields whether an exception was due to a missing FITS keyword.
missing_keyword(ex::FitsError) = ex.code == AstroFITS.CFITSIO.KEY_NO_EXIST
missing_keyword(ex::KeyError) = true
missing_keyword(ex::Exception) = false

# Yields whether an exception was due to a missing FITS column.
missing_column(ex::FitsError) = ex.code == AstroFITS.CFITSIO.COL_NOT_FOUND
missing_column(ex::Exception) = false

"""
    OIFITS.read_column(T=Array, hdu, col, def=OIFITS.unspecified) -> val

Return the content of column `col` in FITS table `hdu` and converted to array type `T`. If
the column is not part of the table, then the default value `def` is returned if specified,
otherwise a [`OIFITS.MissingColumn`](@ref)` exception is thrown. This method provides some
type-stability and add missing leading dimensions as needed.

"""
read_column(hdu::FitsTableHDU, col::String, def = unspecified) =
    read_column(Array, hdu, col, def)

function read_column(::Type{T}, hdu::FitsTableHDU, col::String,
                     def = unspecified) where {T<:Array}
    try
        return read(T, hdu, col; case=false)
    catch ex
        missing_column(ex) || rethrow()
        def === unspecified && throw(MissingColumn(col, hdu))
        return def
    end
end

# Convert column data, provide missing dimensions.
function _convert_column(::Type{T}, col::String, val::T) where {T}
    return val
end

for T in (Bool, Integer, AbstractFloat, Complex{<:AbstractFloat},
          AbstractString)
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

# To read into an existing data-set, first read the OI-FITS file (to make sure it is a
# consistent data-set), then merge.
function read!(ds::OIDataSet, args::Union{AbstractString,FitsFile}...; kwds...)
    for arg in args
        merge!(ds, read(OIDataSet, arg; kwds...))
    end
    return ds
end

OIDataSet(args::Union{AbstractString,FitsFile}...; kwds...) =
    read(OIDataSet, args...; kwds...)

function read(::Type{OIDataSet}, arg::Union{AbstractString,FitsFile},
              args::Union{AbstractString,FitsFile}...; kwds...)
    read!(read(OIDataSet, arg; kwds...),  args...; kwds...)
end

read(::Type{OIDataSet}, filename::AbstractString; kwds...) =
    FitsFile(filename, "r") do f
        read(OIDataSet, f; kwds...)
    end

# Thanks to the implemented methods, reading an OI-FITS file is not too difficult. The
# dependencies must however be read first and the OI-FITS file must be a consistent data-set
# in itself.
function read(::Type{OIDataSet}, f::FitsFile; kwds...)
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

# Type-stable methods to read a single HDU.
for T in (:OI_TARGET, :OI_ARRAY, :OI_WAVELENGTH, :OI_CORR, :OI_VIS,
          :OI_VIS2, :OI_T3, :OI_FLUX, :OI_INSPOL)
    @eval begin
        $T(hdu::FitsHDU; kwds...) = read($T, hdu; kwds...)
        function read(::Type{$T}, hdu::FitsTableHDU; kwds...)
            extname(hdu) == extname($T) || error(
                "FITS table is not an ", extname($T), " extension")
            return _read($T, hdu; kwds...)
        end
    end
end
read(T::Type{<:OIDataBlock}, hdu::FitsHDU; kwds...) =
    error("FITS HDU is not an ", extname(T), " extension")

# Type-instable methods to read a single HDU. For debugging purposes.
OIDataBlock(hdu::FitsHDU; kwds...) = read(OIDataBlock, hdu; kwds...)
function read(::Type{OIDataBlock}, hdu::FitsTableHDU; kwds...)
    extn = extname(hdu)
    extn == "" && error("FITS HDU has no keyword \"EXTNAME\"")
    return _read(datablock_type(extn), hdu; kwds...)
end
read(::Type{OIDataBlock}, hdu::FitsHDU; kwds...) =
    error("FITS HDU is not a FITS table")

datablock_type(extn::AbstractString) =
    (extn == "OI_TARGET"     ? OI_TARGET     :
     extn == "OI_ARRAY"      ? OI_ARRAY      :
     extn == "OI_WAVELENGTH" ? OI_WAVELENGTH :
     extn == "OI_CORR"       ? OI_CORR       :
     extn == "OI_VIS"        ? OI_VIS        :
     extn == "OI_VIS2"       ? OI_VIS2       :
     extn == "OI_T3"         ? OI_T3         :
     extn == "OI_FLUX"       ? OI_FLUX       :
     extn == "OI_INSPOL"     ? OI_INSPOL     :
     error("\"", extn, "\" is not the name of an OI-FITS extension"))

# skip non-table HDUs
_read!(::OIDataSet, ::FitsHDU, ::Integer; kwds...) = nothing

function _read!(ds::OIDataSet, hdu::FitsTableHDU, pass::Integer; kwds...)
    extn = extname(hdu)
    if extn != ""
        if pass == 1
            # Read dependencies.
            if extn == "OI_TARGET"
                isempty(ds.target.list) || error(
                    "only one OI_TARGET data-block is allowed in ",
                    "a compliant OI-FITS file")
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

function _read(T::Type{<:OIDataBlock}, hdu::FitsTableHDU; hack_revn = undef)
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

function _read(T::Type{<:OI_TARGET}, hdu::FitsTableHDU; hack_revn = undef)
    # Read keywords.
    revn = _read_revn(T, hdu, hack_revn)
    nrows = get(Int, hdu, "NAXIS2")

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
            (revn ≥ 2 ? category[i] : ""),
        )
    end
    return OI_TARGET(list; revn=revn)
end

# Methods to allow for hacking the revision number.
_read_revn(T::Type{<:OIDataBlock}, hdu::FitsTableHDU, hack::Integer) = hack
_read_revn(T::Type{<:OIDataBlock}, hdu::FitsTableHDU, ::typeof(undef)) =
   get(Int, hdu, "OI_REVN")
_read_revn(T::Type{<:OIDataBlock}, hdu::FitsTableHDU, hack) =
    if applicable(hack, hdu)
        hack(hdu)
    else
        hack(T, get(Int, hdu, "OI_REVN"))
    end

function _read_keyword!(db::T, ::Val{S}, hdu::FitsTableHDU,
                        spec::FieldDefinition) where {T<:OIDataBlock,S}
    try
        val = get(fieldtype(T, S), hdu, spec.name)
        setfield!(db, S, val)
    catch ex
        (spec.optional && isa(ex, KeyError)) || rethrow(ex)
    end
    nothing
end

function _read_column!(db::T, ::Val{S}, hdu::FitsTableHDU,
                       spec::FieldDefinition) where {T<:OIDataBlock,S}
    try
        val = read_column(fieldtype(T, S), hdu, spec.name)
        setfield!(db, S, val)
    catch ex
        (spec.optional && isa(ex, MissingColumn)) || rethrow()
    end
    nothing
end
