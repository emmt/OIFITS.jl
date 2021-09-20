#
# methods.jl -
#
# Methods for OI-FITS package.
#
#------------------------------------------------------------------------------

"""
    same_name(A, B)

yields whether strings `A` and `B` are the same according to FITS conventions
that case of letters and trailing spaces are irrelevant.

"""
function same_name(A::AbstractString, B::AbstractString)
    i_last = lastindex(A)
    j_last = lastindex(B)
    i = firstindex(A)
    j = firstindex(B)
    @inbounds begin
        while ((i ≤ i_last)&(j ≤ j_last)) && (uppercase(A[i]) == uppercase(B[j]))
            i = nextind(A, i)
            j = nextind(B, j)
        end
        while i ≤ i_last
            if !isspace(A[i])
                return false
            end
            i = nextind(A, i)
        end
        while j ≤ j_last
            if !isspace(B[j])
                return false
            end
            j = nextind(B, j)
        end
    end
    return true
end

"""
    fix_name(str)

yields string `str` converted to uppercase letters and with trailing spaces
removed.  This is to follows FITS conventions that names comparisons should be
done ignoring of the case of letters and trailing spaces.

"""
fix_name(str::AbstractString) = uppercase(rstrip(str))

# Extend A[name] syntax to get an OI-ARRAY, OI-WAVELENGTH, or OI-CORR by its
# name (accoding to FITS conventions).
for (T, field) in ((OIArray,      :arrname),
                   (OIWavelength, :insname),
                   (OICorr,       :corrname))
    @eval function Base.getindex(A::AbstractArray{<:$T}, name::AbstractString)
        @inbounds for i in eachindex(A)
            if same_name(A[i].$field, name)
                return A[i]
            end
        end
        return nothing
    end
end

#------------------------------------------------------------------------------
# ELEMENT TYPE CONVERSION

eltype(x::Union{OIDataBlock,OIData}) = eltype(typeof(x))
eltype(::Type{<:OIDataBlock{T}}) where {T} = T
eltype(::Type{<:OIData{T}}) where {T} = T

"""
    change_eltype(S, T)

yields type similar to `S` but with element type `T`.

""" change_eltype

for D in (:OIDataBlock, :OITarget, :OIArray, :OIWavelength, :OICorr,
          :OIVis, :OIVis2, :OIT3, :OIFlux, :OIInsPol)
    @eval begin
        change_eltype(::Type{<:$D}, T::Type{<:AbstractFloat}) = $D{T}
        convert(::Type{$D}, A::$D) = A
        convert(::Type{$D{T}}, A::$D{T}) where {T<:AbstractFloat} = A
    end
    if D === :OIDataBlock
        @eval begin
            convert(::Type{$D{T}}, A::$D) where {T<:AbstractFloat} =
                convert(change_eltype(typeof(A), T), A)
            $D(A::$D) = A
            $D{T}(A::$D{T}) where {T<:AbstractFloat} = A
            $D{T}(A::$D) where {T<:AbstractFloat} =
                convert($D{T}, A)
        end
    else
        @eval begin
            function convert(::Type{$D{T}}, A::$D) where {T<:AbstractFloat}
                B = $D{T}()
                for sym in fieldnames($D{T})
                    if isdefined(A, sym)
                        _set_field!(B, sym, getfield(A, sym))
                    end
                end
                return B
            end
        end
    end
end

convert(::Type{T}, A::T) where {T<:OIData} = A
convert(::Type{OIData{T}}, A::OIData{T}) where {T<:AbstractFloat} = A
function convert(::Type{OIData{T}}, A::OIData) where {T<:AbstractFloat}
    # Create an empty data-set and first push independent data-blocks.
    B = OIData{T}()
    if isdefined(A, :target)
        push!(B, A.target)
    end
    for db in A.array
        push!(B, db)
    end
    for db in A.instr
        push!(B, db)
    end
    for db in A.correl
        push!(B, db)
    end
    # For other data-blocks, calling push! automatically set dependencies.
    for db in A.vis
        push!(B, db)
    end
    for db in A.vis2
        push!(B, db)
    end
    for db in A.t3
        push!(B, db)
    end
    for db in A.flux
        push!(B, db)
    end
    for db in A.inspol
        push!(B, db)
    end
    return B
end

#------------------------------------------------------------------------------
# PROPERTIES

propertynames(db::OIDataBlock) = (fieldnames(typeof(db))..., :extname)

getproperty(db::OIDataBlock, sym::Symbol) = getproperty(db, Val(sym))
getproperty(db::OIDataBlock, ::Val{S}) where {S} = getfield(db, S)
getproperty(db::OIDataBlock, ::Val{:extname}) = extname(typeof(db))

# Indirections to instrument (OI_WAVELENGTH) for OI_VIS, OI_VIS2, OI_T3, and OI_FLUX.
for type in (:OIVis, :OIVis2, :OIT3, :OIFlux)
    @eval begin
        @eval propertynames(db::$type) =
            (fieldnames(typeof(db))..., :extname, :eff_wave, :eff_band)
        @eval getproperty(db::$type, ::Val{:eff_wave}) =
            getfield(db, :instr).eff_wave
        @eval getproperty(db::$type, ::Val{:eff_band}) =
            getfield(db, :instr).eff_band
    end
end

setproperty!(x::OIDataBlock, f::Symbol, v) = _set_field!(x, f, v)

@inline _get_field(A, sym::Symbol, def=undef) =
    isdefined(A, sym) ? getfield(A, sym) : def

@inline _convert_field(obj::Any, sym::Symbol, val) =
    _convert_field(typeof(obj), sym, val)

@inline _convert_field(T::Type, sym::Symbol, val) =
    try
        return convert(fieldtype(T, sym), val)
    catch ex
        _rethrow_convert_field_error(ex, T, sym, val)
    end

@noinline function _rethrow_convert_field_error(ex::Exception, T::Type,
                                                sym::Symbol, val)
    if isa(ex, MethodError) && ex.f === convert
        rethrow(ArgumentError(string(
            "Cannot `convert` an object of type `", typeof(val),
            "` to an object of type `", fieldtype(T, sym), "` for field `",
            sym, "` of structure of type `", T, "`")))
    else
        rethrow(ex)
    end
end

# `_set_field!(obj,sym,val)` is similar to `setfield!(obj,sym,val)` except that
# conversion of value `val` is automatically done and that error message is
# more informative if conversion is not possible.
@inline _set_field!(obj, sym::Symbol, val) =
    setfield!(obj, sym, _convert_field(obj, sym, val))

@inline function _define_field!(dst, src, sym::Symbol, val=undef)
    if val !== undef
         _set_field!(dst, sym, val)
    elseif isdefined(src, sym)
         _set_field!(dst, sym, getfield(src, sym))
    end
    nothing
end

#------------------------------------------------------------------------------
# BUILDING OF DATA-SETS

function OIData(args::OIDataBlock...)
    T = promote_type(map(eltype, args)...)
    return OIData{T}(args...)
end

function OIData{T}(args::OIDataBlock...) where {T<:AbstractFloat}
    # Create empty instance, then push the dependecy-less data-blocks and,
    # finally, the data-blocks with dependencies.
    data = OIData{T}()
    for db in args
        if !isa(db, DataBlocksWithDependencies)
            push!(data, db)
        end
    end
    for db in args
        if isa(db, DataBlocksWithDependencies)
            push!(data, db)
        end
    end
    return data
end

# Copy the outer structure of an OI-FITS data-block.
function copy(A::T) where {T<:OIDataBlock}
    B = T()
    for sym in fieldnames(T)
        if isdefined(A, sym)
            setfield!(B, sym, getfield(A, sym))
        end
    end
    return B
end

function push!(data::OIData, args::OIDataBlock...)
    for db in args
        push!(data, db)
    end
    return data
end

# This version is to convert element type.
push!(data::OIData{T}, db::OIDataBlock) where {T<:AbstractFloat} =
    push!(data, OIDataBlock{T}(db))

function push!(data::OIData{T}, db::OITarget{T}) where {T<:AbstractFloat}
    isdefined(data, :target) && error("OI_TARGET already defined")
    setfield!(data, :target, db)
    return data
end

# Extend push! for OI_ARRAY, OI_WAVELENGTH, and OI_CORR.
for (type, name, field) in ((:OIArray,      QuoteNode(:arrname),  QuoteNode(:array)),
                            (:OIWavelength, QuoteNode(:insname),  QuoteNode(:instr)),
                            (:OICorr,       QuoteNode(:corrname), QuoteNode(:correl)))
    @eval begin
        function push!(data::OIData{T}, db::$type{T}) where {T<:AbstractFloat}
            isdefined(db, $name) || throw_undefined_field($name)
            other = getfield(data, $field)[getfield(db, $name)]
            if other === nothing
                push!(getfield(data, $field), db)
            elseif other != db # FIXME:
                error(string("attempt to replace already existing ",
                             extname($type), " named \"", name, "\""))
            end
            return data
        end
    end
end

for (type, field) in ((:OIVis,    :vis),
                      (:OIVis2,   :vis2),
                      (:OIT3,     :t3),
                      (:OIFlux,   :flux),
                      (:OIInsPol, :inspol))
    @eval begin
        function push!(data::OIData{T}, db::$type{T}) where {T<:AbstractFloat}
            # Datablock must be copied on write to avoid side-effects.
            copy_on_write = true

            # Set/fix OI_WAVELENGTH dependency.
            isdefined(db, :insname) || throw_undefined_field(:insname)
            db, copy_on_write = _set_dependency!(
                data.instr, db, copy_on_write, Val(:insname), Val(:instr))

            # Set/fix OI_ARRAY dependency.
            if isdefined(db, :arrname)
                db, copy_on_write = _set_dependency!(
                    data.array, db, copy_on_write, Val(:arrname), Val(:array))
            end

            # Set/fix OI_CORREL dependency.
            if $type !== OIInsPol && isdefined(db, :corrname)
                db, copy_on_write = _set_dependency!(
                    data.correl, db, copy_on_write, Val(:corrname), Val(:correl))
            end

            push!(data.$field, db)
            return data
        end
    end
end

@noinline throw_undefined_field(name::Union{AbstractString,Symbol}) =
    throw(ErrorException("undefined field `$name`"))

function _set_dependency!(deps::AbstractVector, # list of dependencies known by data-set
                          db::T,                # datablock to be pushed on data-set
                          copy_on_write::Bool,  # copy datablock on write?
                          ::Val{name},          # field name of dependency name
                          ::Val{value}          # field name of dependency value
                          ) where {T<:OIDataBlock,name,value}

    other = deps[getfield(db, name)]
    if other === nothing
        # Dependency not yet in list.
        if isdefined(db, value)
            push!(deps, getfield(db, value))
        else
            error("no ", extname(eltype(deps)), " named \"",
                  fix_name(getfield(db, name)), "\" found")
        end
    else
        # Another dependency with same name already in data-set.  Check that the
        # two versions are identical.
        setfield = true
        if isdefined(db, value)
            if getfield(db, value) === other
                sefield = false
            elseif getfield(db, value) == other
                error("multiple ", extname(eltype(deps)), " named \"",
                      fix_name(getfield(db, name)), "\" found")
            end
        end
        if setfield
            if copy_on_write
                db = copy(db)
                copy_on_write = false
            end
            setfield!(db, value, other)
        end
    end
    return db, copy_on_write
end

#------------------------------------------------------------------------------
# READING OF OI-FITS FILES

const KeywordTypes = Union{Bool,Int,Cdouble,String}

"""
    read_keyword(T, hdu, key, def=nothing) -> val

yields the value of keyword `key` in FITS header of `hdu` converted to type
`T`.  If the keyword is not part of the header, `def` is returned
(unconverted).  This method provides some type-stability.

"""
function read_keyword(::Type{T}, hdu::HDU, key::String,
                      def=nothing) where {T<:KeywordTypes}
    try
        val, com = FITSIO.read_key(hdu, key)
        return _convert_keyword(T, key, val)
    catch ex
        if isa(ex, CFITSIO.CFITSIOError) && ex.errcode == 202
            # Keyword does not exist in header.
            return def
        end
        rethrow(ex)
    end
end

for (T, S) in ((Bool, Bool),
               (Int, Integer),
               (Cdouble, AbstractFloat),
               (String, AbstractString))
    @eval _convert_keyword(::Type{T}, key::String, val::$S) where {T<:$T} = convert(T, val)
end

@noinline _convert_keyword(::Type{T}, key::String, val::S) where {T,S} =
    error("FITS keyword \"", key, "\" of type ", S,
          " cannot be converted to type ", T)

"""
    read_column(T, hdu, col, def=nothing) -> val

yields the contents of column `col` in FITS table `hdu` converted to array type
`T`.  If the column is not part of the table, `def` is returned (unconverted).
This method provides some type-stability and add missing leading dimensions as
needed.

"""
function read_column(::Type{T}, hdu::TableHDU, col::String,
                     def=nothing) where {S,T<:Union{Vector{S},Matrix{S}}}
    try
        val = read(hdu, col; case_sensitive=false)
        return _convert_column(T, col, val)
    catch ex
        if isa(ex, CFITSIO.CFITSIOError) && ex.errcode == 219
            # Column does not exist in table.
            return def
        end
        rethrow(ex)
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

"""
    extname(db)

yields the extension name of OI-FITS datablock instance or type `db`.  This
method is not exported, the same result for a data-block instance is obtained
by `db.extname`.

"""
exname(db::OIDataBlock) = extname(typeof(db))
extname(::Type{<:OITarget}) = "OI_TARGET"
extname(::Type{<:OIArray}) = "OI_ARRAY"
extname(::Type{<:OIWavelength}) = "OI_WAVELENGTH"
extname(::Type{<:OICorr}) = "OI_CORR"
extname(::Type{<:OIVis}) = "OI_VIS"
extname(::Type{<:OIVis2}) = "OI_VIS2"
extname(::Type{<:OIT3}) = "OI_T3"
extname(::Type{<:OIFlux}) = "OI_FLUX"
extname(::Type{<:OIInsPol}) = "OI_INSPOL"

"""
    symbolic_extname(db)

yields the extension name of OI-FITS datablock instance or type `db` a
symbol.

"""
symbolic_extname(db::OIDataBlock) = symbolic_extname(typeof(db))
symbolic_extname(::Type{<:OITarget}) = :OI_TARGET
symbolic_extname(::Type{<:OIArray}) = :OI_ARRAY
symbolic_extname(::Type{<:OIWavelength}) = :OI_WAVELENGTH
symbolic_extname(::Type{<:OICorr}) = :OI_CORR
symbolic_extname(::Type{<:OIVis}) = :OI_VIS
symbolic_extname(::Type{<:OIVis2}) = :OI_VIS2
symbolic_extname(::Type{<:OIT3}) = :OI_T3
symbolic_extname(::Type{<:OIFlux}) = :OI_FLUX
symbolic_extname(::Type{<:OIInsPol}) = :OI_INSPOL

get_format(::Union{T,Type{T}}, rev::Integer; kwds...) where {T<:OIDataBlock} =
    get_format(symbolic_extname(T), rev; kwds...)
get_format(extname::AbstractString, rev::Integer; kwds...) =
    get_format(Symbol(extname), rev; kwds...)

const ReadInputs = Union{AbstractString,FITS}

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
        extn = read_keyword(String, hdu, "EXTNAME")
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
    if revn === nothing
        error("mandatory FITS \"OI_REVN\" keyword is missing for ", T)
    end
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
    val = read_keyword(fieldtype(T, S), hdu, spec.name)
    if val === nothing
        if !spec.optional
            error("missing mandatory keyword \"", spec.name,
                  "in FITS HDU #", hdu.ext, " \"", extname(db), "\"")
        end
    else
        setfield!(db, S, val)
    end
    nothing
end

function _read_column!(db::T, ::Val{S}, hdu::TableHDU,
                       spec::FieldDefinition) where {T<:OIDataBlock,S}
    val = read_column(fieldtype(T, S), hdu, spec.name)
    if val === nothing
        if !spec.optional
            error("missing mandatory column \"", spec.name,
                  "in FITS HDU #", hdu.ext, " \"", extname(db), "\"")
        end
    else
        setfield!(db, S, val)
    end
    nothing
end
