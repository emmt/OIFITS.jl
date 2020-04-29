#
# objects.jl --
#
# Implement most part of OI-FITS objects behavior.
#
#------------------------------------------------------------------------------

"""
    OIFITS.float_type(obj)

yields the floating-point type used for OI-FITS object `obj`.

"""
float_type(::OIDataBlock{T}) where {T} = T
float_type(::OIMaster{T}) where {T} = T

"""
    OIFITS.contents(obj)

yields the vector of the data-blocks owned by OI-FITS master object `obj`.

"""
contents(obj::OIMaster) = getfield(obj, :all)

"""
    OIFITS.is_attached(db)

yields whether OI-FITS data-block `db` is attached to an OI-FITS master object.

    OIFITS.is_attached(db, owner)

yields whether OI-FITS data-block `db` is attached to the OI-FITS master object
`owner`.

"""
is_attached(db::OIDataBlock) = isdefined(db, :owner)
is_attached(db::OIDataBlock{T}, owner::OIMaster{T}) where {T} =
    (is_attached(db) && db.owner === owner)

# By default, create data in double precision.
OIMaster(args...) = OIMaster{Float64}(args...)
OITarget(; kwds...) = OITarget{Float64}(; kwds...)
OIArray(; kwds...) = OIArray{Float64}(; kwds...)
OIWavelength(; kwds...) = OIWavelength{Float64}(; kwds...)
OICorrelation(; kwds...) = OICorrelation{Float64}(; kwds...)
OIVis(; kwds...) = OIVis{Float64}(; kwds...)
OIVis2(; kwds...) = OIVis2{Float64}()
OIT3(; kwds...) = OIT3{Float64}(; kwds...)
OISpectrum(; kwds...) = OISpectrum{Float64}(; kwds...)
OIPolarization(; kwds...) = OIPolarization{Float64}(; kwds...)

OIMaster{T}(args::OIDataBlock...) where {T} = OIMaster{T}(args)
function OIMaster{T}(args::Union{AbstractVector{<:OIDataBlock},
                                 Tuple{Vararg{OIDataBlock}}}) where {T}
    push!(OIMaster{T}(), args)
end

"""
    OIFITS.convert_float_type(T, obj)

converts the floating-point type of array elements stored by OI-FITS object
`obj` to be `T`.  If floating-point type of `obj` is already `T`, `obj` itself
is returned; otherwise an unlinked object is returned.

""" convert_float_type

# Constructors for copy/conversion.
for T in (:OITarget, :OIArray, :OIWavelength, :OICorrelation,
          :OIVis, :OIVis2, :OIT3, :OISpectrum, :OIPolarization)
    @eval begin
        # Constructor from an object of same kind yields same object if exactly
        # same type, a copy if floating-point type is different.
        $T{T}(obj::$T{T}) where {T} = obj
        $T{T}(obj::$T) where {T} = copyto!($T{T}(), obj)

        # Method to convert the floating-point type.
        convert_float_type(::Type{T}, obj::$T) where {T<:AbstractFloat} =
            $T{T}(obj)

        # Extend `convert` to call constructor.
        Base.convert(::Type{$T{T}}, obj::$T) where {T} = $T{T}(obj)

        # Copy to object of same kind.  Only defined and writable fields are
        # copied.
        function Base.copyto!(dst::$T{T}, src::$T) where {T}
            for sym in fieldnames($T{T})
                if iswritable($T{T}, sym) && isdefined(src, sym)
                    _setproperty!(dst, sym, getfield(src, sym))
                end
            end
            return dst
        end
    end
end

"""
    copy(db)

yields an unlinked copy of OI-FITS data-block `db`.  Only the fields of `db`
that are defined and are not links are initialized in the returned instance.
Defined data fields will be shared between `db` and the returned copy.  Links
such as `db.array`, `db.instr`, etc.  are left undefined in the returned object
which can be safely pushed into another `OIMaster` instance.  The `copyto!`
method has the same behavior of only copying writable and defined fields.

"""
Base.copy(obj::T) where {T<:OIDataBlock} = copyto!(T(), obj)

#------------------------------------------------------------------------------
# INDEXING

@inline Base.iswritable(::T, sym::Symbol) where {T<:OIDataBlock} =
    iswritable(T, sym)
@inline Base.iswritable(::Type{T}, sym::Symbol) where {T<:OIDataBlock} =
    isreadonly(T, sym) == false

@inline Base.isreadonly(::T, sym::Symbol) where {T<:OIDataBlock} =
    isreadonly(T, sym)
@inline Base.isreadonly(::Type{T}, sym::Symbol) where {T<:OIDataBlock} =
    (sym == :extname || sym == :owner)
@inline Base.isreadonly(::Type{T}, sym::Symbol) where {T<:OIData} =
    (sym == :extname || sym == :owner || sym == :array || sym == :instr ||
     sym == :correl)
@inline Base.isreadonly(::Type{T}, sym::Symbol) where {T<:OIPolarization} =
    (sym == :extname || sym == :owner || sym == :array)

# Extend `getproperty` and `setproperty!` to implement consistent `obj.field`
# syntax.  The behavior depends on the type of data-block.  Note that
# `setfield!` does no convert given value, it is the job of `setproperty!` to
# do that if needed.
@inline Base.getproperty(obj::OIDataBlock, sym::Symbol) =
    (sym == :extname ? get_extname(obj) : _getfield(obj, sym))

@inline Base.setproperty!(obj::T, sym::Symbol, val) where {T<:OIDataBlock} =
    (isreadonly(T, sym) ? read_only_field(T, sym) :
     _setproperty!(obj, sym, val))

# Like `getfield` but yields `nothing` if undefined.
@inline function _getfield(obj::T, sym::Symbol) where {T<:OIDataBlock}
    hasfield(T, sym) || type_has_no_field(T, sym)
    isdefined(obj, sym) ? getfield(obj, sym) : nothing
end

# Like `setproperty!` but by-pass checking whether field is writable.
@inline _setproperty!(obj::T, sym::Symbol, val) where {T} =
    setfield!(obj, sym, convert(fieldtype(T, sym), val))

@noinline type_has_no_field(T::Type, sym::Symbol) =
    error("type `", T, "` has no field `", sym, "`")

read_only_field(::T, sym::Symbol) where {T} =
    read_only_field(T, sym)

@noinline read_only_field(::Type{T}, sym::Symbol) where {T} =
    error("read-only field `", sym, "` in type `", T, "`")

# Make OIMaster usable as an interator and as a vector of data-blocks.
Base.length(obj::OIMaster) = length(contents(obj))
Base.size(obj::OIMaster) = (length(obj),)
Base.IndexStyle(::Type{<:OIMaster}) = IndexLinear()
Base.getindex(obj::OIMaster, i::Integer) = getindex(contents(obj), i)
Base.iterate(obj::OIMaster) = iterate(contents(obj))
Base.iterate(obj::OIMaster, state) = iterate(contents(obj), state)

# Iterating over a data-block yields (key,val) pairs where keys are the
# symbolic names of defined fields.

Base.iterate(obj::T) where {T<:OIDataBlock} =
    ((:extname, get_extname(obj)), (0, fieldnames(T)))

function Base.iterate(obj::OIDataBlock,
                      state::Tuple{Int,NTuple{N,Symbol}}) where {N}
    k, fields = state
    while k < N
        k += 1
        sym = fields[k]
        if isdefined(obj, sym)
            return ((sym, getfield(obj, sym)), (k, fields))
        end
    end
    nothing
end

#------------------------------------------------------------------------------

# Known OI-FITS extension names.
const EXTNAMES = ("OI_TARGET", "OI_WAVELENGTH", "OI_ARRAY", "OI_VIS",
                  "OI_VIS2", "OI_T3", "OI_SPECTRUM", "OI_CORR", "OI_INSPOL")

"""
    OIFITS.get_datablock_type([T = Float64,] extname)

yields the data-block type `<:OIDataBlock{T}` associated to the OI-FITS
extension `extname`.

"""
get_datablock_type(extname::Union{AbstractString,Symbol}) =
    get_datablock_type(Float64, extname)

get_datablock_type(::Type{T}, extname::Symbol) where {T<:AbstractFloat} =
    (extname === :OI_TARGET     ? OITarget{T}       :
     extname === :OI_WAVELENGTH ? OIWavelength{T}   :
     extname === :OI_ARRAY      ? OIArray{T}        :
     extname === :OI_VIS        ? OIVis{T}          :
     extname === :OI_VIS2       ? OIVis2{T}         :
     extname === :OI_T3         ? OIT3{T}           :
     extname === :OI_SPECTRUM   ? OISpectrum{T}     :
     extname === :OI_CORR       ? OICorrelation{T}  :
     extname === :OI_INSPOL     ? OIPolarization{T} :
     bad_extname(extname))

get_datablock_type(::Type{T}, extname::AbstractString) where {T<:AbstractFloat} =
    (extname == "OI_TARGET"     ? OITarget{T}       :
     extname == "OI_WAVELENGTH" ? OIWavelength{T}   :
     extname == "OI_ARRAY"      ? OIArray{T}        :
     extname == "OI_VIS"        ? OIVis{T}          :
     extname == "OI_VIS2"       ? OIVis2{T}         :
     extname == "OI_T3"         ? OIT3{T}           :
     extname == "OI_SPECTRUM"   ? OISpectrum{T}     :
     extname == "OI_CORR"       ? OICorrelation{T}  :
     extname == "OI_INSPOL"     ? OIPolarization{T} :
     bad_extname(extname))

@noinline bad_extname(extname::Union{Symbol,AbstractString}) =
    error("bad OI-FITS extension name \"", extname, "\"")


"""
    OIFITS.get_extname(arg) -> str

yields the FITS extension name of an OI-FITS data-block corresponding to `arg`.
Argument can be an OI-FITS data-block instance, an OI-FITS data-block type, or
a FITS Header Data Unit.  If `arg` is an OI-FITS data-block instance,
`arg.extname` yields the same result.

"""
get_extname(::T) where {T<:OIDataBlock} = get_extname(T)
get_extname(::Type{<:OITarget})       = "OI_TARGET"
get_extname(::Type{<:OIWavelength})   = "OI_WAVELENGTH"
get_extname(::Type{<:OIArray})        = "OI_ARRAY"
get_extname(::Type{<:OIVis})          = "OI_VIS"
get_extname(::Type{<:OIVis2})         = "OI_VIS2"
get_extname(::Type{<:OIT3})           = "OI_T3"
get_extname(::Type{<:OISpectrum})     = "OI_SPECTRUM"
get_extname(::Type{<:OICorrelation})  = "OI_CORR"
get_extname(::Type{<:OIPolarization}) = "OI_INSPOL"
get_extname(hdr::FITSHeader) =
    (get_hdu_type(hdr) !== :binary_table ? "" :
     fix_name(get_string(hdr, "EXTNAME", "")))

function show(io::IO, db::OITarget)
    print(io, "OI_TARGET: ")
    if (tgt = db.target) !== nothing
        n = length(tgt)
        if n > 0
            print(io, "target=")
            for i in 1:n
                print(io, (i == 1 ? "[\"" : ", \""),
                      tgt[i], (i == n ? "\"]" : "\""))
            end
        end
    else
        n = 0
    end
    n == 0 && print(io, "<empty>")
end

function show(io::IO, db::OIArray)
    print(io, "OI_ARRAY: ")
    ntels =  (db.sta_index === nothing ? 0 : length(db.sta_index))
    if db.arrname !== nothing
        print(io, "arrname=\"", db.arrname, "\", ", ntels, " telescope(s)")
    else
        print(io, ntels, " telescope(s)")
    end
end

function show(io::IO, db::OIWavelength)
    print(io, "OI_WAVELENGTH: ")
    if db.insname === nothing
        print(io, "<empty>")
    else
        print(io, "insname=\"", db.insname, "\"")
        if db.eff_wave !== nothing
            wave = db.eff_wave
            n = length(db.eff_wave)
            if n < 1
                print(io, " with 0 spectral channels")
            elseif n == 1
                print(io, " with 1 spectral channel at ",
                      round(wave[1]*1e6, digits=3), " µm")
            else
                print(io, " with ", n, " spectral channels from ",
                      round(minimum(wave)*1e6, digits=3), " µm to ",
                      round(maximum(wave)*1e6, digits=3), " µm")
            end
        end
    end
end

function show(io::IO, db::OISpectrum)
    print(io, "OI_SPECTRUM")
end

function _show(io::IO, db::OIDataBlock,
               name::AbstractString, nwaves::Integer, ntimes::Integer)
    print(io, name, ": ")
    if nwaves > 0 && ntimes > 0
        print(io, nwaves*ntimes, " measurements in ", nwaves,
              " spectral channel(s) and ", ntimes, " exposure(s)")
    else
        print(io, "<empty>")
    end
end

function show(io::IO, db::OIVis)
    if db.visamp !== nothing
        dims = size(db.visamp)
        nwaves = dims[1]
        ntimes = dims[2]
    elseif db.visphi !== nothing
        dims = size(db.visphi)
        nwaves = dims[1]
        ntimes = dims[2]
    else
        nwaves = 0
        ntimes = 0
    end
    _show(io, db, "OI_VIS", nwaves, ntimes)
end

function show(io::IO, db::OIVis2)
    if db.vis2data !== nothing
        dims = size(db.vis2data)
        nwaves = dims[1]
        ntimes = dims[2]
    else
        nwaves = 0
        ntimes = 0
    end
    _show(io, db, "OI_VIS2", nwaves, ntimes)
end

function show(io::IO, db::OIT3)
    if db.t3amp !== nothing
        dims = size(db.t3amp)
        nwaves = dims[1]
        ntimes = dims[2]
    elseif db.t3phi !== nothing
        dims = size(db.visphi)
        nwaves = dims[1]
        ntimes = dims[2]
    else
        nwaves = 0
        ntimes = 0
    end
    _show(io, db, "OI_T3", nwaves, ntimes)
end

function show(io::IO, master::OIMaster)
    print(io, "OI_MASTER: (", length(master), " data-block(s))")
    for db in master
        print(io, "\n    ")
        show(io, db)
    end
end

##############################
# METHODS FOR DATA SELECTION #
##############################

"""
     OIFITS.select(obj, args)

yields an iterable over the data-blocks of an OI-FITS master object `obj` whose
extension names are specifed in `args`.

"""
function select(obj::OIMaster{T}, args::AbstractString...) where {T}
    datablocks = Array{OIDataBlock{T}}(undef, 0)
    for db in obj
        if db.extname ∈ args
            push!(datablocks, db)
        end
    end
    return datablocks
end

"""
    OIFITS.select_target(src, id)

selects a subset of data from source `src` corresponding to the target
identified by `id` (an integer or a name).  The source can be an instance of
`OIMaster` or of any `OIDataBlock` sub-types.  The result is of the same type
as the source or `nothing` if the source contains no data for the given target.

The result may share part of its contents with the source `src`.

"""
select_target(obj::Union{OIMaster,OIDataBlock}, id::AbstractString) =
    _select_target(obj, fix_name(id))
select_target(obj::Union{OIMaster,OIDataBlock}, id::Integer) =
    _select_target(obj, to_integer(id))

function _select_target(obj::Union{OIMaster,OIDataBlock}, str::String)
    tgt = obj.target
    tgt === nothing && return nothing
    idx = findfirst(x -> same_name(x, str), tgt.target)
    idx === nothing && return nothing
    select_target(obj, tgt.target_id[idx])
end

function _select_target(src::OIMaster, id::Int)
    dst = OIMaster()
    for db in src
        Builder._push!(dst, select_target(db, id))
    end
    Builder._update_links!(dst)
end

_select_target(src::OIDataBlock, id::Int) = copy(src)

function _select_target(src::T, id::Int) where {T<:OITarget}
    idx = findfirst(x -> x == id, src.target_id)
    idx == nothing && return nothing
    (name, revn, defn) = Builder.get_description(src)
    dst = T()
    for (key, val) in src
        let spec = get(defn.spec, key, nothing)
            spec === nothing && continue
            if spec.iskeyword
                setfield!(dst, key, val)
            else
                setfield!(dst, key, val[idx:idx])
            end
        end
    end
    return dst
end

function _select_target(src::T, id::Int) where {T<:Union{OIVis,OIVis2,OIT3,
                                                        OISpectrum,
                                                        OIPolarization}}
    # Pre-select target.
    target_id = src.target_id
    sel = findall(x -> x == id, target_id)
    length(sel) > 0 || return nothing

    # Build a data-block with selected data starting with an empty structure.
    dst = T()
    if length(sel) == length(target_id)
        # Just copy everythings.
        copyto!(dst, src)
    else
        (name, revn, defn) = Builder.get_description(src)
        for (key, val) in src
            let spec = get(defn.spec, key, nothing)
                spec != nothing || continue
                if spec.iskeyword
                    setfield!(dst, key, val)
                else
                    # Copy a sub-array corresponding to the selection.  (The target
                    # is specified for each last index.)
                    @assert last(size(val)) == length(target_id)
                    setfield!(dst, key, _select_last(arr, sel))
                end
            end
        end
    end
    return dst
end

function _select_last(arr::AbstractArray{T,N}, sel::Vector{Int}) where {T,N}
    siz = ntuple(i -> i == N ? length(sel) : size(arr, i), Val(N))
    sub = Array{T}(undef, siz)
    if N == 1
        for j in 1:length(sel)
            sub[j] = arr[sel[j]]
        end
    elseif N == 2
        for j in 1:length(sel)
            sub[:,j] = arr[:,sel[j]]
        end
    else
        error("unexpected rank ", N)
    end
    return sub
end

function _select_first(arr::AbstractArray{T,N}, sel::Vector{Int}) where {T,N}
    siz = ntuple(i -> i == 1 ? length(sel) : size(val, i), Val(N))
    sub = Array{T}(undef, siz)
    if N == 1
        for j in 1:length(sel)
            sub[j] = val[sel[j]]
        end
    elseif N == 2
        for j in 1:length(sel)
            sub[j,:] = val[sel[j],:]
        end
    else
        error("unexpected rank ", N)
    end
    return sub
end

"""
    OIFITS.select_wavelength(src, wmin, wmax)

selects a subset of data from source `src` corresponding to the measurements at
wavelengths `w` such that `wmin ≤ w ≤ wmax`.  The wavelength bounds are in the
same units as assumed by OI-FITS convention, that is in meters.  The source can
be an instance of `OIMaster` or of any `OIDataBlock` sub-types.

An alternative is to use a selector function:

    OIFITS.select_wavelength(src, sel)

where `sel` is callable and such that `sel(w)` returns whether wavelength `w`
is to be selected.

In any cases, the result is of the same type as the source or `nothing` if the
source contains no data in the given wavelength range.

The result may share part of its contents with the source `src`.

"""
function select_wavelength(src::Union{OIMaster,OIDataBlock},
                           wavemin::Real, wavemax::Real)
    wmin = convert(Float64, wavemin)
    wmax = convert(Float64, wavemax)
    select_wavelength(src, w -> wmin ≤ w ≤ wmax)
end

function select_wavelength(src::OIMaster, select::Function)
    dst = OIMaster()
    for db in src
        Builder._push!(dst, select_wavelength(db, select))
    end
    Builder._update_links!(dst)
end

select_wavelength(src::OIDataBlock, select::Function) = copy(src)

function select_wavelength(src::T,
                           select::Function) where {T<:Union{OIWavelength,
                                                             OIVis,OIVis2,OIT3,
                                                             OISpectrum,
                                                             OIPolarization}}
    # Pre-select wavelengths.
    isintr = isa(src, OIWavelength) # source is instrument?
    wave = (isintr ? src.eff_wave : src.instr.eff_wave)
    wave === nothing && return nothing
    sel = findall(select, wave)
    length(sel) > 0 || return nothing

    # Build a data-block with selected data starting with an empty structure.
    dst = T()
    if length(sel) == length(wave)
        # Just copy everythings.
        copyto!(dst, src)
    else
        (name, revn, defn) = Builder.get_description(src)
        for (key, val) in src
            # Get field specifications.  Do nothing if field.
            let spec = get(defn.spec, key, nothing)
                spec === nothing && continue
                if spec.iskeyword || (spec.multiplier ≥ 0 && isintr)
                    # Field is a keyword or is not a wavelength-wise array.
                    setfield!(dst, key, val)
                else
                    # Copy a sub-array corresponding to the selection.  (The
                    # wavelength corresponds to the first index.)
                    @assert first(size(val)) == length(wave)
                    setfield!(dst, key, _select_first(val, sel))
                end
            end
        end
    end
    return dst
end
