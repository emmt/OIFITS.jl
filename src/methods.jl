#
# methods.jl -
#
# Methods for OI-FITS package.
#
#------------------------------------------------------------------------------

"""
    OIFITS.is_same(a, b)

yields whether `a` and `b` are sufficiently identical for OI-FITS data. If `a`
and `b` are strings, trailing spaces and case of letters are irrelevant. If `a`
and `b` are arrays, they are compared element-wise. This operator is used in
OI-FITS to merge contents. For an `OITargetEntry`, the target identifier is
ignored.

"""
is_same(a, b) = (a === a)
is_same(a::Integer, b::Integer) = (a == b)
is_same(a::AbstractFloat, b::AbstractFloat) = (a == b)
is_same(a::Complex{<:AbstractFloat}, b::Complex{<:AbstractFloat}) = (a == b)

function is_same(A::AbstractArray{<:Any,N},
                 B::AbstractArray{<:Any,N}) where {N}
    axes(A) == axes(B) || return false
    @inbounds for i in eachindex(A, B)
        is_same(A[i], B[i]) || return false
    end
    return true
end

function is_same(A::AbstractString, B::AbstractString)
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

function is_same(A::OITargetEntry, B::OITargetEntry)
    # Private helper.
    @inline _is_same(A::OITargetEntry, B::OITargetEntry, sym::Symbol) =
        is_same(getfield(A, sym), getfield(B, sym))

    # Purposely ignore target_id.
    return (_is_same(A, B, :target   ) &&
            _is_same(A, B, :raep0    ) &&
            _is_same(A, B, :decep0   ) &&
            _is_same(A, B, :equinox  ) &&
            _is_same(A, B, :ra_err   ) &&
            _is_same(A, B, :dec_err  ) &&
            _is_same(A, B, :sysvel   ) &&
            _is_same(A, B, :veltyp   ) &&
            _is_same(A, B, :veldef   ) &&
            _is_same(A, B, :pmra     ) &&
            _is_same(A, B, :pmdec    ) &&
            _is_same(A, B, :pmra_err ) &&
            _is_same(A, B, :pmdec_err) &&
            _is_same(A, B, :parallax ) &&
            _is_same(A, B, :para_err ) &&
            _is_same(A, B, :spectyp  ) &&
            _is_same(A, B, :category ))
end

"""
    OIFITS.is_same(A, B, sym)

yields whether field `sym` of data-blocks `A` and `B` are sufficiently
identical for OI-FITS data.

"""
@inline is_same(A::T, B::T, sym::Symbol) where {T<:OIDataBlock} =
    is_same(A, B, Val(sym))

is_same(A::T, B::T, ::Val{S}) where {T<:OIDataBlock,S} = begin
    if (defined = isdefined(A, S)) != isdefined(B, S)
        return false
    elseif defined
        return is_same(getfield(A, S), getfield(B, S))
    else
        return true
    end
end

"""
    OIFITS.fix_name(str)

yields string `str` converted to uppercase letters and with trailing spaces
removed. This is to follow FITS conventions that names comparisons should be
done ignoring the case of letters and trailing spaces.

"""
fix_name(str::AbstractString) = uppercase(rstrip(str))

# Extend A[name] syntax to get an OI_ARRAY, OI_WAVELENGTH, or OI_CORR by its
# name (according to FITS conventions). The syntax keys(String,A) can be used
# to get an iterator over the string keys of A.
for (type, name) in ((:OI_ARRAY,      :arrname),
                     (:OI_WAVELENGTH, :insname),
                     (:OI_CORR,       :corrname))
    @eval begin
        getindex(A::AbstractArray{<:$type}, key::AbstractString) = begin
            @inbounds for i in eachindex(A)
                if is_same(A[i].$name, key)
                    return A[i]
                end
            end
            throw(KeyError(fix_name(key)))
        end
        keys(::Type{String}, A::AbstractArray{<:$type}) =
            Iterators.map(x -> fix_name(x.$name), A)
        haskey(A::AbstractArray{<:$type}, key::AbstractString) = begin
            @inbounds for i in eachindex(A)
                if is_same(A[i].$name, key)
                    return true
                end
            end
            return false
        end
        get(A::AbstractArray{<:$type}, key::AbstractString, def) = begin
            @inbounds for i in eachindex(A)
                if is_same(A[i].$name, key)
                    return A[i]
                end
            end
            return def
        end
    end
end

"""
    OIFITS.extname(arg)

yields the name of the FITS extension for argument `arg` which can be an
OI-FITS datablock instance or type or a FITS table HDU. In this latter case, an
empty string is returned if the `EXTNAME` keyword is not found and the
extension name is converted to upper case letters and trailing spaces discarded
otherwise. This method is not exported, the same result for a data-block
instance `db` is obtained by `db.extname`.

For a data-block instance or type `db`, an optional first argument, `Symbol` or
`String`, may be used to specify the type of the result:

    OIFITS.extname(Symbol, db) -> ext::Symbol
    OIFITS.extname(String, db) -> ext::String

"""
extname(hdu::FitsTableHDU) = fix_name(read_keyword(String, hdu, "EXTNAME", ""))
extname(db::OIDataBlock) = extname(typeof(db))
extname(T::Type{<:OIDataBlock}) = extname(String, T)
extname(S::Type{<:Union{String,Symbol}}, db::OIDataBlock) =
    extname(S, typeof(db))

for (N, T) in ((:OI_TARGET,     OI_TARGET    ),
               (:OI_ARRAY,      OI_ARRAY     ),
               (:OI_WAVELENGTH, OI_WAVELENGTH),
               (:OI_CORR,       OI_CORR      ),
               (:OI_VIS,        OI_VIS       ),
               (:OI_VIS2,       OI_VIS2      ),
               (:OI_T3,         OI_T3        ),
               (:OI_FLUX,       OI_FLUX      ),
               (:OI_INSPOL,     OI_INSPOL    ),)
    @eval begin
        extname(::Type{String}, ::Type{<:$T}) = $(String(N))
        extname(::Type{Symbol}, ::Type{<:$T}) = $(QuoteNode(N))
    end
end

"""
    OIFITS.get_format(db,  revn=db.revn; throw_errors=false)
    OIFITS.get_format(T,   revn;         throw_errors=false)
    OIFITS.get_format(ext, revn;         throw_errors=false)

all yield OI-FITS definitions for data-block `db`, of data-block type `T`, or
of OI-FITS extension `ext` (a string or a symbol).

"""
get_format(::Type{T}, revn::Integer; kwds...) where {T<:OIDataBlock} =
    get_format(extname(Symbol, T), revn; kwds...)
get_format(db::OIDataBlock, revn::Integer=db.revn; kwds...) =
    get_format(extname(Symbol, db), revn; kwds...)
get_format(extname::Union{Symbol,AbstractString}, revn::Integer; kwds...) =
    get_format(Symbol(extname), revn; kwds...)

#------------------------------------------------------------------------------
# PROPERTIES OF DATA-BLOCKS
#
# We use the `db.key` syntax to give access to other properties (e.g.,
# `extname`, `name`, `eff_wave`, or `eff_band`) than the fields of the
# data-block `db`.
#
# To allow for efficient handling of the `db.key` syntax, we wrap the `key`
# symbol as a `Val(key)` so that each `getproperty` or `setproperty!` method
# can be specific to the value of `key`.

# _properties(T) yields public properties for type `T`.
for T in (OI_TARGET, OI_ARRAY, OI_WAVELENGTH, OI_CORR,
          OI_VIS, OI_VIS2, OI_T3, OI_FLUX, OI_INSPOL)
    isconcretetype(T) || continue
    extra_fields = if T <: Union{OI_ARRAY,OI_WAVELENGTH,OI_CORR}
        (:extname, :name,)
    elseif T <: Union{OI_VIS,OI_VIS2,OI_T3,OI_FLUX}
        (:extname, :eff_wave, :eff_band,)
    else
        (:extname,)
    end
    @eval _properties(::Type{$T}) = $((fieldnames(T)..., extra_fields...))
end

propertynames(db::OIDataBlock) = _properties(typeof(db))

getproperty(db::OIDataBlock, sym::Symbol) = getproperty(db, Val(sym))
getproperty(db::OIDataBlock, ::Val{S}) where {S} = getfield(db, S)
getproperty(db::OIDataBlock, ::Val{:extname}) = extname(typeof(db))

setproperty!(db::OIDataBlock, sym::Symbol, val) =
    setproperty!(db, Val(sym), val)

setproperty!(db::T, ::Val{S}, val) where {S,T<:OIDataBlock} = begin
    try
        setfield!(db, S, convert(fieldtype(T, S), val))
    catch ex
        rethrow_convert_field_error(ex, T, sym, val)
    end
end

for (type, name) in ((:OI_ARRAY,      :arrname),
                     (:OI_WAVELENGTH, :insname),
                     (:OI_CORR,       :corrname))
    @eval begin
        getproperty(db::$type, ::Val{:name}) = db.$name
        setproperty!(db::$type, ::Val{:name}, val) = (db.$name = val)
    end
end

# NOTE: Directly calling `getfield` is necessary to optimize the indirection
#       (speedup: 1.4ns instead of 370ns).
_instr(db::Union{OI_VIS,OI_VIS2,OI_T3,OI_FLUX}) = getfield(db, :instr)
_instr(db::OI_WAVELENGTH) = db
for type in (:OI_VIS, :OI_VIS2, :OI_T3, :OI_FLUX)
    @eval begin
        getproperty(db::$type, ::Val{:eff_wave}) =
            _instr(db).eff_wave
        getproperty(db::$type, ::Val{:eff_band}) =
            _instr(db).eff_band
        setproperty!(db::$type, ::Val{:eff_wave}, val) =
            _instr(db).eff_wave = val
        setproperty!(db::$type, ::Val{:eff_band}, val) =
            _instr(db).eff_band = val
    end
end

@noinline function rethrow_convert_field_error(ex::Exception, T::Type,
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

#------------------------------------------------------------------------------
# PROPERTIES OF DATA-SETS
#
# Here, we want to make the "private" fields of data-sets (the dictionaries)
# "hidden" to the user. They can only be retrieved by getfield(ds,sym).

# Public properties for data-sets.
_properties(::Type{OIDataSet}) = (:target, :array, :instr, :correl,
                                  :vis, :vis2, :t3, :flux, :inspol)

# Only the "public" fields of data-sets are part of the properties.
propertynames(ds::OIDataSet) = _properties(OIDataSet)

# An OIDataSet is immutable so explicitely forbid writing a field is not
# necessary but the following is to provide a more explicit message.
@noinline setproperty!(db::OIDataSet, sym::Symbol, val) =
    error("attempt to ", (sym ∈ propertynames(db) ? "set read-only" :
                          "access non-existing or private"),
          " field `", sym, "` in ", nameof(typeof(db)), " instance")

getproperty(db::OIDataSet, sym::Symbol) = getproperty(db, Val(sym))

# By default, accessing any field by the dot notation is forbidden.
getproperty(db::OIDataSet, ::Val{S}) where {S} =
    error("attempt to access non-existing or private field `", S,
          "` in ", nameof(typeof(db)), " instance")

# Each public field is explicitely allowed.
for sym in _properties(OIDataSet)
    quoted_sym = QuoteNode(sym)
    @eval begin
        getproperty(db::OIDataSet, ::Val{$quoted_sym}) =
            getfield(db, $quoted_sym)
    end
end

#------------------------------------------------------------------------------
# COPY DATA-BLOCKS AND DATA-SETS

# Copy the outer structure of an OI-FITS data-block.
function copy(A::T) where {T<:OIDataBlock}
    B = T(undef)
    for sym in fieldnames(T)
        if isdefined(A, sym)
            setfield!(B, sym, getfield(A, sym))
        end
    end
    return B
end

function copy(A::OIDataSet)
    # Create new instance and copy contents. Including the dictionary for
    # re-writing target identifiers. The copy is thus in the same "state" as
    # the original. Note that using `map` is ~ 3.6 times slower here.
    B = OIDataSet()
    _copy!(B, A, :target)
    _copy!(B, A, :target_dict)
    _copy!(B, A, :target_id_map)
    _copy!(B, A, :array)
    _copy!(B, A, :array_dict)
    _copy!(B, A, :instr)
    _copy!(B, A, :instr_dict)
    _copy!(B, A, :correl)
    _copy!(B, A, :correl_dict)
    _copy!(B, A, :vis)
    _copy!(B, A, :vis2)
    _copy!(B, A, :t3)
    _copy!(B, A, :flux)
    _copy!(B, A, :inspol)
    return B
end

function _copy!(dest::Vector{T}, src::Vector{T}) where {T}
    if dest !== src && length(src) > 0
        copyto!(resize!(dest, length(src)), src)
    end
    return dest
end

function _copy!(dest::Dict{K,V}, src::Dict{K,V}) where {K,V}
    empty!(dest)
    for (key, val) in src
        dest[key] = val
    end
    return dest
end

function _copy!(dest::OI_TARGET, src::OI_TARGET)
    dest.revn = src.revn
    _copy!(dest.list, src.list)
    return dest
end

# @inline is needed here to make this as fast as possible.
# An alternative is to use Val(sym) (without @inline).
@inline _copy!(dest::T, src::T, sym::Symbol) where {T} =
    _copy!(getfield(dest, sym), getfield(src, sym))
#_copy!(dest::T, src::T, ::Val{S}) where {T,S} =
#    _copy!(getfield(dest, S), getfield(src, S))

#------------------------------------------------------------------------------
# PUSH DATA-BLOCKS IN DATA-SETS

# Create empty instance, then push all data-blocks which must be ordered.
OIDataSet(args::OIDataBlock...) = push!(OIDataSet(), args...)

function empty!(ds::OIDataSet)
    @inline _empty!(A, sym::Symbol) = empty!(getfield(A, sym))
    empty!(ds.target.list)
    ds.target.revn = 0
    _empty!(ds, :target_dict)
    _empty!(ds, :target_id_map)
    _empty!(ds, :array)
    _empty!(ds, :array_dict)
    _empty!(ds, :instr)
    _empty!(ds, :instr_dict)
    _empty!(ds, :correl)
    _empty!(ds, :correl_dict)
    _empty!(ds, :vis)
    _empty!(ds, :vis2)
    _empty!(ds, :t3)
    _empty!(ds, :flux)
    _empty!(ds, :inspol)
    return ds
end

for (type, pass) in ((:OI_TARGET,     1),
                     (:OI_ARRAY,      1),
                     (:OI_WAVELENGTH, 1),
                     (:OI_CORR,       1),
                     (:OI_VIS,        2),
                     (:OI_VIS2,       2),
                     (:OI_T3,         2),
                     (:OI_FLUX,       2),
                     (:OI_INSPOL,     2))
    @eval function _push!(ds::OIDataSet, db::$type, pass::Integer)
        if pass == $pass
            push!(ds, db)
        end
        nothing
    end
end

"""
    push!(ds::OIDataSet, db) -> ds

adds an OI-FITS data-block `db` in the data-set `ds` and returns `ds`. After
the operation, the contents of `db` is left unchanged but may be partially
shared by `ds`. Depending on type of `db`, different things can happen:

- If `db` is an `OI_TARGET` instance, the targets from `db` not already in the
  data-set `ds` are added to the list of targets in `ds` and the internal
  dictionary in `ds` mapping target identifiers in `db` to those in `ds` is
  reinitialized for subsequent data-blocks.

- If `db` is an `OI_ARRAY`, `OI_WAVELENGTH`, or `OI_CORR` instance, it is added
  to the corresponding list of data-blocks in the data-set `ds` unless an entry
  with a name matching that of `db` already exists in `ds`. In this latter
  case, nothing is done except checking that the two data-blocks with matching
  names have the same contents. This is to ensure the consistency of the
  data-set `ds`.

- If `db` is an `OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, or `OI_INSPOL`
  instance, it is added to the corresponding list of data-blocks in the
  data-set `ds` after having rewritten its target identifiers according to the
  mapping set by the last `OI_TARGET` pushed into the data-set `ds`. Add
  keyword `rewrite_target_id=false` to avoid rewritting target identifiers.

"""
function push!(ds::OIDataSet, db::OI_TARGET)
    # Check data-block.
    check(db)

    # Get current list of targets in data-set and corresponding mapping of
    # target names to identifiers.
    target_list = ds.target.list
    target_dict = getfield(ds, :target_dict)
    length(target_list) == length(target_dict) || error(
        "inconsistent list and dictionary of targets")

    # Get and reset dictionary to reindex target identifiers in this and
    # subsequent data-blocks.
    target_id_map = getfield(ds, :target_id_map)
    empty!(target_id_map)

    # Insert all targets from the data-block to the data-set.
    for tgt in db.list
        name = fix_name(tgt.target)
        id = Int(tgt.target_id)
        if haskey(target_id_map, id)
            is_same(target_list[target_id_map[id]], tgt) || error(
                "duplicate target identifier ", id, " in ", db.extname,
                " but for different targets")
        end
        if haskey(target_dict, name)
            # The same target already existed in destination. Check that they
            # have the same parameters.
            index = target_dict[name]
            is_same(target_list[index], tgt) || error(
                "target \"", tgt.target, "\" in ", db.extname, " has ",
                "different parameters than a previous one matching this name")
        else
            # This target is new.
            index = length(target_list) + 1
            push!(target_list, OITargetEntry(tgt; target_id = index))
            target_dict[name] = index
        end
        target_id_map[id] = index
    end

    # Finally fix revision number.
    ds.target.revn = max(ds.target.revn, db.revn)
    return ds
end

function push!(ds::OIDataSet,
               db::T) where {T<:Union{OI_ARRAY,OI_WAVELENGTH,OI_CORR}}
    check(db)
    dict = _dict(T, ds)
    list = _list(T, ds)
    name = fix_name(db.name)
    if haskey(dict, name)
        # Make sure the two data-blocks with the same name are
        # identical.
        assert_identical(db, list[dict[name]])
    else
        # Push new data-block in destination.
        dict[name] = length(push!(list, db))
    end
    return ds
end

function push!(ds::OIDataSet, db::T;
               rewrite_target_id::Bool = true) where {
                   T<:Union{OI_VIS,OI_VIS2,OI_T3,OI_FLUX,OI_INSPOL}}
    # Define private helper.
    function find_depencency(dest_list::Vector{T},
                             dest_dict::Dict{String,Int},
                             dep::AbstractString,
                             optional::Bool = false) where {T<:OIDataBlock}
        name = fix_name(dep)
        if haskey(dest_dict, name)
            return dest_list[dest_dict[name]]
        else
            optional || error(
                "no extensions ", extname(T), " match name \"", name, "\"")
            return nothing
        end
    end

    check(db)
    new_db = copy(db) # make a copy to avoid side effects
    if rewrite_target_id
        new_db.target_id = rewrite_indices(db.target_id,
                                           getfield(ds, :target_id_map))
    else
        id_min, id_max = extrema(db.target_id)
        if id_min < 1 || id_max > length(db.target)
            error("out of range target identifier(s) (", id_min, ":", id_max,
                  " ⊈ 1:", length(db.target), ")")
        end
    end
    if isdefined(new_db, :arrname)
        let dep = find_depencency(ds.array, getfield(ds, :array_dict),
                                  new_db.arrname, true)
            if dep !== nothing
                new_db.array = dep
            end
        end
    end
    if isdefined(new_db, :insname)
        new_db.instr = find_depencency(
            ds.instr, getfield(ds, :instr_dict), new_db.insname)
    end
    if isdefined(new_db, :corrname)
        new_db.correl = find_depencency(
            ds.correl, getfield(ds, :correl_dict), new_db.corrname)
    end
    push!(_list(T, ds), new_db)
    return ds
end

# `_list(T,ds)` yields the list of entries of type `T` stored in `ds`.
_list(::Type{OI_ARRAY},      ds::OIDataSet) = ds.array
_list(::Type{OI_WAVELENGTH}, ds::OIDataSet) = ds.instr
_list(::Type{OI_CORR},       ds::OIDataSet) = ds.correl
_list(::Type{OI_VIS},        ds::OIDataSet) = ds.vis
_list(::Type{OI_VIS2},       ds::OIDataSet) = ds.vis2
_list(::Type{OI_T3},         ds::OIDataSet) = ds.t3
_list(::Type{OI_FLUX},       ds::OIDataSet) = ds.flux
_list(::Type{OI_INSPOL},     ds::OIDataSet) = ds.inspol
#_list(::Type{OITargetEntry},ds::OIDataSet) = ds.target.list

# `_dict(T,ds)` yields the dictionary associated to entries of type `T` stored
# by the OI-FITS data-set `ds`.
_dict(::Type{OI_ARRAY},      ds::OIDataSet) = getfield(ds, :array_dict)
_dict(::Type{OI_WAVELENGTH}, ds::OIDataSet) = getfield(ds, :instr_dict)
_dict(::Type{OI_CORR},       ds::OIDataSet) = getfield(ds, :correl_dict)
#_dict(::Type{OITargetEntry},ds::OIDataSet) = getfield(ds, :target_dict)

"""
    OIFITS.rewrite_indices(inds, dict) -> inds′

yields array of indices similar to `inds` and rewritten according to the
dictionary `dict` used as a mapping from integers to integers. If the indices
are left unchanged, the input array is returned.

"""
function rewrite_indices(inds::AbstractArray{<:Integer,N},
                         dict::Dict{<:Integer,<:Integer}) where {N}
    new_inds = similar(inds)
    unchanged = true
    for i in eachindex(new_inds, inds)
        idx = inds[i]
        new_idx = dict[idx]
        new_inds[i] = new_idx
        unchanged &= (new_idx == idx)
    end
    return (unchanged ? inds : new_inds)
end

"""
    OIFITS.assert_identical(A, B)

throws an exception if *named* data-blocks `A` and `B` are not identical. A
*named* data-block is an `OI_ARRAY`, `OI_WAVELENGTH`, or `OI_CORR` extension
that can be uniquely identifed by its name.

The variant

    OIFITS.assert_identical(A, B, sym)

is to assert that fields `A.sym` and `B.sym` are identical, it shall only be
called if `A` and `B` have the same type and the same name.

""" assert_identical

function assert_identical(A::T, B::T) where {T<:Union{OI_ARRAY,OI_WAVELENGTH,OI_CORR}}
    if !(A === B)
        is_same(A.name, B.name) || throw_assertion_error(
            "two ", extname(T), " data-blocks named \"", A.name,
            "\" and \"", B.name, "\" cannot be identical")
        assert_identical(A, B, *)
    end
end

function assert_identical(A::OI_ARRAY, B::OI_ARRAY, ::typeof(*))
    assert_identical(A, B, :revn)
    assert_identical(A, B, :frame)
    assert_identical(A, B, :arrayx)
    assert_identical(A, B, :arrayy)
    assert_identical(A, B, :arrayz)
    assert_identical(A, B, :tel_name)
    assert_identical(A, B, :sta_name)
    assert_identical(A, B, :sta_index)
    assert_identical(A, B, :diameter)
    assert_identical(A, B, :staxyz)
    assert_identical(A, B, :fov)
    assert_identical(A, B, :fovtype)
end

function assert_identical(A::OI_WAVELENGTH, B::OI_WAVELENGTH, ::typeof(*))
    assert_identical(A, B, :revn)
    assert_identical(A, B, :eff_wave)
    assert_identical(A, B, :eff_band)
end

function assert_identical(A::OI_CORR, B::OI_CORR, ::typeof(*))
    assert_identical(A, B, :revn)
    assert_identical(A, B, :ndata)
    assert_identical(A, B, :iindx)
    assert_identical(A, B, :jindx)
    assert_identical(A, B, :corr)
end

@inline assert_identical(A::T, B::T, sym::Symbol) where {T<:NamedDataBlock} =
    assert_identical(A, B, Val(sym))

assert_identical(A::T, B::T, val::Val{sym}) where {T<:NamedDataBlock,sym} =
    is_same(A, B, val) || throw_assertion_error(
        "two ", extname(T), " extensions both named \"", fix_name(A.name),
        " have different field `", sym, "`")

throw_assertion_error(msg::AbstractString) = throw(AssertionError(msg))
@noinline throw_assertion_error(args...) =
    throw_assertion_error(string(args...))

"""
    OIFITS.check(db) -> (nchns, nfrms)

checks that mandatory fields of data-block `db` are defined and that defined
fields have correct sizes. The number of spectral channels `nchns` and of
temporal frames `nfrms` are returned (with a value of `-1` if not relevant for
the type of `db`).

    OIFITS.check(ds)

checks all data-blocks of data-set `ds`.

""" check

function check(ds::OIDataSet)
    check(ds.target)
    for db in ds.array
        check(db)
    end
    for db in ds.instr
        check(db)
    end
    for db in ds.correl
        check(db)
    end
    for db in ds.vis
        check(db)
    end
    for db in ds.vis2
        check(db)
    end
    for db in ds.t3
        check(db)
    end
    for db in ds.flux
        check(db)
    end
    for db in ds.inspol
        check(db)
    end
    nothing
end

function check(db::OIDataBlock)
    nchns = if isdefined(db, :instr)
        length(getfield(db, :instr).eff_wave)
    else
        -1 # number of spectral channels not yet known
    end
    nfrms = -1 # number of temporal frames not yet known
    for spec in get_format(db; throw_errors=true)
        nchns, nfrms = _check(db, Val(spec.symb), spec, nchns, nfrms)
    end
    return nchns, nfrms
end

function check(db::OI_TARGET)
    isdefined(db, :list) || throw_undefined_field(db, :list)
    return -1, -1
end

function _check(db::OIDataBlock, ::Val{S}, spec::FieldDefinition,
                nchns::Int, nfrms::Int) where {S}
    if !isdefined(db, S)
        spec.optional || throw_undefined_field(db, S)
    else
        val = getfield(db, S)
        dims = (isa(val, AbstractString) ? () : size(val))
        ndims_spec = ndims(spec)
        length(dims) == ndims_spec || error(
            "field `", S, "` of ", db.extname, " must be a ",
            (ndims_spec == 0 ? "scalar" : "$ndims_spec-D array"))
        rank = spec.rank
        if ndims_spec ≥ 1
            # Last dimension is number of temporal frames.
            val_nfrms = dims[end]
            if nfrms == -1
                nfrms = val_nfrms
            else
                nfrms == val_nfrms || error(
                    "incompatible number of temporal frames in field `", S,
                    "` of ", db.extname, " (got ", val_nfrms, ", should be ",
                    nfrms, ")")
            end
        end
        if rank > 1 &&  ndims_spec == 2 && spec.type !== :A
            dims[1] == rank || error(
                "invalid 1st dimension in field `", S,
                "` of ", db.extname, " (got ", dims[1], ", should be ",
                rank, ")")
        elseif rank < 0
            # All leading dimensions are equal to the number of spectral
            # channels.
            val_nchns = dims[1]
            for i in 2:length(dims)-1
                dims[i] == val_nchns || error(
                    "leading dimensions in field `", S, "` of ", db.extname,
                    " must be all equal to the numeber of spectral channels,",
                    " size is ", dims)
            end
            if nchns == -1
                nchns = val_nchns
            else
                nchns == val_nchns || error(
                    "incompatible number of spectral channels in field `", S,
                    "` of ", db.extname, " (got ", val_nchns, ", should be ",
                    nchns, ")")
            end
        end
    end
    return nchns, nfrms
end

const FieldName = Union{AbstractString,Symbol}

@noinline throw_undefined_field(name::FieldName) =
    throw(ErrorException("undefined field `$name`"))

@noinline throw_undefined_field(obj::T, name::FieldName) where {T} =
    throw(ErrorException(string(
        "undefined field `", name, "` in object of type `", nameof(T), "`")))

@noinline throw_undefined_field(obj::OIDataBlock, name::FieldName) =
    throw(ErrorException(string(
        "mandatory field `", name, "` undefined in `", obj.extname,
        "` data-block")))


#------------------------------------------------------------------------------
# MERGE DATA-SETS

merge(A::OIDataSet) = copy(A)
merge(A::OIDataSet, args::OIDataSet...) = merge!(OIDataSet(), A, args...)

merge!(A::OIDataSet) = A
merge!(A::OIDataSet, others::OIDataSet...) = begin
    for B in others
        merge!(A, B)
    end
    return A
end

function merge!(A::OIDataSet, B::OIDataSet)
    # First push targets and dependencies.
    push!(A, B.target)
    for db in B.array
        push!(A, db)
    end
    for db in B.instr
        push!(A, db)
    end
    for db in B.correl
        push!(A, db)
    end

    # Then push all other data-blocks.
    for db in B.vis
        push!(A, db)
    end
    for db in B.vis2
        push!(A, db)
    end
    for db in B.t3
        push!(A, db)
    end
    for db in B.flux
        push!(A, db)
    end
    for db in B.inspol
        push!(A, db)
    end

    # Return updated data-set.
    return A
end

#------------------------------------------------------------------------------
# OI_TARGET DATA-BLOCKS

"""
    OITargetEntry([def::OITargetEntry]; kwds...)

yields an entry of the `OI_TARGET` table in an OI-FITS file. All fields are
specified by keywords. If template entry `def` is specified, keywords have
default values taken from `def`; otherwise all keywords are mandatory but
`category` which is assumed to be `""` if unspecified.

"""
function OITargetEntry(;
                       target_id ::Integer,
                       target    ::AbstractString,
                       raep0     ::AbstractFloat,
                       decep0    ::AbstractFloat,
                       equinox   ::AbstractFloat,
                       ra_err    ::AbstractFloat,
                       dec_err   ::AbstractFloat,
                       sysvel    ::AbstractFloat,
                       veltyp    ::AbstractString,
                       veldef    ::AbstractString,
                       pmra      ::AbstractFloat,
                       pmdec     ::AbstractFloat,
                       pmra_err  ::AbstractFloat,
                       pmdec_err ::AbstractFloat,
                       parallax  ::AbstractFloat,
                       para_err  ::AbstractFloat,
                       spectyp   ::AbstractString,
                       category  ::AbstractString = empty_string)
    return OITargetEntry(target_id,
                         target,
                         raep0,
                         decep0,
                         equinox,
                         ra_err,
                         dec_err,
                         sysvel,
                         veltyp,
                         veldef,
                         pmra,
                         pmdec,
                         pmra_err,
                         pmdec_err,
                         parallax,
                         para_err,
                         spectyp,
                         category)
end

function OITargetEntry(def::OITargetEntry;
                       target_id ::Integer        = def.target_id,
                       target    ::AbstractString = def.target,
                       raep0     ::AbstractFloat  = def.raep0,
                       decep0    ::AbstractFloat  = def.decep0,
                       equinox   ::AbstractFloat  = def.equinox,
                       ra_err    ::AbstractFloat  = def.ra_err,
                       dec_err   ::AbstractFloat  = def.dec_err,
                       sysvel    ::AbstractFloat  = def.sysvel,
                       veltyp    ::AbstractString = def.veltyp,
                       veldef    ::AbstractString = def.veldef,
                       pmra      ::AbstractFloat  = def.pmra,
                       pmdec     ::AbstractFloat  = def.pmdec,
                       pmra_err  ::AbstractFloat  = def.pmra_err,
                       pmdec_err ::AbstractFloat  = def.pmdec_err,
                       parallax  ::AbstractFloat  = def.parallax,
                       para_err  ::AbstractFloat  = def.para_err,
                       spectyp   ::AbstractString = def.spectyp,
                       category  ::AbstractString = def.category)
    return OITargetEntry(target_id,
                         target,
                         raep0,
                         decep0,
                         equinox,
                         ra_err,
                         dec_err,
                         sysvel,
                         veltyp,
                         veldef,
                         pmra,
                         pmdec,
                         pmra_err,
                         pmdec_err,
                         parallax,
                         para_err,
                         spectyp,
                         category)
end

"""
    OI_TARGET(lst=OITargetEntry[]; revn=0)

yields an `OI_TARGET` data-block. Optional argument `lst` is a vector of
`OITargetEntry` specifying the targets (none by default). Keyword `revn`
specifies the revision number.

An instance, say `db`, of `OI_TARGET` has the following properties:

    db.extname # OI-FITS extension name
    db.list    # list of targets
    db.revn    # revision number

and can be used as an iterable or as an array to access the targets
individually by their index or by their name:

    length(db) # the number of targets
    db[i]      # the i-th target
    db[key]    # the target whose name matches string `key`

    for tgt in db; ...; end # loop over all targets

Note that, when an `OI_TARGET` instance is pushed in a data-set, target
identifiers (field `target_id`) are automatically rewritten to be identical to
the index in the list of targets of the data-set.

Standard methods `get` and `haskey` work as expected and according to the type
(integer or string) of the key. For the `keys` method, the default is to return
an iterator over the target names, but the type of the expected keys can be
specified:

    keys(db)          # iterator over target names
    keys(String, db)  # idem
    keys(Integer, db) # iterator over target indices
    keys(Int, db)     # idem

Call [`OIFITS.get_column`](@ref) to retrieve a given target field for all
targets of `OI_TARGET` data-block in the form of a vector.

""" OI_TARGET

# Make OI_TARGET instances iterable and behave more or less like vectors
# (with properties).
ndims(db::OI_TARGET)   =  ndims(db.list)
size(db::OI_TARGET)    =   size(db.list)
size(db::OI_TARGET, i) =   size(db.list, i)
axes(db::OI_TARGET)    =   axes(db.list)
axes(db::OI_TARGET, i) =   axes(db.list, i)
length(db::OI_TARGET)  = length(db.list)

IndexStyle(db::OI_TARGET) = IndexStyle(OI_TARGET)
IndexStyle(::Type{OI_TARGET}) = IndexStyle(fieldtype(OI_TARGET, :list))

eltype(db::OI_TARGET) = eltype(OI_TARGET)
eltype(::Type{OI_TARGET}) = OITargetEntry

@inline @propagate_inbounds getindex(db::OI_TARGET, i::Integer) = db.list[i]

getindex(A::AbstractVector{OITargetEntry}, name::AbstractString) = begin
    x = get(A, name, undef)
    x === undef && error("no targets named \"", fix_name(name), "\" found")
    return x
end

getindex(db::OI_TARGET, name::AbstractString) = begin
    x = get(db.list, name, undef)
    x === undef && error( "no targets named \"", fix_name(name),
                          "\" found in ", db.extname, "data-block")
    return x
end

firstindex(db::OI_TARGET) = firstindex(db.list)
lastindex(db::OI_TARGET) = lastindex(db.list)
eachindex(db::OI_TARGET) = eachindex(db.list)

function iterate(db::OI_TARGET,
                 state::Tuple{Int,Int} = (firstindex(db),
                                          lastindex(db)))
    i, n = state
    if i > n
        return nothing
    else
        return db.list[i], (i + 1, n)
    end
end

haskey(db::OI_TARGET, i::Integer) = (1 ≤ i ≤ length(db))

haskey(db::OI_TARGET, key::AbstractString) = haskey(db.list, key)

haskey(A::AbstractVector{OITargetEntry}, key::AbstractString) = begin
    for tgt in A
        if is_same(tgt.target, key)
            return true
        end
    end
    return false
end

get(db::OI_TARGET, i::Integer, def) =
    if 1 ≤ i ≤ length(db)
        return @inbounds db[i]
    else
        return def
    end

get(db::OI_TARGET, key::AbstractString, def) = get(db.list, key, def)

get(A::AbstractVector{OITargetEntry}, key::AbstractString, def) = begin
    for tgt in A
        if is_same(tgt.target, key)
            return tgt
        end
    end
    return def
end

# Note: this is the default behavior.
values(db::OI_TARGET) = db

keys(db::OI_TARGET) = Iterators.map(x -> fix_name(x.target), db)
keys(::Type{String}, db::OI_TARGET) = keys(db)
keys(::Type{Int}, db::OI_TARGET) = eachindex(db)
keys(::Type{Integer}, db::OI_TARGET) = eachindex(db)

"""
    OIFITS.get_column([T,] db::OI_TARGET, col)

yields the column `col` of an OI-FITS data-block `db`. Column is identified by
`col` which is either `sym` or `Val(sym)` where `sym` is the symbolic name of
the corresponding field in `OITargetEntry`.

Optional argument `T` is to specify the element type of the returned array.

""" get_column

@inline get_column(db::OI_TARGET, sym::Symbol) =
    get_column(db, Val(sym))

@inline get_column(T::Type, db::OI_TARGET, sym::Symbol) =
    get_column(T, db, Val(sym))

get_column(db::OI_TARGET, ::Val{sym}) where {sym} =
    get_column(fieldtype(OITargetEntry, sym), db, sym)

function get_column(::Type{T}, db::OI_TARGET, ::Val{sym}) where {T,sym}
    len = length(db)
    col = Vector{T}(undef, len)
    for i in 1:len
        col[i] = getfield(db[i], sym)
    end
    return col
end
