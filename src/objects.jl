#
# objects.jl --
#
# Implement most part of OI-FITS objects behavior.
#
#------------------------------------------------------------------------------

"""
    OIFITS.contents(obj)

yields the contents of data-block or master `obj`.

"""
contents(obj::OIDataBlock) = getfield(obj, :contents)
contents(obj::OIMaster) = getfield(obj, :all)

"""
    OIFITS.is_attached(db)

yields whether data-block `db` is attached to a master structure.

"""
is_attached(db::OIDataBlock) = getfield(db, :attached)

#------------------------------------------------------------------------------
# INDEXING

# Extend getproperty/setproperty! to implement consistent db.field syntax.  The
# behavior depends on the type of data-block.
Base.getproperty(obj::OIDataBlock, sym::Symbol) =
    (sym == :extname ? get_extname(obj) :
     getindex(obj, sym))
Base.setproperty!(obj::OIDataBlock, sym::Symbol, val) =
    (sym == :extname ? read_only_field(sym) :
     setindex!(obj, val, sym))

Base.getproperty(obj::OIPolarization, sym::Symbol) =
    (sym == :array   ? getfield(obj, :arr) :
     sym == :extname ? get_extname(obj) :
     getindex(obj, sym))
Base.setproperty!(obj::OIPolarization, sym::Symbol, val) =
    (sym == :array   ? setfield!(obj, :arr, val) :
     sym == :extname ? read_only_field(sym) :
     setindex!(obj, sym, val))

Base.getproperty(obj::Union{OIVis,OIVis2,OIT3,OISpectrum}, sym::Symbol) =
    (sym == :array   ? getfield(obj, :arr) :
     sym == :instr   ? getfield(obj, :ins) :
     sym == :correl  ? getfield(obj, :corr) :
     sym == :extname ? get_extname(obj) :
     getindex(obj, sym))
Base.setproperty!(obj::Union{OIVis,OIVis2,OIT3,OISpectrum}, sym::Symbol, val) =
    (sym == :array   ? setfield!(obj, :arr, val) :
     sym == :instr   ? setfield!(obj, :ins, val) :
     sym == :correl  ? setfield!(obj, :corr, val) :
     sym == :extname ? read_only_field(sym) :
     setindex!(obj, sym, val))

Base.getproperty(obj::OIMaster, sym::Symbol) =
    (sym == :array   ? getfield(obj, :arr) :
     sym == :instr   ? getfield(obj, :ins) :
     sym == :correl  ? getfield(obj, :corr) :
     sym == :target  ? getfield(obj, :tgt) :
     unknown_field(sym))
Base.setproperty!(obj::OIMaster, sym::Symbol, val) =
    read_only_field(sym)

@noinline read_only_field(name::Union{AbstractString,Symbol}) =
    error("field `", name, "` is read-only")

@noinline unknown_field(name::Union{AbstractString,Symbol}) =
    error("unknown field `", name, "`")

# Make OIMaster usable as an interator and as a vector of data-blocks.
Base.length(obj::OIMaster) = length(contents(obj))
Base.size(obj::OIMaster) = (length(obj),)
Base.IndexStyle(::Type{<:OIMaster}) = IndexLinear()
Base.getindex(obj::OIMaster, i::Integer) =
    getindex(contents(obj), i)
Base.iterate(obj::OIMaster) = iterate(contents(obj))
Base.iterate(obj::OIMaster, state) = iterate(contents(obj), state)
# FIXME: Base.IteratorSize(::Type{<:OIMaster}) = IteratorSize(Vector{OIDataBlock})
# FIXME: Base.IteratorEltype(::Type{<:OIMaster}) = IteratorEltype(Vector{OIDataBlock})

# OIDataBlock can be indexed by the name (either as a string or as a
# symbol) of the field and can be used as an iterator.

Base.length(obj::OIDataBlock) = length(contents(obj))

Base.iterate(obj::OIDataBlock) = iterate(contents(obj))
Base.iterate(obj::OIDataBlock, state) = iterate(contents(obj), state)

Base.IteratorSize(::Type{<:OIDataBlock}) = IteratorSize(OIContents)
Base.IteratorEltype(::Type{<:OIDataBlock}) = IteratorEltype(OIContents)

Base.getindex(obj::OIDataBlock, key::Symbol) = get(obj, key, nothing)
Base.getindex(obj::OIDataBlock, key::AbstractString) =
    getindex(obj, Symbol(key))

Base.setindex!(obj::OIDataBlock, val, key::Symbol) =
    setindex!(contents(obj), val, key)
Base.setindex!(obj::OIDataBlock, val, key::AbstractString) =
    setindex!(obj, val, Symbol(key))

Base.haskey(obj::OIDataBlock, key::Symbol) = haskey(contents(obj), key)
Base.haskey(obj::OIDataBlock, key::AbstractString) = haskey(obj, Symbol(key))
Base.haskey(obj::OIDataBlock, key) = false

Base.get(obj::OIDataBlock, key::Symbol, def) = get(contents(obj), key, def)
Base.get(obj::OIDataBlock, key::AbstractString, def) =
    get(obj, Symbol(key), def)

Base.get(f::Function, obj::OIDataBlock, key::Symbol) = get(f, contents(obj), key)
Base.get(f::Function, obj::OIDataBlock, key::AbstractString) =
    get(f, obj, Symbol(key))

Base.get!(obj::OIDataBlock, key::Symbol, def) = get!(contents(obj), key, def)
Base.get!(obj::OIDataBlock, key::AbstractString, def) =
    get!(obj, Symbol(key), def)

Base.get!(f::Function, obj::OIDataBlock, key::Symbol) =
    get!(f, contents(obj), key)
Base.get!(f::Function, obj::OIDataBlock, key::AbstractString) =
    get!(f, obj, Symbol(key))

Base.keys(obj::OIDataBlock) = keys(contents(obj))
Base.values(obj::OIDataBlock) = values(contents(obj))

Base.getkey(obj::OIDataBlock, key::Symbol, def) =
    getkey(contents(obj), key, def)
Base.getkey(obj::OIDataBlock, key::AbstractString, def) =
    getkey(obj, Symbol(key), def)

Base.delete!(obj::OIDataBlock, key::Symbol) = delete!(contents(obj), key)
Base.delete!(obj::OIDataBlock, key::AbstractString) = delete!(obj, Symbol(key))

Base.pop!(obj::OIDataBlock, key::Symbol) = pop!(obj, key, nothing)
Base.pop!(obj::OIDataBlock, key::AbstractString) = pop!(obj, Symbol(key))

Base.pop!(obj::OIDataBlock, key::Symbol, def) = pop!(contents(obj), key, def)
Base.pop!(obj::OIDataBlock, key::AbstractString, def) =
    pop!(obj, Symbol(key), def)

Base.pairs(obj::OIDataBlock) = pairs(contents(obj))
Base.pairs(sty::IndexStyle, obj::OIDataBlock) = pairs(sty, contents(obj))

Base.merge(obj::OIDataBlock, others::AbstractDict...) =
    OIDataBlock(merge(contents(obj), others...), nothing)
Base.merge!(obj::OIDataBlock, others::AbstractDict...) =
    (merge!(contents(obj), others...); obj)

Base.sizehint!(obj::OIDataBlock, n) = (sizehint!(contents(obj), n); obj)

Base.keytype(::Type{<:OIDataBlock}) = keytype(OIContents)

Base.valtype(::Type{<:OIDataBlock}) = valtype(OIContents)

#------------------------------------------------------------------------------

# Known OI-FITS extension names.
const EXTNAMES = ("OI_TARGET", "OI_WAVELENGTH", "OI_ARRAY", "OI_VIS",
                  "OI_VIS2", "OI_T3", "OI_SPECTRUM", "OI_CORR", "OI_INSPOL")

"""
    OIFITS.get_datablock_type(extname) -> T

yields the data-block type `T <: OIDataBlock` associated to the OI-FITS
extension `extname`.

"""
get_datablock_type(extname::Symbol) =
    (extname === :OI_TARGET     ? OITarget       :
     extname === :OI_WAVELENGTH ? OIWavelength   :
     extname === :OI_ARRAY      ? OIArray        :
     extname === :OI_VIS        ? OIVis          :
     extname === :OI_VIS2       ? OIVis2         :
     extname === :OI_T3         ? OIT3           :
     extname === :OI_SPECTRUM   ? OISpectrum     :
     extname === :OI_CORR       ? OICorrelation  :
     extname === :OI_INSPOL     ? OIPolarization :
     bad_extname(extname))

get_datablock_type(extname::AbstractString) =
    (extname == "OI_TARGET"     ? OITarget       :
     extname == "OI_WAVELENGTH" ? OIWavelength   :
     extname == "OI_ARRAY"      ? OIArray        :
     extname == "OI_VIS"        ? OIVis          :
     extname == "OI_VIS2"       ? OIVis2         :
     extname == "OI_T3"         ? OIT3           :
     extname == "OI_SPECTRUM"   ? OISpectrum     :
     extname == "OI_CORR"       ? OICorrelation  :
     extname == "OI_INSPOL"     ? OIPolarization :
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
get_extname(::Type{OITarget})       = "OI_TARGET"
get_extname(::Type{OIWavelength})   = "OI_WAVELENGTH"
get_extname(::Type{OIArray})        = "OI_ARRAY"
get_extname(::Type{OIVis})          = "OI_VIS"
get_extname(::Type{OIVis2})         = "OI_VIS2"
get_extname(::Type{OIT3})           = "OI_T3"
get_extname(::Type{OISpectrum})     = "OI_SPECTRUM"
get_extname(::Type{OICorrelation})  = "OI_CORR"
get_extname(::Type{OIPolarization}) = "OI_INSPOL"
get_extname(hdr::FITSHeader) =
    (get_hdu_type(hdr) !== :binary_table ? "" :
     fix_name(get_string(hdr, "EXTNAME", "")))

function show(io::IO, db::OITarget)
    print(io, "OI_TARGET: ")
    if haskey(db, :target)
        tgt = db[:target]
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
    ntels = (haskey(db, :sta_index) ? length(db[:sta_index]) : 0)
    if haskey(db, :arrname)
        print(io, "arrname=\"", db[:arrname], "\", ", ntels, " telescope(s)")
    else
        print(io, ntels, " telescope(s)")
    end
end

function show(io::IO, db::OIWavelength)
    print(io, "OI_WAVELENGTH: ")
    if ! haskey(db, :insname)
        print(io, "<empty>")
    else
        print(io, "insname=\"", db[:insname], "\"")
        if haskey(db, :eff_wave)
            wave = db[:eff_wave]
            n = length(db[:eff_wave])
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

function show(io::IO, db::OIDataBlock,
              name::AbstractString, nwaves::Integer, ntimes::Integer)
    print(io, name, ": ")
    if nwaves > 0 && ntimes > 0
        print(io, nwaves*ntimes, " measurements in ", nwaves,
              " spectral channel(s) and ", ntimes, " exposure(s)")
    else
        print(io, name, "<empty>")
    end
end

function show(io::IO, db::OIVis)
    if haskey(db, :visamp)
        dims = size(db[:visamp])
        nwaves = dims[1]
        ntimes = dims[2]
    else
        nwaves = 0
        ntimes = 0
    end
    show(io, db, "OI_VIS", nwaves, ntimes)
end

function show(io::IO, db::OIVis2)
    if haskey(db, :vis2data)
        dims = size(db[:vis2data])
        nwaves = dims[1]
        ntimes = dims[2]
    else
        nwaves = 0
        ntimes = 0
    end
    show(io, db, "OI_VIS2", nwaves, ntimes)
end

function show(io::IO, db::OIT3)
    if haskey(db, :t3amp)
        dims = size(db[:t3amp])
        nwaves = dims[1]
        ntimes = dims[2]
    else
        nwaves = 0
        ntimes = 0
    end
    show(io, db, "OI_T3", nwaves, ntimes)
end

function show(io::IO, master::OIMaster)
    print(io, "OI_MASTER: (", length(master), " data-block(s))")
    for db in master
        print(io, "\n    ")
        show(io, db)
    end
end


##########################
# BUILDING OF DATABLOCKS #
##########################

# Build constructors for the various data-blocks types.
for (func, name) in ((:new_target,     "OI_TARGET"),
                     (:new_array,      "OI_ARRAY"),
                     (:new_wavelength, "OI_WAVELENGTH"),
                     (:new_vis,        "OI_VIS"),
                     (:new_vis2,       "OI_VIS2"),
                     (:new_t3,         "OI_T3"),
                     (:new_spectrum,   "OI_SPECTRUM"),
                     (:new_corr,       "OI_CORR"),
                     (:new_inspol,     "OI_INSPOL"))
    @eval begin
        function $func(master::Union{OIMaster, Nothing}=nothing;
                       revn::Integer = 2, kwds...)
            db = build_datablock($name, revn, kwds)
            if master != nothing
                push!(master, db)
            end
            return db
        end
    end
end

build_datablock(extname::AbstractString, revn::Integer, kwds) =
    build_datablock(to_string(extname), to_integer(revn), kwds)

function build_datablock(extname::String, revn::Int, kwds)
    def = Parser.get_definition(extname, revn)
    contents = OIContents()
    contents[:revn] = revn
    rows = -1         # number of rows in the OI-FITS table
    channels = -1     # number of spectral channels
    for (field, value) in kwds
        # Check whether this field exists.
        spec = get(def.spec, field, nothing)
        spec === nothing && error("OI-FITS extension ", extname,
                                  " has no field `", field, "`")

        # Check value type.
        if spec.type === :LOGICAL
            is_logical(value) || error("expecting boolean value for `", field,
                                       "` field of OI-FITS extension ",
                                       extname)
        elseif spec.type === :INTEGER
            is_integer(value) || error("expecting integer value for `", field,
                                       "` field of OI-FITS extension ",
                                       extname)
            value = to_integer(value)
        elseif spec.type === :FLOAT
            is_float(value) || error("expecting floating-point value for `",
                                     field, "` field of OI-FITS extension ",
                                     extname)
            value = to_float(value)
        elseif spec.type === :COMPLEX
            is_complex(value) || error("expecting complex value for `",
                                       field, "` field of OI-FITS extension ",
                                       extname)
            value = to_complex(value)
        elseif spec.type === :STRING
            is_string(value) || error("expecting string value for `", field,
                                      "` field of OI-FITS extension ", extname)
        else
            error("*** CORRUPTED WORKSPACE ***")
        end

        # Check value dimensions.
        dims = (isa(value, Array) ? size(value) : ())
        rank = length(dims)
        if spec.iskeyword
            rank == 0 || error("expecting a scalar value for `", field,
                               "` field in OI-FITS extension ", extname)
        else
            # Check array rank.
            mult = spec.multiplier
            maxrank = (spec.type === :STRING || mult == 1 ? 1 :
                       mult == -2 ? 3 : 2)
            rank > maxrank && error("bad number of dimensions for `", field,
                                    "` field of OI-FITS extension ", extname)

            # Fields may have up to 3 dimensions.  For now, the last dimension
            # is the number of rows.  FIXME: The ordering of dimensions must
            # be changed: the number of rows should be the first dimension.
            dim0 = (rank >= 1 ? dims[end] : 1)
            dim1 = (rank >= 2 ? dims[1] : 1)
            dim2 = (rank >= 3 ? dims[2] : 1)

            if rows == -1
                rows = dim0
            else
                rows == dim0 || error("incompatible number of rows for `",
                                      field, "` field of OI-FITS extension ",
                                      extname)
            end

            if mult < 0
                # Expecting an N-by-W or N-by-W-by-W array (N is the number of
                # rows and W the number of channels).
                (mult == -2 && dim1 != dim2) &&
                    error("bad dimensions for `", field,
                          "` field of OI-FITS extension ", extname)
                if channels == -1
                    channels = dim1
                else
                    channels == dim1 || error("incompatible number of ",
                                              "spectral channels for `",
                                              field,
                                              "` field of OI-FITS extension ",
                                              extname)
                end
            elseif spec.type != :STRING && dim1 != mult
                error("bad dimensions for `", field,
                      "` field of OI-FITS extension ", extname)
            end
        end

        # Store value.
        contents[field] = value
    end

    # Check that all mandatory fields have been given.
    nerrs = 0
    for field in def.fields
        spec = def.spec[field]
        if ! haskey(contents, field) && ! spec.isoptional
            warn("missing value for `", field, "` field of OI-FITS extension ",
                 extname)
            nerrs += 1
        end
    end
    if nerrs > 0
        error("some mandatory fields are missing in OI-FITS extension ",
              extname)
    end
    return get_datablock_type(extname)(contents)
end


###################################################
# MANAGING THE DATABLOCK CONTAINER (THE "MASTER") #
###################################################

new_master(datablocks::OIDataBlock...) = new_master(datablocks)
function new_master(datablocks::Union{AbstractVector{<:OIDataBlock},
                                      Tuple{Vararg{OIDataBlock}}})
    push!(OIMaster(), datablocks)
end

Base.push!(master::OIMaster) = master

Base.push!(master::OIMaster, datablocks::OIDataBlock...) =
    push!(master, datablocks)

function Base.push!(master::OIMaster,
                    data::Union{AbstractVector{<:OIDataBlock},
                                Tuple{Vararg{OIDataBlock}}})
    # First push each new data-block.  Second, update links of all stored
    # data-blocks.
    for i in eachindex(data)
        _push!(master, data[i])
    end
    _update_links!(master)
end

# Just push one data-block.
function _push!(master::OIMaster, db::OIDataBlock)
    is_attached(db) && error("data-block already attached")
    if isa(db, OITarget)
        getfield(master, :tgt) == nothing ||
            error("only one OI_TARGET data-block can be attached")
        setfield!(master, :tgt, db)
    elseif isa(db, OIWavelength)
        insname = fix_name(db.insname)
        ins = getfield(master, :ins)
        haskey(ins, insname) && error("master already have an ",
                                      "OI_WAVELENGTH data-block with ",
                                      "INSNAME=\"", insname, "\"")
        ins[insname] = db
    elseif isa(db, OIArray)
        arrname = fix_name(db.arrname)
        arr = getfield(master, :arr)
        haskey(arr, arrname) && error("master already have an ",
                                             "OI_ARRAY data-block with ",
                                             "ARRNAME=\"", arrname, "\"")
        arr[arrname] = db
    elseif isa(db, OICorrelation)
        corrname = fix_name(db.corrname)
        corr = getfield(master, :corr)
        haskey(corr, corrname) && error("master already have an ",
                                        "OI_CORR data-block with ",
                                        "CORRNAME=\"", corrname, "\"")
        corr[corrname] = db
    end
    push!(contents(master), db)
    setfield!(db, :attached, true)
end

function _update_links!(master::OIMaster)
    for i in eachindex(master)
        _update_links!(master, master[i])
    end
    master
end

_update_links!(master::OIMaster, db::OITarget) = nothing
_update_links!(master::OIMaster, db::OIArray) = nothing
_update_links!(master::OIMaster, db::OICorrelation) = nothing
_update_links!(master::OIMaster, db::OIWavelength) = nothing
_update_links!(master::OIMaster, db::OIPolarization) = _link_array!(master, db)
function _update_links!(master::OIMaster,
                        db::Union{OIVis,OIVis2,OIT3,OISpectrum})
    _link_array!(master, db)
    _link_instr!(master, db)
    _link_correl!(master, db)
end

function _link_array!(master::OIMaster, db::OIDataBlock)
    if db.arrname !== nothing
        setfield!(db, :arr, get(master.array, db.arrname, nothing))
    end
end

function _link_instr!(master::OIMaster, db::OIDataBlock)
    if db.insname !== nothing
        setfield!(db, :ins, get(master.instr, db.insname, nothing))
    end
end

function _link_correl!(master::OIMaster, db::OIDataBlock)
    if db.corrname !== nothing
        setfield!(db, :corr, get(master.correl, db.corrname, nothing))
    end
end

"""
    check_structure(obj)

check the consistency to the structure of the OI-FITS master or OI-FITS
data-block `obj`.

"""
function check_structure(master::OIMaster)
    master.target === nothing && error("missing mandatory OI_TARGET extension")
    for i in eachindex(master)
        check_structure(master[i])
    end
end

check_structure(db::OITarget) = nothing
check_structure(db::OIArray) = nothing
check_structure(db::OICorrelation) = nothing
check_structure(db::OIWavelength) = nothing

function check_structure(db::OIPolarization)
    _check_array(db)
    nothing
end

function check_structure(db::Union{OIVis,OIVis2,OIT3,OISpectrum})
    _check_array(db)
    _check_instr(db)
    _check_correl(db)
    nothing
end

function _check_array(db::OIDataBlock)
    db.arrname !== nothing && db.array === nothing &&
        warn("OI_ARRAY data-block with ARRNAME=\"",
             db.arrname, "\" not found in master")
end

function _check_instr(db::OIDataBlock)
    db.instr === nothing &&
        error("OI_WAVELENGTH data-block with INSNAME=\"",
              db.insname, "\" not found in master")
end

function _check_correl(db::OIDataBlock)
    db.corrname !== nothing && db.correl === nothing &&
        warn("OI_CORR data-block with CORRNAME=\"",
             db.corrname, "\" not found in master")
end

##############################
# METHODS FOR DATA SELECTION #
##############################

"""
    OIFITS.clone(db)

clones OI-FITS data-block `db`.  Most objects referenced by `db` will be shared
by the clone but neither the links (`db.array`, `db.instr`, etc.) nor the
fields not defined by OI-FITS standard.  This method is intended when a
data-block from a given `OIMaster` instance is to be pushed to another
`OIMaster` instance.

"""
function clone(db::OIDataBlock)
    (name, revn, defn) = Parser.get_description(db)
    dict = OIContents()
    for (key, val) in contents(db)
        if haskey(defn.spec, key)
            dict[key] = val
        end
    end
    return build_datablock(name, revn, dict)
end

"""
     OIFITS.select(obj, args)

yields an iterable over the data-blocks of an OI-FITS master object `obj` whose
extension names are specifed in `args`.

"""
function select(obj::OIMaster, args::AbstractString...)
    datablocks = Array{OIDataBlock}(undef, 0)
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
    idx = findfirst(x -> fix_name(x) == str, tgt.target)
    idx === nothing && return nothing
    select_target(obj, tgt.target_id[idx])
end

function _select_target(src::OIMaster, id::Int)
    dst = OIMaster()
    for db in src
        let sel = select_target(db, id)
            sel === nothing || _push!(dst, sel)
        end
    end
    _update_links!(dst)
end

_select_target(db::OIDataBlock, id::Int) = clone(db)

function _select_target(db::OITarget, id::Int)
    (name, revn, defn) = Parser.get_description(db)
    idx = findfirst(x -> x == id, db.target_id)
    idx == nothing && return nothing
    dict = OIContents()
    for (key, val) in contents(db)
        let spec = get(defn.spec, key, nothing)
            spec === nothing && continue
            if spec.iskeyword
                dict[key] = db[key]
            else
                dict[key] = db[key][idx:idx]
            end
        end
    end
    return build_datablock(name, revn, dict)
end

function _select_target(db::Union{OIVis,OIVis2,OIT3,
                                  OISpectrum,OIPolarization}, id::Int)
    (name, revn, defn) = Parser.get_description(db)
    target_id = db.target_id
    sel = findall(x -> x == id, target_id)
    length(sel) > 0 || return nothing
    dict = OIContents()
    cpy = (length(sel) == length(target_id)) # just copy?
    for (key, val) in contents(db)
        let spec = get(defn.spec, key, nothing)
            spec != nothing || continue
            if cpy || spec.iskeyword
                dict[key] = db[key]
            else
                # Copy a sub-array corresponding to the selection.  (The target
                # is specified for each last index.)
                src = db[key]
                @assert last(size(src)) == length(target_id)
                rank = ndims(src)
                dims = ntuple(i -> i == rank ? length(sel) : size(src, i), rank)
                dst = Array{eltype(src)}(undef, dims)
                dict[key] = dst
                if rank == 1
                    for j in 1:length(sel)
                        dst[j] = src[sel[j]]
                    end
                elseif rank == 2
                    for j in 1:length(sel)
                        dst[:,j] = src[:,sel[j]]
                    end
                else
                    error("unexpected rank ", rank)
                end
            end
        end
    end
    return build_datablock(name, revn, dict)
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

function select_wavelength(src::OIMaster, selector::Function)
    dst = OIMaster()
    for db in src
        sel = select_wavelength(db, selector)
        sel === nothing || _push!(dst, sel)
    end
    return _update_links!(dst)
end

select_wavelength(db::OIDataBlock, selector::Function) = clone(db)

function select_wavelength(db::Union{OIWavelength,OIVis,OIVis2,OIT3,
                                     OISpectrum,OIPolarization},
                           selector::Function)
    # Pre-select wavelengths.
    isintr = isa(db, OIWavelength) # source is instrument?
    wave = (isintr ? db.eff_wave : db.instr.eff_wave)
    wave === nothing && return nothing
    sel = findall(selector, wave)
    length(sel) > 0 || return nothing

    # Build a data-block with selected data.
    (name, revn, defn) = Parser.get_description(db)
    dict = OIContents()
    iscopy = (length(sel) == length(wave)) # just copy?
    for (key, val) in contents(db)
        spec = get(defn.spec, key, nothing)
        spec === nothing && continue
        if iscopy || spec.iskeyword || (spec.multiplier ≥ 0 && isintr == false)
            dict[key] = db[key]
        else
            # Copy a sub-array corresponding to the selection.  (The wavelength
            # corresponds to the first index.)
            let src = db[key]
                @assert first(size(src)) == length(wave)
                rank = ndims(src)
                dims = ntuple(i -> i == 1 ? length(sel) : size(src, i), rank)
                dst = Array{eltype(src)}(undef, dims)
                dict[key] = dst
                if rank == 1
                    for j in 1:length(sel)
                        dst[j] = src[sel[j]]
                    end
                elseif rank == 2
                    for j in 1:length(sel)
                        dst[j,:] = src[sel[j],:]
                    end
                else
                    error("unexpected rank ", rank)
                end
            end
        end
    end
    return build_datablock(name, revn, dict)
end
