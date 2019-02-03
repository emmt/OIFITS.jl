#
# oidata.jl --
#
# Implement OI-FITS data types.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015-2019: Éric Thiébaut.
#
#------------------------------------------------------------------------------

###############
# BASIC TYPES #
###############

# Conversion to real type (no-op. if already of the correct type).
to_real(x::Cdouble) = x
to_real(x::Array{Cdouble}) = x
to_real(x::Array{T}) where {T<:Real} = convert(Array{Cdouble}, x)
to_real(x::Real) = convert(Cdouble, x)

to_integer(x::Int) = x
to_integer(x::Array{Int}) = x
to_integer(x::Array{T}) where {T<:Integer} = convert(Array{Int}, x)
to_integer(x::Integer) = convert(Int, x)

to_complex(x::Complex{Cdouble}) = x
to_complex(x::Array{Complex{Cdouble}}) = x
to_complex(x::AbstractArray{T}) where {T<:Union{Real,Complex}} =
    convert(Array{Complex{Cdouble}}, x)
to_complex(x::Union{Real,Complex}) = convert(Complex{Cdouble}, x)

# OI-FITS files stores the following 4 different data types:
const _DTYPE_LOGICAL = 1 # for format letter 'L'
const _DTYPE_INTEGER = 2 # for format letters 'I' or 'J'
const _DTYPE_REAL    = 3 # for format letters 'D' or 'E'
const _DTYPE_COMPLEX = 4 # for format letter 'C'
const _DTYPE_STRING  = 5 # for format letter 'A'

# The following dictionary is used for quick conversion of FITS format
# letter to data type.
const _DATATYPES = Dict('l' =>  _DTYPE_LOGICAL,
                        'L' =>  _DTYPE_LOGICAL,
                        'i' =>  _DTYPE_INTEGER,
                        'I' =>  _DTYPE_INTEGER,
                        'j' =>  _DTYPE_INTEGER,
                        'J' =>  _DTYPE_INTEGER,
                        'e' =>  _DTYPE_REAL,
                        'E' =>  _DTYPE_REAL,
                        'd' =>  _DTYPE_REAL,
                        'D' =>  _DTYPE_REAL,
                        'c' =>  _DTYPE_COMPLEX,
                        'C' =>  _DTYPE_COMPLEX,
                        'a' =>  _DTYPE_STRING,
                        'A' =>  _DTYPE_STRING)

is_logical(::Any) = false
is_logical(::Bool) = true
is_logical(::Array{Bool}) = true

is_string(::Any) = false
is_string(::AbstractString) = true
is_string(::Array{T}) where {T<:AbstractString} = true

is_integer(::Any) = false
is_integer(::Integer) = true
is_integer(::Array{T}) where {T<:Integer} = true

is_real(::Any) = false
is_real(::Real) = true
is_real(::Array{T}) where {T<:Real} = true

is_complex(::Any) = false
is_complex(::Complex) = true
is_complex(::Array{T}) where {T<:Complex} = true


#################
# OI-FITS TYPES #
#################

const OIContents = Dict{Symbol,Any}

@compat abstract type OIDataBlock end

mutable struct OITarget <: OIDataBlock
    owner
    contents::OIContents
    OITarget(contents::OIContents) = new(nothing, contents)
end

mutable struct OIArray <: OIDataBlock
    owner
    contents::OIContents
    OIArray(contents::OIContents) = new(nothing, contents)
end

mutable struct OIWavelength <: OIDataBlock
    owner
    contents::OIContents
    OIWavelength(contents::OIContents) = new(nothing, contents)
end

mutable struct OICorrelation <: OIDataBlock
    owner
    contents::OIContents
    OICorrelation(contents::OIContents) = new(nothing, contents)
end

mutable struct OIPolarization <: OIDataBlock
    owner
    arr::Union{OIArray,Nothing}
    contents::OIContents
    OIPolarization(contents::OIContents) = new(nothing, nothing, contents)
end

mutable struct OIVis <: OIDataBlock
    owner
    arr::Union{OIArray,Nothing}
    ins::Union{OIWavelength,Nothing}
    corr::Union{OICorrelation,Nothing}
    contents::OIContents
    OIVis(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

mutable struct OIVis2 <: OIDataBlock
    owner
    arr::Union{OIArray,Nothing}
    ins::Union{OIWavelength,Nothing}
    corr::Union{OICorrelation,Nothing}
    contents::OIContents
    OIVis2(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

mutable struct OIT3 <: OIDataBlock
    owner
    arr::Union{OIArray,Nothing}
    ins::Union{OIWavelength,Nothing}
    corr::Union{OICorrelation,Nothing}
    contents::OIContents
    OIT3(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

mutable struct OISpectrum <: OIDataBlock
    owner
    arr::Union{OIArray,Nothing}
    ins::Union{OIWavelength,Nothing}
    corr::Union{OICorrelation,Nothing}
    contents::OIContents
    OISpectrum(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

# Correspondance between OI-FITS data-block names and Julia types.
const _DATABLOCKS = Dict("OI_TARGET"     => OITarget,
                         "OI_WAVELENGTH" => OIWavelength,
                         "OI_ARRAY"      => OIArray,
                         "OI_VIS"        => OIVis,
                         "OI_VIS2"       => OIVis2,
                         "OI_T3"         => OIT3,
                         "OI_SPECTRUM"   => OISpectrum,
                         "OI_CORR"       => OICorrelation,
                         "OI_INSPOL"     => OIPolarization)

const _EXTNAMES = Dict{valtype(_DATABLOCKS), keytype(_DATABLOCKS)}()

function __init__()
    # Dictionaries whose keys are types may not be precompiled (at least for
    # Julia versions anterior to 0.5).
    for (key, val) in _DATABLOCKS
        haskey(_EXTNAMES, val) && error("values must be unique")
        _EXTNAMES[val] = key
    end
end

get_dbname(db::OIDataBlock) = _EXTNAMES[typeof(db)]

# OIData is any OI-FITS data-block which contains interferometric data.
const OIData = Union{OIVis,OIVis2,OIT3}

# OIDataBlock can be indexed by the name (either as a string or as a
# symbol) of the field.
getindex(db::OIDataBlock, key::Symbol) = get(db.contents, key, nothing)
getindex(db::OIDataBlock, key::AbstractString) = getindex(db, Symbol(key))
haskey(db::OIDataBlock, key::Symbol) = haskey(db.contents, key)
haskey(db::OIDataBlock, key::AbstractString) = haskey(db.contents, Symbol(key))
keys(db::OIDataBlock) = keys(db.contents)

# OIMaster stores the contents of an OI-FITS file.  Data-blocks containing
# measurements (OI_VIS, OI_VIS2 and OI_T3) are stored into a vector and
# thus indexed by an integer.  Named data-blocks (OI_ARRAY and
# OI_WAVELENGTH) are indexed by their names (converted to upper case
# letters, with leading and trailing spaces stripped, multiple spaces
# replaced by a single ordinary space).
mutable struct OIMaster
    all::Vector{OIDataBlock}            # All data-blocks
    update_pending::Bool                # Update is needed?
    tgt::Union{OITarget,Nothing}
    arr::Dict{Name,OIArray}
    ins::Dict{Name,OIWavelength}
    corr::Dict{Name,OICorrelation}
    function OIMaster()
        new(Array{OIDataBlock}(undef, 0),
            false,
            nothing,
            Dict{Name,OIArray}(),
            Dict{Name,OIWavelength}(),
            Dict{Name,OICorrelation}())
    end
end

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
    print(io, "OI_MASTER: (", length(master.all), " data-block(s))")
    for db in master.all
        print(io, "\n    ")
        show(io, db)
    end
end


####################################
# PARSING OF DATABLOCK DEFINITIONS #
####################################

# OIFieldDef is used to store the definition of a keyword/column field.
mutable struct OIFieldDef
    name::Name      # Keyword/column name as a string.
    symb::Symbol    # Keyword/column symbolic name.
    keyword::Bool   # Is keyword? (otherwise column)
    optional::Bool  # Optional field?
    multiplier::Int # Multiplier: 1 for keywords, number of cells for columns
                    # (a negative number -N means an array of N dimensions each
                    # equal to the number of spectral channels.
    dtype::Int      # Data type.
    units::Name     # Units.
    descr::Name     # Description.
end

# OIDataBlockDef is used to store the definition of data-block.
mutable struct OIDataBlockDef
    dbname::Name
    fields::Vector{Symbol}        # ordered field symbolic names
    spec::Dict{Symbol,OIFieldDef} # dictionary of field specifications
    function OIDataBlockDef(dbname::AbstractString, vect::Vector{OIFieldDef})
        spec = Dict{Symbol,OIFieldDef}()
        fields = Array{Symbol}(undef, length(vect))
        for j in 1:length(vect)
            entry = vect[j]
            fields[j] = entry.symb
            spec[entry.symb] = entry
        end
        new(dbname, fields, spec)
    end
end

# OIFormatDef is used to store all the data-block definitions for a given
# revision number.
const OIFormatDef = Dict{Name,OIDataBlockDef}

# _FORMATS table is indexed by the datablock name and then by the revision
# number of the corresponding OI-FITS table.
const _FORMATS = Dict{Name,Vector{Union{OIDataBlockDef,Nothing}}}()

# The default format version number.
default_revision() = 2

# _FIELDS is a dictionary indexed by the data-block name (e,g., "OI_VIS"), each
# entry stores a set of its fields.
const _FIELDS = Dict{Name,Set{Symbol}}()

"""
    get_def(db, revn = default_revision())

yields the `OIDataBlockDef` for datablock of type `db` (e.g., "OI_TARGET") in
revison `revn` of OI-FITS standard.

"""
function get_def(dbname::AbstractString, revn::Integer = default_revision())
    if ! haskey(_FORMATS, dbname)
        error("unknown data-block: $db")
    end
    v = _FORMATS[dbname]
    if revn < 1 || revn > length(v) || v[revn] === nothing
        error("unsupported revision number: $revn")
    end
    return v[revn]
end

function get_def(db::OIDataBlock)
    get_def(get_dbname(db), get_revn(db))
end

function get_descr(db::OIDataBlock)
    name = get_dbname(db)
    revn = get_revn(db)
    return (name, revn, get_def(name, revn))
end

function add_def(dbname::AbstractString,
                 revn::Integer,
                 tbl::Vector{S}) where {S<:AbstractString}
    if (! startswith(dbname, "OI_") || dbname != uppercase(dbname)
        || occursin(" ", dbname))
        error("invalid data-block name: \"$dbname\"")
    end
    if revn < 1
        error("invalid revision number: $revn")
    end
    if haskey(_FORMATS, dbname) && revn <= length(_FORMATS[dbname]) &&
        _FORMATS[dbname][revn] != nothing
        error("data-block \"$dbname\" version $revn already defined")
    end

    fields = get(_FIELDS, dbname, Set{Symbol}())
    def = Array{OIFieldDef}(undef, 0)
    keyword = true
    for j in 1:length(tbl)
        row = strip(tbl[j])
        m = match(r"^([^ ]+) +([^ ]+) +(.*)$", row)
        if m === nothing
            if match(r"^-+$", row) === nothing
                error("syntax error in OI_FITS definition: \"$row\"")
            end
            keyword = false
            continue
        end
        name = uppercase(m.captures[1])
        symb = symbolicname(name)
        format = m.captures[2]
        descr = m.captures[3]
        optional = (format[1] == '?')
        i = (optional ? 2 : 1)
        dtype = get(_DATATYPES, format[i], nothing)
        if dtype === nothing
            error("invalid format type in OI_FITS definition: \"$row\"")
        end
        if keyword
            if length(format) != i
                error("invalid keyword format in OI_FITS definition: \"$row\"")
            end
            multiplier = 1
        else
            # Very naive code to parse the dimension list of the column format.
            if length(format) < i + 3 || format[i+1] != '(' ||
                format[end] != ')'
                error("missing column dimensions in OI_FITS definition: \"$row\"")
            end
            format = uppercase(format[i+2:end-1])
            if format == "W"
                multiplier = -1
            elseif format == "W,W"
                multiplier = -2
            else
                multiplier = parse(Int, format)
                if multiplier <= 0
                    error("invalid multiplier in OI_FITS definition: \"$row\"")
                end
            end
        end
        mp = match(r"^(.*[^ ]) +\[([^\]])\]$", descr)
        if mp === nothing
            units = ""
        else
            descr = mp.captures[1]
            units = mp.captures[2]
        end
        push!(def, OIFieldDef(name, symb, keyword, optional, multiplier,
                              dtype, units, descr))
        push!(fields, symb)
    end

    # Insert the data-block definition in the global table.
    if haskey(_FORMATS, dbname)
        v = _FORMATS[dbname]
    else
        v = Array{Union{OIDataBlockDef,Nothing}}(undef, 0)
        _FORMATS[dbname] = v
    end
    while revn > length(v)
        push!(v, nothing)
    end
    v[revn] = OIDataBlockDef(dbname, def)
    _FIELDS[dbname] = fields
end

# Convert the name of an OI-FITS keyword/column into a valid symbol.
function symbolicname(name::AbstractString)
    key = lowercase(name)
    if key == "oi_revn"
        return :revn
    else
        return Symbol(replace(key, r"[^a-z0-9_]" => '_'))
    end
end


#########################
# BUILDING OF DATABLOCK #
#########################

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
                       revn::Integer = default_revision(), kwds...)
            db = build_datablock($name, revn, kwds)
            if master != nothing
                attach!(master, db)
            end
            return db
        end
    end
end

function build_datablock(dbname::AbstractString, revn::Integer, kwds)
    def = get_def(dbname, revn)
    contents = OIContents()
    contents[:revn] = to_integer(revn)
    rows = -1         # number of rows in the OI-FITS table
    channels = -1     # number of spectral channels
    for (field, value) in kwds
        # Check whether this field exists.
        spec = get(def.spec, field, nothing)
        if spec === nothing
            error("data-block $dbname has no field \"$field\"")
        end

        # Check value type.
        if spec.dtype == _DTYPE_LOGICAL
            if ! is_logical(value)
                error("expecting boolean value for field \"$field\" in $dbname")
            end
        elseif spec.dtype == _DTYPE_INTEGER
            if ! is_integer(value)
                error("expecting integer value for field \"$field\" in $dbname")
            end
            value = to_integer(value)
        elseif spec.dtype == _DTYPE_REAL
            if ! is_real(value)
                error("expecting real value for field \"$field\" in $dbname")
            end
            value = to_real(value)
        elseif spec.dtype == _DTYPE_COMPLEX
            if ! is_complex(value)
                error("expecting complex value for field \"$field\" in $dbname")
            end
            value = to_complex(value)
        elseif spec.dtype == _DTYPE_STRING
            if ! is_string(value)
                error("expecting string value for field \"$field\" in $dbname")
            end
        else
            error("*** CORRUPTED WORKSPACE ***")
        end

        # Check value dimensions.
        dims = (isa(value, Array) ? size(value) : ())
        rank = length(dims)
        if spec.keyword
            if rank != 0
                error("expecting a scalar value for field \"$field\" in $dbname")
            end
        else
            # Check array rank.
            mult = spec.multiplier
            maxrank = (spec.dtype == _DTYPE_STRING || mult == 1 ? 1 :
                       mult == -2 ? 3 :
                       2)
            if rank > maxrank
                error("bad number of dimensions for field \"$field\" in $dbname")
            end

            # Fields may have up to 3 dimensions.  For now, the last dimension
            # is the number of rows.  FIXME: The ordering of dimensions must
            # be changed: the number of rows should be the frist dimension.
            dim0 = (rank >= 1 ? dims[end] : 1)
            dim1 = (rank >= 2 ? dims[1] : 1)
            dim2 = (rank >= 3 ? dims[2] : 1)

            if rows == -1
                rows = dim0
            elseif rows != dim0
                error("incompatible number of rows for field \"$field\" in $dbname")
            end

            if mult < 0
                # Expecting an N-by-W or N-by-W-by-W array (N is the number of
                # rows and W the number of channels).
                if mult == -2 && dim1 != dim2
                    error("bad dimensions for field \"$field\" in $dbname")
                end
                if channels == -1
                    channels = dim1
                elseif channels != dim1
                    error("incompatible number of spectral channels for field \"$field\" in $dbname")
                end
            elseif spec.dtype != _DTYPE_STRING && dim1 != mult
                error("bad dimensions for field \"$field\" in $dbname")
            end
        end

        # Store value.
        contents[field] = value
    end

    # Check that all mandatory fields have been given.
    nerrs = 0
    for field in def.fields
        spec = def.spec[field]
        if ! haskey(contents, field) && ! spec.optional
            @warn("missing value for field \"$field\" in $dbname")
            nerrs += 1
        end
    end
    if nerrs > 0
        error("some fields are missing in $dbname")
    end
    return _DATABLOCKS[dbname](contents)
end


###################################################
# MANAGING THE DATABLOCK CONTAINER (THE "MASTER") #
###################################################

# Letter case and trailing spaces are insignificant according to FITS
# conventions.
fixname(name::AbstractString) = uppercase(rstrip(name))

# OIDataBlock and OIMaster can be used as iterators.
@static if isdefined(Base, :iterate)
    import Base: iterate

    iterate(iter::OIDataBlock) = iterate(iter.contents)
    iterate(iter::OIDataBlock, state) = iterate(iter.contents, state)

    iterate(iter::OIMaster) = iterate(iter.all)
    iterate(iter::OIMaster, state) = iterate(iter.all, state)
else
    import Base: start, done, next

    start(iter::OIDataBlock) = start(iter.contents)
    done(iter::OIDataBlock, state) = done(iter.contents, state)
    next(iter::OIDataBlock, state) = next(iter.contents, state)

    start(iter::OIMaster) = start(update!(iter).all)
    done(iter::OIMaster, state) = done(iter.all, state)
    next(iter::OIMaster, state) = next(iter.all, state)
end

function select(master::OIMaster, args::AbstractString...)
    datablocks = Array{OIDataBlock}(undef, 0)
    for db in master
        if get_dbname(db) ∈ args
            push!(datablocks, db)
        end
    end
    return datablocks
end

"""
Assuming `mst` is an instance of `OIMaster`, then:
```
    get_target(mst)
```
yields the `OITarget` data-block of `mst` if any, `nothing` otherwise.

Assuming `tgt` is an instance of `OITarget`, then:
```
    get_target(tgt)
```
yields the "TARGET" column of `tgt` which is an array of target names.
"""
get_target(master::OIMaster) = update!(master).tgt


"""
Assuming `mst` is an instance of `OIMaster`, then:
```
    get_array(mst, arrname)
```
yields the `OIArray` data-block of `mst` whose name is `arrname` if any,
`nothing` otherwise.

Assuming `db` is an instance of a sub-type of `OIDataBlock`, then:
```
    get_array(db)
```
yields corresponding instance of `OIArray` if any, `nothing` otherwise.
"""
function get_array(master::OIMaster, arrname::AbstractString)
    get(update!(master).arr, fixname(arrname), nothing)
end
get_array(db::OIDataBlock) = nothing
get_array(db::Union{OIVis,OIVis2,OIT3}) = db.arr
get_array(db::OIArray) = db

"""
Assuming `mst` is an instance of `OIMaster`, then:
```
    get_instrument(mst, insname)
```
yields the `OIWavelength` data-block of `mst` whose name is `insname` if any,
`nothing` otherwise.

Assuming `db` is an instance of a sub-type of `OIDataBlock`, then:
```
    get_instrument(db)
```
yields corresponding instance of `OIWavelength` if any, `nothing` otherwise.
"""
function get_instrument(master::OIMaster, insname::AbstractString)
    get(update!(master).ins, fixname(insname), nothing)
end
get_instrument(db::OIDataBlock) = nothing
get_instrument(db::Union{OIVis,OIVis2,OIT3}) = db.ins
get_instrument(db::OIWavelength) = db

"""
Assuming `mst` is an instance of `OIMaster`, then:
```
    get_targets(mst)
```
yields the names of the targets defined in `mst`.
"""
function get_targets(master::OIMaster)
    tgt = get_target(master)
    return tgt === nothing ? Array{String}(undef, 0) : get_target(tgt)
end

"""
Assuming `mst` is an instance of `OIMaster`, then:
```
    get_arrays(mst)
```
yields the names of the interferometric arrays defined in `mst`.
"""
get_arrays(master::OIMaster) = collect(keys(update!(master).arr))

"""
Assuming `mst` is an instance of `OIMaster`, then:
```
    get_instruments(mst)
```
yields the names of the instruments defined in `mst`.
"""
get_instruments(master::OIMaster) = collect(keys(update!(master).ins))

function new_master(datablocks::OIDataBlock...)
    master = OIMaster()
    for db in datablocks
        attach!(master, db)
    end
    return update!(master)
end

function new_master(datablocks::Array{OIDataBlock})
    master = OIMaster()
    for db in datablocks
        attach!(master, db)
    end
    return update!(master)
end

function attach!(master::OIMaster, db::OIDataBlock)
    db.owner === nothing || error("data-block already attached")
    if isa(db, OITarget)
        if master.tgt != nothing
            error("only one OI_TARGET data-block can be attached")
        end
        master.tgt = db
    elseif isa(db, OIWavelength)
        insname = fixname(db[:insname])
        if haskey(master.ins, insname)
            error("master already have an OI_WAVELENGTH data-block with INSNAME=\"$insname\"")
        end
        master.ins[insname] = db
    elseif isa(db, OIArray)
        arrname = fixname(db[:arrname])
        if haskey(master.arr, arrname)
            error("master already have an OI_ARRAY data-block with ARRNAME=\"$arrname\"")
        end
        master.arr[arrname] = db
    elseif isa(db, OICorrelation)
        corrname = fixname(db[:corrname])
        if haskey(master.corr, corrname)
            error("master already have an OI_CORR data-block with CORRNAME=\"$corrname\"")
        end
        master.corr[corrname] = db
    end
    push!(master.all, db)
    master.update_pending = true
    db.owner = master
    nothing
end

function update!(master::OIMaster)
    if master.update_pending
        if master.tgt === nothing
            error("missing mandatory OI_TARGET data-block")
        end
        for db in master.all
            if haskey(db, :insname) && ! isa(db, Union{OIWavelength,OIPolarization})
                insname = fixname(db[:insname])
                db.ins = get(master.ins, insname, nothing)
                if db.ins === nothing
                    error("OI_WAVELENGTH data-block with INSNAME=\"$insname\" not found in master")
                end
            end
            if haskey(db, :arrname) && ! isa(db, OIArray)
                arrname = fixname(db[:arrname])
                db.arr = get(master.arr, arrname, nothing)
                if db.arr === nothing
                    @warn("OI_ARRAY data-block with ARRNAME=\"$arrname\" not found in master")
                end
            end
            if haskey(db, :corrname) && ! isa(db, OICorrelation)
                corrname = fixname(db[:corrname])
                db.corr = get(master.corr, corrname, nothing)
                if db.corr === nothing
                    @warn("OI_CORR data-block with CORRNAME=\"$corrname\" not found in master")
                end
            end
        end
        master.update_pending = false
    end
    return master
end

"""
### Clone a data-block

This method clones an existing data-block.  The array data are shared between
the clones but not the links (owner, target, instrument, etc.).  Only the
fields defined by OI-FITS standard are cloned.  This method is intended when
a data-block from a given `OIMaster` instance is to be attached to another
`OIMaster` instance.
"""
function clone(db::OIDataBlock)
    (name, revn, defn) = get_descr(db)
    data = Dict{Symbol,Any}()
    for (key, val) in db.contents
        if haskey(defn.spec, key)
            data[key] = val
        end
    end
    return build_datablock(name, revn, data)
end


"""
### Select data for a given target

The method `select_target` selects a subset of data corresponding to a given
target.   The general syntax is:
```
    out = select_target(inp, tgt)
```
where `inp` is the input data (can be an instance of `OIMaster` or of any
`OIDataBlock` sub-types), `tgt` is the target number or name.

The result `out` may share part of its contents with the input data `inp`.

The result may be `nothing` if the input contains no data for the given target.
"""
select_target(db::OIDataBlock, target::Integer) = clone(db)

function select_target(db::OITarget, target::Integer)
    (name, revn, defn) = get_descr(db)
    k = findfirst(id -> id == target, get_target_id(db))
    k == 0 && return
    data = Dict{Symbol,Any}()
    for (key, val) in db.contents
        spec = get(defn.spec, key, nothing)
        spec != nothing || continue
        if spec.keyword
            data[key] = db[key]
        else
            data[key] = db[key][k:k]
        end
    end
    return build_datablock(name, revn, data)
end

function select_target(db::OIData, target::Integer)
    (name, revn, defn) = get_descr(db)
    target_id = get_target_id(db)
    sel = find(id -> id == target, target_id)
    length(sel) > 0 || return
    data = Dict{Symbol,Any}()
    cpy = (length(sel) == length(target_id)) # just copy?
    for (key, val) in db.contents
        spec = get(defn.spec, key, nothing)
        spec != nothing || continue
        if cpy || spec.keyword
            data[key] = db[key]
        else
            # Copy a sub-array corresponding to the selection.  (The target is
            # specified for each last index.)
            src = db[key]
            @assert(size(src)[end] == length(target_id))
            rank = ndims(src)
            dims = ntuple(i -> i == rank ? length(sel) : size(src, i), rank)
            dst = Array{eltype(src)}(undef, dims)
            data[key] = dst
            if rank == 1
                for j in 1:length(sel)
                    dst[j] = src[sel[j]]
                end
            elseif rank == 2
                for j in 1:length(sel)
                    dst[:,j] = src[:,sel[j]]
                end
            else
                error("unexpected rank $rank")
            end
        end
    end
    return build_datablock(name, revn, data)
end

function select_target(master::OIMaster, target::Integer)
    result = OIMaster()
    for db in master
        sel = select_target(db, target)
        sel != nothing && attach!(result, sel)
    end
    return update!(result)
end

function select_target(master::OIMaster, target::AbstractString)
    db = get_target(master)
    if db != nothing
        k = findfirst(name -> name == target, get_target(db))
        k != 0 && return select_target(master, get_target_id(db)[k])
    end
    error("target name not found")
end

"""
### Select data at given wavelengths

The method `select_wavelength` selects a subset of data on the basis of their
wavelength.   The general syntax is:
```
    out = select_wavelength(inp, sel)
```
where `inp` is the input data (can be an instance of `OIMaster` or of any
`OIDataBlock` sub-types), `sel` is a predicate function which takes a
wavelength (in meters) as argument and returns true if this wavelength is to be
selected and false otherwise.

An alternative is:
```
    out = select_wavelength(inp, wmin, wmax)
```
to select wavelengths `w` such that `wmin ≤ w ≤ wmax`.  The wavelength
bounds are in meters.

The result `out` may share part of its contents with the input data `inp`.

The result may be `nothing` if the input contains no data at selected
wavelengths.
"""
function select_wavelength(inp::Union{OIMaster,OIDataBlock},
                           wavemin::Real, wavemax::Real)
    wmin = convert(Cdouble, wavemin)
    wmax = convert(Cdouble, wavemax)
    select_wavelength(inp, w -> wmin ≤ w ≤ wmax)
end

select_wavelength(db::OIDataBlock, selector::Function) = clone(db)

function select_wavelength(db::Union{OIWavelength,OIVis,OIVis2,OIT3,OISpectrum},
                           selector::Function)
    (name, revn, defn) = get_descr(db)
    wave = get_eff_wave(db)
    sel = find(selector, wave)
    length(sel) > 0 || return
    data = Dict{Symbol,Any}()
    cpy = (length(sel) == length(wave)) # just copy?
    iswave = (typeof(db) == OIWavelength)
    for (key, val) in db.contents
        spec = get(defn.spec, key, nothing)
        spec != nothing || continue
        if cpy || spec.keyword || (spec.multiplier >= 0 && ! iswave)
            data[key] = db[key]
        else
            # Copy a sub-array corresponding to the selection.  (The wavelength
            # corresponds to the first index.)
            src = db[key]
            @assert(size(src, 1) == length(wave))
            rank = ndims(src)
            dims = ntuple(i -> i == 1 ? length(sel) : size(src, i), rank)
            dst = Array{eltype(src)}(undef, dims)
            data[key] = dst
            if rank == 1
                for j in 1:length(sel)
                    dst[j] = src[sel[j]]
                end
            elseif rank == 2
                for j in 1:length(sel)
                    dst[j,:] = src[sel[j],:]
                end
            else
                error("unexpected rank $rank")
            end
        end
    end
    return build_datablock(name, revn, data)
end

function select_wavelength(master::OIMaster, selector::Function)
    result = OIMaster()
    for db in master
        sel = select_wavelength(db, selector)
        sel != nothing && attach!(result, sel)
    end
    return update!(result)
end
