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
# Copyright (C) 2015, Éric Thiébaut.
#
#------------------------------------------------------------------------------

###############
# BASIC TYPES #
###############

# The default floating point type corresponds to C double precision.
typealias OIReal Cdouble

# Conversion to real type (no-op. if already of the correct type).
to_real(x::OIReal) = x
to_real(x::Array{OIReal}) = x
to_real{T<:Real}(x::Array{T}) = convert(Array{OIReal}, x)
to_real(x::Real) = convert(OIReal, x)

# OI-FITS files stores the following 4 different data types:
const _DTYPE_LOGICAL = 1 # for format letter 'L'
const _DTYPE_INTEGER = 2 # for format letters 'I' or 'J'
const _DTYPE_REAL    = 3 # for format letters 'D' or 'E'
const _DTYPE_STRING  = 4 # for format letter 'A'

# The following dictionary is used for quick conversion of FITS format
# letter to data type.
const _DATATYPES = ['l' =>  _DTYPE_LOGICAL,
                    'L' =>  _DTYPE_LOGICAL,
                    'i' =>  _DTYPE_INTEGER,
                    'I' =>  _DTYPE_INTEGER,
                    'j' =>  _DTYPE_INTEGER,
                    'J' =>  _DTYPE_INTEGER,
                    'e' =>  _DTYPE_REAL,
                    'E' =>  _DTYPE_REAL,
                    'd' =>  _DTYPE_REAL,
                    'D' =>  _DTYPE_REAL,
                    'a' =>  _DTYPE_STRING,
                    'A' =>  _DTYPE_STRING]

is_logical(::Any) = false
is_logical(::Bool) = true
is_logical(::Array{Bool}) = true

is_string(::Any) = false
is_string(::String) = true
is_string{T<:String}(::Array{T}) = true

is_integer(::Any) = false
is_integer(::Integer) = true
is_integer{T<:Integer}(::Array{T}) = true

is_real(::Any) = false
is_real(::Real) = true
is_real{T<:Real}(::Array{T}) = true


#################
# OI-FITS TYPES #
#################

typealias OIContents Dict{Symbol,Any}

abstract OIDataBlock
type OITarget <: OIDataBlock
    owner
    contents::OIContents
    OITarget(contents::OIContents) = new(nothing, contents)
end

type OIArray <: OIDataBlock
    owner
    contents::OIContents
    OIArray(contents::OIContents) = new(nothing, contents)
end

type OIWavelength <: OIDataBlock
    owner
    contents::OIContents
    OIWavelength(contents::OIContents) = new(nothing, contents)
end

type OISpectrum <: OIDataBlock
    owner
    contents::OIContents
    OISpectrum(contents::OIContents) = new(nothing, contents)
end

type OIVis <: OIDataBlock
    owner
    arr::Union(OIArray,Nothing)
    ins::Union(OIWavelength,Nothing)
    contents::OIContents
    OIVis(contents::OIContents) = new(nothing, nothing, nothing, contents)
end

type OIVis2 <: OIDataBlock
    owner
    arr::Union(OIArray,Nothing)
    ins::Union(OIWavelength,Nothing)
    contents::OIContents
    OIVis2(contents::OIContents) = new(nothing, nothing, nothing, contents)
end

type OIT3 <: OIDataBlock
    owner
    arr::Union(OIArray,Nothing)
    ins::Union(OIWavelength,Nothing)
    contents::OIContents
    OIT3(contents::OIContents) = new(nothing, nothing, nothing, contents)
end

# Correspondance between OI-FITS data-block names and Julia types.
_DATABLOCKS = Dict{ASCIIString,DataType}(["OI_TARGET"     => OITarget,
                                          "OI_WAVELENGTH" => OIWavelength,
                                          "OI_ARRAY"      => OIArray,
                                          "OI_SPECTRUM"   => OISpectrum,
                                          "OI_VIS"        => OIVis,
                                          "OI_VIS2"       => OIVis2,
                                          "OI_T3"         => OIT3])
_EXTNAMES = Dict{DataType,ASCIIString}()

for (key, val) in _DATABLOCKS
    _EXTNAMES[val] = key
end

get_dbname(db::OIDataBlock) = _EXTNAMES[typeof(db)]

# OIData is any OI-FITS data-block which contains interferometric data.
typealias OIData Union(OIVis,OIVis2,OIT3)

# OIDataBlock can be indexed by the name (either as a string or as a
# symbol) of the field.
getindex(db::OIDataBlock, key::Symbol) = get(db.contents, key, nothing)
getindex(db::OIDataBlock, key::String) = getindex(db, symbol(key))
haskey(db::OIDataBlock, key::Symbol) = haskey(db.contents, key)
haskey(db::OIDataBlock, key::String) = haskey(db.contents, symbol(key))
keys(db::OIDataBlock) = keys(db.contents)

# OIDataBlock can be used as iterators.
start(db::OIDataBlock) = start(db.contents)
done(db::OIDataBlock, state) = done(db.contents, state)
next(db::OIDataBlock, state) = next(db.contents, state)

# OIMaster stores the contents of an OI-FITS file.  Data-blocks containing
# measurements (OI_VIS, OI_VIS2 and OI_T3) are stored into a vector and
# thus indexed by an integer.  Named data-blocks (OI_ARRAY and
# OI_WAVELENGTH) are indexed by their names (converted to upper case
# letters, with leading and trailing spaces stripped, multiple spaces
# replaced by a single ordinary space).
type OIMaster
    all::Vector{OIDataBlock}            # All data-blocks
    update_pending::Bool                # Update is needed?
    target::Union(OITarget,Nothing)
    arr::Dict{ASCIIString,OIArray}
    ins::Dict{ASCIIString,OIWavelength}
    function OIMaster()
        new(Array(OIDataBlock, 0),
            false,
            nothing,
            Dict{ASCIIString,OIArray}(),
            Dict{ASCIIString,OIWavelength}())
    end
end


####################################
# PARSING OF DATABLOCK DEFINITIONS #
####################################

# OIFieldDef is used to store the definition of a keyword/column field.
type OIFieldDef
    name::ASCIIString    # keyword/column name as a string
    symb::Symbol         # keyword/column symbolic name
    keyword::Bool        # is keyword? (otherwise column)
    multiplier::Int      # multiplier; for keywords, 0 means optional and 1
                         # means required; for columns, a negative number
                         # means abs(n) times the number of spectral
                         # channels;
    dtype::Int           # data type
    units::ASCIIString   # units
    descr::ASCIIString   # description
end

# OIDataBlockDef is used to store the definition of data-block.
type OIDataBlockDef
    dbname::ASCIIString
    fields::Vector{Symbol}        # ordered field symbolic names
    spec::Dict{Symbol,OIFieldDef} # dictionary of field specifications
    function OIDataBlockDef(dbname::ASCIIString, vect::Vector{OIFieldDef})
        spec = Dict{Symbol,OIFieldDef}()
        fields = Array(Symbol, length(vect))
        for j in 1:length(vect)
            entry = vect[j]
            fields[j] = entry.symb
            spec[entry.symb] = entry
        end
        new(dbname, fields, spec)
    end
end

# OIFormatDef is used to store all the data-block definitions
# for a given revision number.
typealias OIFormatDef Dict{ASCIIString,OIDataBlockDef}

# _FORMATS array is indexed by the revision number.
_FORMATS = Array(OIFormatDef, 0)

# The default format version number is the highest registered one.
default_revision() = length(_FORMATS)

# _FIELDS is a dictionary indexed by the data-block name (e,g., "OI_VIS"), each
# entry stores a set of its fields.
_FIELDS = Dict{ASCIIString,Set{Symbol}}()

# get_def(db, revn) -- yields the OIDataBlockDef for datablock of type `db`
#                      (e.g., "OI_TARGET") in revison `revn` of OI-FITS
#                      standard.
#
function get_def(dbname::ASCIIString, revn::Integer)
    if revn < 0 || revn > length(_FORMATS)
        error("unsupported revision number: $revn")
    end
    if ! haskey(_FORMATS[revn], dbname)
        error("unknown data-block: $db")
    end
    _FORMATS[revn][dbname]
end

function add_def(dbname::ASCIIString, revn::Integer, tbl::Vector{ASCIIString})
    if ! beginswith(dbname, "OI_") || dbname != uppercase(dbname) || contains(dbname, " ")
        error("invalid data-block name: \"$db\"")
    end
    if revn < 0 || revn > 2
        error("invalid revision number: $revn")
    end

    fields = get(_FIELDS, dbname, Set{Symbol}())
    def = Array(OIFieldDef, 0)
    keyword = true
    for j in 1:length(tbl)
        row = strip(tbl[j])
        m = match(r"^([^ ]+) +([^ ]+) +([^ ]+) +(.*)$", row)
        if m == nothing
            if match(r"^-+$", row) == nothing
                error("syntax error in OI_FITS definition: \"$row\"")
            end
            keyword = false
            continue
        end
        name = uppercase(m.captures[1])
        symb = name2symbol(name)
        dtype = get(_DATATYPES, m.captures[2][end], nothing)
        if dtype == nothing
            error("invalid format type in OI_FITS definition: \"$row\"")
        end
        multiplier = int(m.captures[2][1:end-1])
        units = m.captures[3]
        descr = m.captures[4]
        if keyword
            if ! (multiplier == 0 || multiplier == 1)
                error("multiplier must be 0 or 1 in OI_FITS definition: \"$row\"")
            end
        else
            if ! (multiplier == -1 || multiplier > 0)
                error("invalid multiplier in OI_FITS definition: \"$row\"")
            end
        end
        push!(def, OIFieldDef(name, symb, keyword, multiplier,
                              dtype, units, descr))
        push!(fields, symb)
    end

    # Insert the data-block definition in the global table.
    while revn > length(_FORMATS)
        push!(_FORMATS, OIFormatDef())
    end
    _FORMATS[revn][dbname] = OIDataBlockDef(dbname, def)
    _FIELDS[dbname] = fields
end

# This version takes care of converting the string into a symbol.
get_def(db::String, revn::Integer) =  get_def(symbol(db), revn)

# get_def(db) -- returns the definition for the default format version.
get_def(db) =  get_def(db, default_revision())

# Convert the name of an OI-FITS keyword/column into a valid symbol.
function name2symbol(name::String)
    key = lowercase(name)
    if key == "oi_revn"
        return :revn
    else
        return symbol(replace(key, r"[^a-z0-9_]", '_'))
    end
end


#########################
# BUILDING OF DATABLOCK #
#########################

# Build constructors for the various data-blocks types.
for (func, name) in ((:new_target,     "OI_TARGET"),
                     (:new_array,      "OI_ARRAY"),
                     (:new_wavelength, "OI_WAVELENGTH"),
                     (:new_spectrum,   "OI_SPECTRUM"),
                     (:new_vis,        "OI_VIS"),
                     (:new_vis2,       "OI_VIS2"),
                     (:new_t3,         "OI_T3"))
    @eval begin
        function $func(master::Union(OIMaster, Nothing)=nothing;
                       revn::Integer=default_revision(), args...)
            db = build_datablock($name, revn, args)
            if master != nothing
                attach!(master, db)
            end
            return db
        end
    end
end

function build_datablock(dbname::ASCIIString, revn::Integer, args)
    def = get_def(dbname, revn)
    contents = OIContents([:revn => revn])
    nrows = -1     # number of measurements
    ncols = -1     # number of spectral channels
    nerrs = 0      # number of errors so far
    for (field, value) in args
        # Check whether this field exists.
        spec = get(def.spec, field, nothing)
        if spec == nothing
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
            value = int(value)
        elseif spec.dtype == _DTYPE_REAL
            if ! is_real(value)
                error("expecting real value for field \"$field\" in $dbname")
            end
            value = to_real(value)
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
            n = (spec.multiplier == 1 || spec.dtype == _DTYPE_STRING ? 1 : 2)
            if rank != n
                error("expecting a $n-D array value for field \"$field\" in $dbname")
            end
            if nrows == -1
                nrows = dims[end]
            elseif nrows != dims[end]
                error("incompatible number of rows for value of field \"$field\" in $dbname")
            end
            if spec.multiplier < 0
                if ncols == -1
                    ncols = dims[2]
                elseif ncols != dims[2]
                    error("incompatible number of columns for value of field \"$field\" in $dbname")
                end
            end
        end

        # Store value.
        contents[field] = value
    end

    # Check that all mandatory fields have been given.
    nerrs = 0
    for field in def.fields
        spec = def.spec[field]
        if ! haskey(contents, field) && spec.multiplier != 0
            warn("missing value for field \"$field\" in $dbname")
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
fixname(name::String) = uppercase(rstrip(name))

# Make OIMaster an iterator.
start(master::OIMaster) = start(update(master).all)
done(master::OIMaster, state) = done(master.all, state)
next(master::OIMaster, state) = next(master.all, state)

function select(master::OIMaster, args::String...)
    datablocks = Array(OIDataBlock, 0)
    for db in master
        if get_dbname(db) ∈ args
            push!(datablocks, db)
        end
    end
    return datablocks
end

oifits_target(master::OIMaster) = master.target
oifits_array(master::OIMaster, arrname::String) = master.array[fixname(arrname)]
oifits_wavelength(master::OIMaster, insname::String) = master.wavelength[fixname(insname)]
#oifits_vis(master::OIMaster) = master.vis
#oifits_vis(master::OIMaster, k::Integer) = master.vis[k]
#oifits_vis2(master::OIMaster) = master.vis2
#oifits_vis2(master::OIMaster, k::Integer) = master.vis2[k]
#oifits_t3(master::OIMaster) = master.t3
#oifits_t3(master::OIMaster, k::Integer) = master.t3[k]

function new_master(datablocks::OIDataBlock...)
    master = OIMaster()
    for db in datablocks
        attach!(master, db)
    end
    master
end

function new_master(datablocks::Array{OIDataBlock})
    master = OIMaster()
    for db in datablocks
        attach!(master, db)
    end
    master
end

function attach!(master::OIMaster, db::OIDataBlock)
    db.owner == nothing || error("data-block already attached")
    if isa(db, OITarget)
        if master.target != nothing
            error("only one OI_TARGET data-block can be attached")
        end
        master.target = db
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
    end
    push!(master.all, db)
    master.update_pending = true
    db.owner = master
    nothing
end

function update(master::OIMaster)
    if master.update_pending
        if master.target == nothing
            error("missing mandatory OI_TARGET data-block")
        end
        for db in master.all
            if isa(db, OIData)
                insname = fixname(db[:insname])
                arrname = fixname(db[:arrname])
                db.ins = get(master.ins, insname, nothing)
                db.arr = get(master.arr, arrname, nothing)
                if db.ins == nothing
                    error("OI_WAVELENGTH data-block with INSNAME=\"$insname\" not found in master")
                end
                if db.arr == nothing
                    warn("OI_ARRAY data-block with ARRNAME=\"$arrname\" not found in master")
                end
            end
        end
        master.update_pending = false
    end
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
