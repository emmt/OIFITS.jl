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
# Copyright (C) 2015: Éric Thiébaut.
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
const _DTYPE_COMPLEX = 4 # for format letter 'C'
const _DTYPE_STRING  = 5 # for format letter 'A'

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
                    'c' =>  _DTYPE_COMPLEX,
                    'C' =>  _DTYPE_COMPLEX,
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

type OICorrelation <: OIDataBlock
    owner
    contents::OIContents
    OICorrelation(contents::OIContents) = new(nothing, contents)
end

type OIPolarization <: OIDataBlock
    owner
    contents::OIContents
    OIPolarization(contents::OIContents) = new(nothing, contents)
end

type OIVis <: OIDataBlock
    owner
    arr::Union(OIArray,Nothing)
    ins::Union(OIWavelength,Nothing)
    corr::Union(OICorrelation,Nothing)
    contents::OIContents
    OIVis(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

type OIVis2 <: OIDataBlock
    owner
    arr::Union(OIArray,Nothing)
    ins::Union(OIWavelength,Nothing)
    corr::Union(OICorrelation,Nothing)
    contents::OIContents
    OIVis2(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

type OIT3 <: OIDataBlock
    owner
    arr::Union(OIArray,Nothing)
    ins::Union(OIWavelength,Nothing)
    corr::Union(OICorrelation,Nothing)
    contents::OIContents
    OIT3(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

type OISpectrum <: OIDataBlock
    owner
    arr::Union(OIArray,Nothing)
    ins::Union(OIWavelength,Nothing)
    corr::Union(OICorrelation,Nothing)
    contents::OIContents
    OISpectrum(contents::OIContents) = new(nothing, nothing, nothing, nothing, contents)
end

# Correspondance between OI-FITS data-block names and Julia types.
_DATABLOCKS = Dict{ASCIIString,DataType}(["OI_TARGET"     => OITarget,
                                          "OI_WAVELENGTH" => OIWavelength,
                                          "OI_ARRAY"      => OIArray,
                                          "OI_VIS"        => OIVis,
                                          "OI_VIS2"       => OIVis2,
                                          "OI_T3"         => OIT3,
                                          "OI_SPECTRUM"   => OISpectrum,
                                          "OI_CORR"       => OICorrelation,
                                          "OI_INSPOL"     => OIPolarization])
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
    corr::Dict{ASCIIString,OICorrelation}
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
    optional::Bool       # optional field?
    multiplier::Int      # multiplier: 1 for keywords, number of cells for
                         # columns (a negative number -N means an array of
                         # N dimensions each equal to the number of spectral
                         # channels
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

# _FORMATS table is indexed by the datablock name and then
# by the revision number of the corresponding OI-FITS table.
_FORMATS = Dict{ASCIIString,Vector{Union(OIDataBlockDef,Nothing)}}()

# The default format version number.
default_revision() = 2

# _FIELDS is a dictionary indexed by the data-block name (e,g., "OI_VIS"), each
# entry stores a set of its fields.
_FIELDS = Dict{ASCIIString,Set{Symbol}}()

# get_def(dbname, revn) -- Yields the definition of datablock of type `dbname`
#               (e.g., "OI_TARGET") of OI-FITS standard.  Argument `revn` is
#               the revision number of the definition.  The type of the result
#               is `OIDataBlockDef`.
function get_def(dbname::ASCIIString, revn::Integer)
    if ! haskey(_FORMATS, dbname)
        error("unknown data-block: $db")
    end
    v = _FORMATS[dbname]
    if revn < 1 || revn > length(v) || v[revn] == nothing
        error("unsupported revision number: $revn")
    end
    v[revn]
end

function add_def(dbname::ASCIIString, revn::Integer, tbl::Vector{ASCIIString})
    if ! beginswith(dbname, "OI_") || dbname != uppercase(dbname) || contains(dbname, " ")
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
    def = Array(OIFieldDef, 0)
    keyword = true
    for j in 1:length(tbl)
        row = strip(tbl[j])
        m = match(r"^([^ ]+) +([^ ]+) +(.*)$", row)
        if m == nothing
            if match(r"^-+$", row) == nothing
                error("syntax error in OI_FITS definition: \"$row\"")
            end
            keyword = false
            continue
        end
        name = uppercase(m.captures[1])
        symb = name2symbol(name)
        format = m.captures[2]
        descr = m.captures[3]
        optional = (format[1] == '?')
        i = (optional ? 2 : 1)
        dtype = get(_DATATYPES, format[i], nothing)
        if dtype == nothing
            error("invalid format type in OI_FITS definition: \"$row\"")
        end
        if keyword
            if length(format) != i
                error("invalid keyword format in OI_FITS definition: \"$row\"")
            end
            multiplier = 1
        else
            # Very naive code to parse the dimension list of the column
            # format.
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
                multiplier = int(format)
                if multiplier <= 0
                    error("invalid multiplier in OI_FITS definition: \"$row\"")
                end
            end
        end
        m = match(r"^(.*[^ ]) +\[([^\]])\]$", descr)
        if m == nothing
            units = ""
        else
            descr = m.captures[1]
            units = m.captures[2]
        end
        push!(def, OIFieldDef(name, symb, keyword, optional, multiplier,
                              dtype, units, descr))
        push!(fields, symb)
    end

    # Insert the data-block definition in the global table.
    if haskey(_FORMATS, dbname)
        v = _FORMATS[dbname]
    else
        v = Array(Union(OIDataBlockDef,Nothing), 0)
        _FORMATS[dbname] = v
    end
    while revn > length(v)
        push!(v, nothing)
    end
    v[revn] = OIDataBlockDef(dbname, def)
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
                     (:new_vis,        "OI_VIS"),
                     (:new_vis2,       "OI_VIS2"),
                     (:new_t3,         "OI_T3"),
                     (:new_spectrum,   "OI_SPECTRUM"),
                     (:new_corr,       "OI_CORR"),
                     (:new_inspol,     "OI_INSPOL"))
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
    rows = -1         # number of rows in the OI-FITS table
    channels = -1     # number of spectral channels
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
    nerrs = 0      # number of errors so far
    for field in def.fields
        spec = def.spec[field]
        if ! haskey(contents, field) && ! spec.optional
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

function update(master::OIMaster)
    if master.update_pending
        if master.target == nothing
            error("missing mandatory OI_TARGET data-block")
        end
        for db in master.all
            if haskey(db, :insname) && ! isa(db, OIWavelength)
                insname = fixname(db[:insname])
                db.ins = get(master.ins, insname, nothing)
                if db.ins == nothing
                    error("OI_WAVELENGTH data-block with INSNAME=\"$insname\" not found in master")
                end
            end
            if haskey(db, :arrname) && ! isa(db, OIArray)
                arrname = fixname(db[:arrname])
                db.arr = get(master.arr, arrname, nothing)
                if db.arr == nothing
                    warn("OI_ARRAY data-block with ARRNAME=\"$arrname\" not found in master")
                end
            end
            if haskey(db, :corrname) && ! isa(db, OICorrelation)
                corrname = fixname(db[:corrname])
                db.corr = get(master.corr, corrname, nothing)
                if db.corr == nothing
                    warn("OI_CORR data-block with CORRNAME=\"$corrname\" not found in master")
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
