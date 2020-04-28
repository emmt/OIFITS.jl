#
# builder.jl --
#
# Methods to parse OI-FITS extension definitions and to build data-block
# instances.
#
#------------------------------------------------------------------------------

module Builder

using Compat

using ..OIFITS
using ..OIFITS:
    OIData,
    contents,
    convert_key_type,
    fix_name,
    get_extname,
    is_attached,
    is_float,
    is_integer,
    is_logical,
    is_string,
    to_fieldname,
    to_float,
    to_integer,
    to_logical,
    to_string

"""

Type `OIFITS.Builder.FieldDefinition` is used to store the definition of a
keyword/column field.

"""
struct FieldDefinition
    name::String     # Keyword/column name as a string.
    symb::Symbol     # Keyword/column symbolic name.
    iskeyword::Bool  # Is keyword? (otherwise column)
    isoptional::Bool # Optional field?
    multiplier::Int  # Multiplier: 1 for keywords, number of cells for columns
                     # (a negative number -N means an array of N dimensions
                     # each equal to the number of spectral channels.
    type::Symbol     # Data type.
    units::String    # Units.
    descr::String    # Description.
end

"""

Type `OIFITS.Builder.DataBlockDefinition` is used to store the definition of
data-block.

"""
struct DataBlockDefinition
    extname::String                    # name of OI-FITS extension
    revn::Int                          # revision of OI-FITS extension
    fields::Vector{Symbol}             # ordered symbolic names of fields
    spec::Dict{Symbol,FieldDefinition} # dictionary of field specifications
    function DataBlockDefinition(extname::AbstractString, revn::Integer,
                                 defs::Vector{FieldDefinition})
        @assert revn ≥ 1
        spec = Dict{Symbol,FieldDefinition}()
        fields = Array{Symbol}(undef, length(defs))
        @inbounds for j in 1:length(defs)
            def = defs[j]
            fields[j] = def.symb
            spec[def.symb] = def
        end
        new(extname, revn, fields, spec)
    end
end

OIFITS.get_extname(def::DataBlockDefinition) = def.extname

# FIELDS is a dictionary indexed by the extension name (e,g., "OI_VIS"), each
# entry stores a set of its fields.
const FIELDS = Dict{String,Set{Symbol}}()

# REVISIONS is a dictionary indexed by the extension name (e,g., "OI_VIS") and
# whose value is the last revision number.
const REVISIONS = Dict{String,Int}()

# FORMATS table is indexed by a 2-tuple composed by the datablock name
# and by the revision number of the corresponding OI-FITS table.
const FORMATS = Dict{Tuple{String,Int},DataBlockDefinition}()
format_key(name::AbstractString, revn::Integer) =
    (to_string(name), to_integer(revn))

"""
    OIFITS.Builder.last_revision(arg) -> revn

yields the last revision number of OI-FITS extension given its name, a
data-block type or a data-block instance.  The value `0` is returned if `arg`
is not a known OI-FITS extension name.

"""
last_revision(db::OIDataBlock) = last_revision(db.extname)
last_revision(::Type{T}) where {T<:OIDataBlock} = last_revision(get_extname(T))
last_revision(extname::String) = get(REVISIONS, extname, 0)

"""
    OIFITS.Builder.initialize!(obj, dict) -> obj

set as many fields as possible to zero in structured object `obj`, then
populates the fields of the structured object `obj` with the key-value pairs in
`dict`, finally returns `obj`.  This method is intended for freshly created
objects whose fields have not yet been initialized.  If key `:revn` is
undefined in `dict`, the last revision of the format is assumed.

"""
function initialize!(obj::T, dict::AbstractDict) where {T<:OIDataBlock}
    clear_fields!(obj)
    if haskey(dict, :revn) == false
        obj.revn = last_revision(T)
    end
    length(dict) > 0 && populate!(obj, dict)
    return obj
end

"""
    OIFITS.Builder.populate!(obj, dict) -> obj

populates the fields of the structured object `obj` with the key-value pairs
in `dict` and returns `obj`.

"""
function populate!(obj::T, dict::AbstractDict) where {T}
    @assert isstructtype(T)
    for (key, val) in dict
        sym = Symbol(key)
        setproperty!(obj, sym, val)
    end
    return obj
end

"""
    OIFITS.Builder.clear_fields!(obj) -> obj

set as many fields as possible to zero in structured object `obj` and returns
`obj`.

"""
function clear_fields!(obj::T) where {T}
    @assert isstructtype(T)
    for sym in fieldnames(T)
        clear_field!(obj, sym)
    end
    return obj
end

"""
    OIFITS.Builder.clear_field!(obj, sym)

set field `sym` in structured object `obj` to zero if that makes sense, that is
if the fields is a number.

"""
clear_field!(obj::T, sym::Symbol) where {T} =
    (_clear_field!(obj, sym, fieldtype(T, sym)); nothing)

_clear_field!(obj, sym::Symbol, ::Type{T}) where {T<:Number} =
    setfield!(obj, sym, zero(T))

_clear_field!(obj, sym::Symbol, ::Type{<:Any}) =
    nothing

# Extend Base.push: to add data-blocks to
Base.push!(obj::OIMaster) = obj

Base.push!(obj::OIMaster, others::OIDataBlock...) = push!(obj, others)

function Base.push!(obj::OIMaster{T},
                    others::Union{AbstractVector{<:OIDataBlock},
                                Tuple{Vararg{OIDataBlock}}}) where {T}
    # First push each new data-block.  Second, update links of all stored
    # data-blocks.
    for i in eachindex(others)
        _push!(obj, convert_float_type(T, others[i]))
    end
    _update_links!(obj)
end

# Push nothing.
_push!(master::OIMaster, ::Nothing) = nothing

# Just push one data-block.
function _push!(master::OIMaster{T}, db::OIDataBlock{T}) where {T}
    # Do nothing if data-block already attached to master.
    is_attached(db, master) && return

    # Consider special cases.
    if isa(db, OITarget)
        isdefined(master, :target) &&
            error("only one OI_TARGET data-block can be attached")
        setfield!(master, :target, db)
    elseif isa(db, OIWavelength)
        insname = fix_name(db.insname)
        haskey(master.instr, insname) && error("master already have an ",
                                               "OI_WAVELENGTH data-block with ",
                                               "INSNAME=\"", insname, "\"")
        master.instr[insname] = db
    elseif isa(db, OIArray)
        arrname = fix_name(db.arrname)
        haskey(master.array, arrname) && error("master already have an ",
                                               "OI_ARRAY data-block with ",
                                               "ARRNAME=\"", arrname, "\"")
        master.array[arrname] = db
    elseif isa(db, OICorrelation)
        corrname = fix_name(db.corrname)
        haskey(master.correl, corrname) && error("master already have an ",
                                                 "OI_CORR data-block with ",
                                                 "CORRNAME=\"", corrname, "\"")
        master.correl[corrname] = db
    end

    # Push data-block to the list of data-blocks stored by the master.
    if is_attached(db)
        # Make a copy if data-block attached to any (other) master.
        push!(contents(master), copy(db))
    else
        push!(contents(master), db)
    end
    setfield!(db, :owner, master)
    nothing
end

function _update_links!(master::OIMaster)
    for i in eachindex(master)
        _update_links!(master[i])
    end
    master
end

_update_links!(db::OITarget) = nothing
_update_links!(db::OIArray) = nothing
_update_links!(db::OICorrelation) = nothing
_update_links!(db::OIWavelength) = nothing
_update_links!(db::OIPolarization) = _link_array!(db)
_update_links!(db::OIData) = begin
    _link_array!(db)
    _link_instr!(db)
    _link_correl!(db)
end

_link_array!(db::OIDataBlock) = _link_field!(db, :arrname, :array)
_link_instr!(db::OIDataBlock) = _link_field!(db, :insname, :instr)
_link_correl!(db::OIDataBlock) = _link_field!(db, :corrname, :correl)

function _link_field!(db::OIDataBlock, from::Symbol, to::Symbol)
    if (is_attached(db) && (name = getproperty(db, from)) !== nothing &&
        (lnk = get(getfield(db.owner, to), fix_name(name), nothing)) !== nothing)
        setfield!(db, to, lnk)
    end
    nothing
end

"""
     OIFITS.Builder.get_datatype(c)

yields the symbolic identifier of the data type in a FITS table column defined by
the letter `c`.

| Symbol     | Letters        |
|:-----------|:---------------|
| `:LOGICAL` | `'L'`          |
| `:INTEGER` | `'I'` or `'J'` |
| `:REAL`    | `'D'` or `'E'` |
| `:COMPLEX` | `'C'`          |
| `:STRING`  | `'A'`          |
| `:UNKNOWN` | other letters  |

"""
get_datatype(c::Char) =
    # It is about 2 times faster to search with a series of tests instead of
    # using a dictionary.  Converting to uppercase/lowercase is much longer so
    # we keep upper- and lowercase version of each letter.
    (c == 'A' ? :STRING  :
     c == 'C' ? :COMPLEX :
     c == 'D' ? :FLOAT   :
     c == 'E' ? :FLOAT   :
     c == 'I' ? :INTEGER :
     c == 'J' ? :INTEGER :
     c == 'L' ? :LOGICAL :
     c == 'a' ? :STRING  :
     c == 'c' ? :COMPLEX :
     c == 'd' ? :FLOAT   :
     c == 'e' ? :FLOAT   :
     c == 'i' ? :INTEGER :
     c == 'j' ? :INTEGER :
     c == 'l' ? :LOGICAL : :UNKNOWN)

"""
    OIFITS.Builder.get_description(db) -> (extname, revn, defn)

yields the FITS extension name, revision number and format definitions of the
OI-FITS data-block `db`.

"""
get_description(db::OIDataBlock) = (db.extname, db.revn, get_definition(db))

"""
    OIFITS.Builder.get_definition(db)

or

    OIFITS.Builder.get_definition(extname, revn)

yield the format definition for datablock `db` of for the OI-FITS extension
`extname` (e.g., "OI_TARGET") in revision `revn` of OI-FITS standard, throwing
an error if not found.  The returned value is an instance of `DataBlockDefinition`.

    OIFITS.get_definition(extname, revn, def)

is similar but returns `def` if the format is not found.

"""
function get_definition(extname::AbstractString, revn::Integer)
    val = get_definition(extname, revn, nothing)
    val === nothing && error("unknown OI-FITS extension \"", extname,
                             "\" (revision ", revn, ")")
    return val
end

get_definition(extname::AbstractString, revn::Integer, def) =
    get(FORMATS, format_key(extname, revn), def)

get_definition(db::OIDataBlock) = get_definition(db.extname, db.revn)

"""
    OIFITS.Builder.has_member(obj, sym)

yields whether object `obj` has a defined member of symbolic name `sym`.
Object `obj` can be a dictionary with symbolic keys or an `OIDataBlock`
instance.

"""
@inline has_member(dict::AbstractDict{Symbol}, sym::Symbol) = haskey(dict, sym)
@inline has_member(obj::T, sym::Symbol) where {T<:OIDataBlock} =
    hasfield(T, sym) && isdefined(obj, sym)

"""
    OIFITS.Builder.check_contents(obj)

checks whether object `obj` is a valid instance of an OI-FITS data-block.
This method checks that all mandatory members of `obj` are defined and have
consistent size.

    OIFITS.Builder.check_contents(obj, extname, revn)

checks whether the entries/fields of object `obj` are valid for an OI-FITS
extension `extname`, revison `revn`.  Object `obj` can be a dictionary with
symbolic keys or an `OIDataBlock` instance.

FIXME: check number of wavelengths etc. for a linked object.

"""
function check_contents(dict::AbstractDict{<:Union{Symbol,AbstractString}},
                        extname::AbstractString, revn::Integer)
    check_contents(convert_key_type(Symbol, dict),
                   to_string(extname), to_integer(revn))
end

check_contents(obj::OIDataBlock) = check_contents(obj, obj.extname, obj.revn)

function check_contents(obj::Union{AbstractDict{Symbol},OIDataBlock},
                        extname::String, revn::Int)
    def = get_definition(extname, revn)
    rows = -1         # number of rows in the OI-FITS table
    channels = -1     # number of spectral channels
    for (key, val) in obj
        # Ignore fields not in the specifications.
        spec = get(def.spec, key, nothing)
        spec === nothing && continue

        # Check value type.
        if spec.type === :LOGICAL
            is_logical(val) || error("expecting boolean value for `", key,
                                     "` field of OI-FITS extension ", extname)
        elseif spec.type === :INTEGER
            is_integer(val) || error("expecting integer value for `", key,
                                     "` field of OI-FITS extension ", extname)
        elseif spec.type === :FLOAT
            is_float(val) || error("expecting floating-point value for `", key,
                                   "` field of OI-FITS extension ", extname)
        elseif spec.type === :COMPLEX
            is_complex(val) || error("expecting complex value for `", key,
                                     "` field of OI-FITS extension ", extname)
        elseif spec.type === :STRING
            is_string(val) || error("expecting string value for `", key,
                                    "` field of OI-FITS extension ", extname)
        else
            error("unknown value type `:", spec.type, "`")
        end

        # Check value dimensions.
        dims = (isa(val, Array) ? size(val) : ())
        rank = length(dims)
        if spec.iskeyword
            rank == 0 || error("expecting a scalar value for field `", key,
                               "` in OI-FITS extension ", extname)
        else
            # Check array rank.
            mult = spec.multiplier
            maxrank = (spec.type === :STRING || mult == 1 ? 1 :
                       mult == -2 ? 3 : 2)
            rank > maxrank && error("bad number of dimensions for field `",
                                    key, "` of OI-FITS extension ", extname)

            # Fields may have up to 3 dimensions.  For now, the last dimension
            # is the number of rows.  FIXME: The ordering of dimensions must
            # be changed: the number of rows should be the first dimension.
            dim0 = (rank >= 1 ? dims[end] : 1)
            dim1 = (rank >= 2 ? dims[1] : 1)
            dim2 = (rank >= 3 ? dims[2] : 1)

            if rows == -1
                rows = dim0
            else
                rows == dim0 || error("incompatible number of rows for field `",
                                      key, "` of OI-FITS extension ",
                                      extname)
            end

            if mult < 0
                # Expecting an N-by-W or N-by-W-by-W array (N is the number of
                # rows and W the number of channels).
                (mult == -2 && dim1 != dim2) &&
                    error("bad dimensions for field `", key,
                          "` of OI-FITS extension ", extname)
                if channels == -1
                    channels = dim1
                else
                    channels == dim1 || error("incompatible number of spectral",
                                              "channels for field `", key,
                                              "` of OI-FITS extension ",
                                              extname)
                end
            elseif spec.type != :STRING && dim1 != mult
                error("bad dimensions for field `", key,
                      "` of OI-FITS extension ", extname)
            end
        end
    end

    # Check that all mandatory fields have been given.
    nerrs = 0
    for key in def.fields
        spec = def.spec[key]
        if (spec.isoptional || has_member(obj, key)) == false
            warn("missing value for field `", key, "` of OI-FITS extension ",
                 extname)
            nerrs += 1
        end
    end
    nerrs == 0 || error("some mandatory fields are missing in OI-FITS ",
                        "extension ", extname)
    nothing
end

"""
    OIFITS.Builder.define(extname, revn, defs)

defines revision `revn` of the OI-FITS extension `extname`.  The format of the
OI-FITS extension is described `defs`, a vector of strings like:

    ["KEYWORD FORMAT DESCR",
      ...,
      ...,
     "---------------------------",
     "COLUMN  FORMAT DESCR",
      ...,
      ...,
      ...]

where:

- `KEYWORD` is the keyword name in FITS header.

- `COLUMN` is the column name in FITS table (`TTYPE`).

- `FORMAT` may be prefixed with a `?` to indicate an optional field and is:

  - for keywords, a single letter indicating the type,

  - for columns, a letter indicating the type followed by list of dimensions in
    parenthesis (letter `W` for a dimension means number of wavelengths).

- `DESCR` is a short description; units, if any, are indicated at the end
  between square brackets.

There may be any number of keyword definitions and any number of column
definitions, the two parts are separated by a dash line like

    "--------------".

See also [`get_definition`](@ref).

"""
function define(extname::AbstractString,
                revn::Integer,
                tbl::Vector{<:AbstractString})
    # Check OI-FITS extension name and revision number.
    (startswith(extname, "OI_") && extname == uppercase(extname) &&
     ! occursin(" ", extname)) || error("invalid OI-FITS extension name: \"",
                                        extname, "\"")
    revn ≥ 1 || error("invalid OI-FITS revision number: ", revn)
    fmtkey = format_key(extname, revn)
    haskey(FORMATS, fmtkey) && error("revision ", revn,
                                      " of OI-FITS extension ", extname,
                                      " already defined")

    # Parse table of definitions.
    function bad_definition(reason::AbstractString, extname::AbstractString,
                            revn::Integer, linenum::Integer,
                            code::AbstractString)
        error(reason, " in definition of OI-FITS extension ", extname,
              " (revision ", revn, ", line ", linenum, "): \"", code, "\"")
    end
    fields = get(FIELDS, extname, Set{Symbol}())
    def = Array{FieldDefinition}(undef, 0)
    iskeyword = true
    for rownum in 1:length(tbl)
        row = strip(tbl[rownum])
        m = match(r"^([^ ]+) +([^ ]+) +(.*)$", row)
        if m === nothing
            match(r"^-+$", row) === nothing &&
                bad_definition("syntax error", extname, revn, rownum, row)
            iskeyword = false
            continue
        end
        name = uppercase(m.captures[1])
        symb = to_fieldname(name)
        format = m.captures[2]
        descr = m.captures[3]
        isoptional = (format[1] == '?')
        i = (isoptional ? 2 : 1)
        type = get_datatype(format[i])
        type == :UNKNOWN && bad_definition("invalid type letter",
                                           extname, revn, rownum, row)
        if iskeyword
            length(format) == i ||
                bad_definition("invalid keyword format",
                               extname, revn, rownum, row)
            multiplier = 1
        else
            # Very naive code to parse the dimension list of the column format.
            (length(format) > i + 2 && format[i+1] == '('
             && format[end] == ')') ||
                 bad_definition("missing column dimension(s)",
                                extname, revn, rownum, row)
            format = uppercase(format[i+2:end-1])
            if format == "W"
                multiplier = -1
            elseif format == "W,W"
                multiplier = -2
            else
                multiplier = tryparse(Int, format)
                (multiplier === nothing || multiplier < 1) &&
                    bad_definition("invalid multiplier",
                                   extname, revn, rownum, row)
            end
        end
        mp = match(r"^(.*[^ ]) +\[([^\]])\]$", descr)
        if mp === nothing
            units = ""
        else
            descr = mp.captures[1]
            units = mp.captures[2]
        end
        push!(def, FieldDefinition(name, symb, iskeyword, isoptional,
                                   multiplier, type, units, descr))
        push!(fields, symb)
    end

    # Insert the data-block definition in the global table.
    FORMATS[fmtkey] = DataBlockDefinition(extname, revn, def)
    FIELDS[extname] = fields
    REVISIONS[extname] = max(revn, get(REVISIONS, extname, 0))
end

include("formats.jl")

end # module
