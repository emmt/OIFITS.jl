#
# builder.jl --
#
# Methods to parse OI-FITS extension definitions and to build data-block
# instances.
#
#------------------------------------------------------------------------------

module Builder

#using Compat

using ..OIFITS
using ..OIFITS:
    OIData,
    contents,
    convert_float_type,
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
    name::String    # Keyword/column name as a string.
    symb::Symbol    # Keyword/column symbolic name.
    type::Symbol    # Data type.
    multiplier::Int # Multiplier: 0 for keywords, otherwise number of cells for
                    # columns except that a negative number -N means an array
                    # of N dimensions each equal to the number of spectral
                    # channels.
    optional::Bool  # Optional field?
    descr::String   # Description.
    units::String   # Units.
end

"""

Type `OIFITS.Builder.DataBlockDefinition` is used to store the definition of
a given data-block type.

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

"""

`OIFITS.Builder.FIELDS` is a dictionary indexed by the extension name (e,g.,
"OI_VIS"), each entry stores a set of its fields.

"""
const FIELDS = Dict{String,Set{Symbol}}()

"""

`OIFITS.Builder.REVISIONS` is a dictionary indexed by the extension name (e,g.,
"OI_VIS") and whose value is the last revision number.

"""
const REVISIONS = Dict{String,Int}()

"""

`OIFITS.Builder.FORMATS` is a dictionary indexed by a 2-tuple composed by the
datablock name and by the revision number of the corresponding OI-FITS table.

"""
const FORMATS = Dict{Tuple{String,Int},DataBlockDefinition}()
format_key(name::AbstractString, revn::Integer) =
    (to_string(name), to_integer(revn))

"""
    OIFITS.Builder.last_revision(arg) -> revn

yields the last revision number of OI-FITS extension given its name, a
data-block type or a data-block instance.  The value `0` is returned if `arg`
is not a known OI-FITS extension name.

"""
last_revision(db::OIDataBlock) = last_revision(typeof(db))
last_revision(T::Type{<:OIDataBlock}) = last_revision(get_extname(T))
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
    _clear_field!(obj, sym, fieldtype(T, sym))

_clear_field!(obj, sym::Symbol, ::Type{T}) where {T<:Number} = begin
    setfield!(obj, sym, zero(T))
    nothing
end

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


common_length(x::Int) = x
common_length(x::Int, y::Int) = (x === y ? x : -1)
@inline common_length(x::Int, y::Int, vals::Int...) =
    common_length(common_length(x, y), vals)


checked_length(db::OITarget) = begin
    len = common_length(length(db.target_id),
                        length(db.target),
                        length(db.raep0),
                        length(db.decep0),
                        length(db.equinox),
                        length(db.ra_err),
                        length(db.dec_err),
                        length(db.sysvel),
                        length(db.veltyp),
                        length(db.veldef),
                        length(db.pmra),
                        length(db.pmdec),
                        length(db.pmra_err),
                        length(db.pmdec_err),
                        length(db.parallax),
                        length(db.para_err),
                        length(db.spectyp),
                        length(db.category))
    len ≥ 0 ||
        throw(DimensionMismatch("all members of OITarget must have the same length"))
    return len
end

# Just push one data-block.
function push!(master::OIMaster{T}, db::OITarget{T}) where {T}

    len = checked_length(master.target)
    append!(master.target.target_id, db.target_id)
    append!(master.target.target, db.target)
    append!(master.target.raep0, db.raep0)
    append!(master.target.decep0, db.decep0)
    append!(master.target.equinox, db.equinox)
    append!(master.target.ra_err, db.ra_err)
    append!(master.target.dec_err, db.dec_err)
    append!(master.target.sysvel, db.sysvel)
    append!(master.target.veltyp, db.veltyp)
    append!(master.target.veldef, db.veldef)
    append!(master.target.pmra, db.pmra)
    append!(master.target.pmdec, db.pmdec)
    append!(master.target.pmra_err, db.pmra_err)
    append!(master.target.pmdec_err, db.pmdec_err)
    append!(master.target.parallax, db.parallax)
    append!(master.target.para_err, db.para_err)
    append!(master.target.spectyp, db.spectyp)
    append!(master.target.category, db.category)
    master.checked = false
    return master
end

function push!(master::OIMaster{T}, db::OIWavelength{T}) where {T}
    insname = fix_name(db.insname) # FIXME: should be done as soon as possible
    haskey(master.insname_dict, insname) &&
        error("master already have an OI_WAVELENGTH data-block with INSNAME=\"",
              insname, "\"")
    push!(master.wavelength, db)
    master.insname_dict[insname] = detached_copy(db, insname=insname)
    master.checked = false
    return master
end

function push!(master::OIMaster{T}, db::OIArray{T}) where {T}
    arrname = fix_name(db.arrname) # FIXME: should be done as soon as possible
    haskey(master.arrname_dict, arrname) &&
        error("master already have an OI_ARRAY data-block with ARRNAME=\"",
              arrname, "\"")
    master.arrname_dict[arrname] = detached_copy(db, arrname=arrname)
    master.checked = false
    return master
end

function push!(master::OIMaster{T}, db::OICorr{T}) where {T}
    corrname = fix_name(db.corrname) # FIXME: should be done as soon as possible
    haskey(master.corrname_dict, corrname) &&
        error("master already have an OI_CORR data-block with CORRNAME=\"",
              corrname, "\"")
    master.corrname_dict[corrname] = detached_copy(db, corrname=corrname)
    master.checked = false
    return master
end

for (sym, DB) in ((:vis,    OIVis),
                  (:vis2,   OIVis2),
                  (:t3,     OIT3),
                  (:flux,   OIFlux),
                  (:inspol, OIInsPol))
    @eval function push!(master::OIMaster{T}, db::$DB{T}) where {T}
        push!(master.$sym, detached_copy(db))
        master.checked = false
        return master
    end
end

function _update_links!(master::OIMaster)
    for i in eachindex(master)
        _update_links!(master[i])
    end
    master
end

_update_links!(db::OITarget, master::OIMaster) = nothing
_update_links!(db::OIArray, master::OIMaster) = nothing
_update_links!(db::OICorr, master::OIMaster) = nothing
_update_links!(db::OIWavelength, master::OIMaster) = nothing
_update_links!(db::OIInsPol, master::OIMaster) = begin
    _link_array!(db, master)
    _link_instr!(db, master)
    nothing
end
_update_links!(db::OIDataBlock, master::OIMaster) = begin
    _link_array!(db, master)
    _link_instr!(db, master)
    _link_correl!(db, master)
    nothing
end

_link_array!(db::OIDataBlock, master::OIMaster) =
    db.array = master.array_dict[db.arrname]

_link_instr!(db::OIDataBlock, master::OIMaster) =
    db.instr = master.instr_dict[db.insname]

_link_correl!(db::OIDataBlock, master::OIMaster) =
    db.correl = master.correl_dict[db.corrname]


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
    # we keep the upper- and lower-case versions of each letter.
    (((c === 'A')||(c === 'a')) ? :STRING  :
     ((c === 'C')||(c === 'c')) ? :COMPLEX :
     ((c === 'D')||(c === 'd')||
      (c === 'E')||(c === 'e')) ? :FLOAT   :
     ((c === 'I')||(c === 'i')||
      (c === 'J')||(c === 'j')) ? :INTEGER :
     ((c === 'L')||(c === 'l')) ? :LOGICAL : :UNKNOWN)

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


is_keyword(def::FieldDefinition) = (def.multiplier == 0)
is_optional(def::FieldDefinition) = def.optional

function check_field_value(nrows::Ref{Int}, nwaves::Ref{Int},
                           val, spec::FieldDefinition, extname)
    # Check value type.
    if spec.type === :L
        is_logical(val) || error("expecting boolean value for `", spec.name,
                                 "` field of OI-FITS extension ", extname)
    elseif spec.type ∈ (:I, :J)
        is_integer(val) || error("expecting integer value for `", spec.name,
                                 "` field of OI-FITS extension ", extname)
    elseif spec.typ ∈ (:D, :E)
        is_float(val) || error("expecting floating-point value for `", spec.name,
                               "` field of OI-FITS extension ", extname)
    elseif spec.type === :C
        is_complex(val) || error("expecting complex value for `", spec.name,
                                 "` field of OI-FITS extension ", extname)
    elseif spec.type === :A
        is_string(val) || error("expecting string value for `", spec.name,
                                "` field of OI-FITS extension ", extname)
    else
        error("unknown value type `:", spec.type, "`")
    end
end

function check_contents(db::OIDataBlock)
    def = get_definition(db)
    extname = get_extname(db)
    nrows = Ref(-1)   # number of rows in the OI-FITS table
    nwaves = Ref(-1)  # number of spectral channels
    for name in fieldnames(typeof(db))
        # Ignore fields not in the specifications.
        spec = get(def.spec, name, nothing)
        spec === nothing && continue

        # Check field value.
        if !isdefined(db, name)
            spec.optional || error("missing mandatory ",
                                   (is_keyword(spec) ? "keyword" : "column"),
                                   " `", name, "` in OI-FITS extension ",
                                   extname)
        else
            check_field_value(nrows, nwaves, getfield(db, name), spec, extname)
        end
    end

    for (key, val) in obj
        # Ignore fields not in the specifications.
        spec = get(def.spec, key, nothing)
        spec === nothing && continue

        # Check value type.
        if spec.type === :L
            is_logical(val) || error("expecting boolean value for `", key,
                                     "` field of OI-FITS extension ", extname)
        elseif spec.type ∈ (:I, :J)
            is_integer(val) || error("expecting integer value for `", key,
                                     "` field of OI-FITS extension ", extname)
        elseif spec.typ ∈ (:D, :E)
            is_float(val) || error("expecting floating-point value for `", key,
                                   "` field of OI-FITS extension ", extname)
        elseif spec.type === :C
            is_complex(val) || error("expecting complex value for `", key,
                                     "` field of OI-FITS extension ", extname)
        elseif spec.type === :A
            is_string(val) || error("expecting string value for `", key,
                                    "` field of OI-FITS extension ", extname)
        else
            error("unknown value type `:", spec.type, "`")
        end

        # Check value dimensions.
        dims = size(val)
        rank = length(dims)
        if is_keyword(spec)
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
                revns::Tuple{Vararg{Integer}},
                tbl::Vector{<:AbstractString})
    for revn in revns
        define(extname, revn, tbl)
    end
end

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
    @noinline function bad_definition(reason::AbstractString,
                                      extname::AbstractString,
                                      revn::Integer, linenum::Integer,
                                      code::AbstractString)
        error(reason, " in definition of OI-FITS extension ", extname,
              " (revision ", revn, ", line ", linenum, "): \"", code, "\"")
    end
    fields = get(FIELDS, extname, Set{Symbol}())
    def = FieldDefinition[]
    keyword = true
    for rownum in 1:length(tbl)
        row = strip(tbl[rownum])
        m = match(r"^([^ ]+) +([^ ]+) +(.*)$", row)
        if m === nothing
            match(r"^-+$", row) === nothing &&
                bad_definition("syntax error", extname, revn, rownum, row)
            keyword = false
            continue
        end
        name = uppercase(m.captures[1])
        symb = to_fieldname(name)
        format = m.captures[2]
        descr = m.captures[3]
        optional = (format[1] == '?')
        i = (optional ? 2 : 1)
        type = get_datatype(format[i])
        type == :UNKNOWN && bad_definition("invalid type letter",
                                           extname, revn, rownum, row)
        if keyword
            length(format) == i ||
                bad_definition("invalid keyword format",
                               extname, revn, rownum, row)
            multiplier = 0
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
                if multiplier === nothing || multiplier < 1
                    bad_definition("invalid multiplier",
                                   extname, revn, rownum, row)
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
        push!(def, FieldDefinition(name, symb, type, multiplier,
                                   optional, descr, units))
        push!(fields, symb)
    end

    # Insert the data-block definition in the global table.
    FORMATS[fmtkey] = DataBlockDefinition(extname, revn, def)
    FIELDS[extname] = fields
    REVISIONS[extname] = max(revn, get(REVISIONS, extname, 0))
    nothing
end

include("formats.jl")

end # module
