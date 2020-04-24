#
# parser.jl --
#
# Parser of OI-FITS data formats.
#
#------------------------------------------------------------------------------

module Parser

using ..OIFITS
using ..OIFITS: to_string, to_integer, to_fieldname

# OIFieldDef is used to store the definition of a keyword/column field.
mutable struct OIFieldDef
    name::String     # Keyword/column name as a string.
    symb::Symbol     # Keyword/column symbolic name.
    iskeyword::Bool  # Is keyword? (otherwise column)
    isoptional::Bool # Optional field?
    multiplier::Int  # Multiplier: 1 for keywords, number of cells for columns
                     # (a negative number -N means an array of N dimensions
                     # each equal to the number of spectral channels.
    type::Symbol    # Data type.
    units::String    # Units.
    descr::String    # Description.
end


# OIDataBlockDef is used to store the definition of data-block.
mutable struct OIDataBlockDef
    extname::String
    fields::Vector{Symbol}        # ordered field symbolic names
    spec::Dict{Symbol,OIFieldDef} # dictionary of field specifications
    function OIDataBlockDef(extname::AbstractString, vect::Vector{OIFieldDef})
        spec = Dict{Symbol,OIFieldDef}()
        fields = Array{Symbol}(undef, length(vect))
        for j in 1:length(vect)
            entry = vect[j]
            fields[j] = entry.symb
            spec[entry.symb] = entry
        end
        new(extname, fields, spec)
    end
end

# OIFormatDef is used to store all the data-block definitions for a given
# revision number.
const OIFormatDef = Dict{String,OIDataBlockDef}

# FIELDS is a dictionary indexed by the data-block name (e,g., "OI_VIS"), each
# entry stores a set of its fields.
const FIELDS = Dict{String,Set{Symbol}}()

# FORMATS table is indexed by a 2-tuple composed by the datablock name
# and by the revision number of the corresponding OI-FITS table.
const FORMATS = Dict{Tuple{String,Int},OIDataBlockDef}()
format_key(name::AbstractString, revn::Integer) =
    (to_string(name), to_integer(revn))

"""
    get_datatype(c)

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
    OIFITS.Parser.get_definition(db)

or

    OIFITS.Parser.get_definition(extname, revn)

yield the format definition for datablock `db` of for the OI-FITS extension
`extname` (e.g., "OI_TARGET") in revision `revn` of OI-FITS standard, throwing
an error if not found.  The returned value is an instance of `OIDataBlockDef`.

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

get_definition(db::OIDataBlock) =
    get_definition(get_extname(db), get_revn(db))

"""
    OIFITS.define(extname, revn, defs)

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
    def = Array{OIFieldDef}(undef, 0)
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
        push!(def, OIFieldDef(name, symb, iskeyword, isoptional, multiplier,
                              type, units, descr))
        push!(fields, symb)
    end

    # Insert the data-block definition in the global table.
    FORMATS[fmtkey] = OIDataBlockDef(extname, def)
    FIELDS[extname] = fields
end

include("formats.jl")

end # module
