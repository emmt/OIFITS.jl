#
# accessors.jl --
#
# Accessor functions for OI-FITS objects.  Most of these are kept for backward
# compatibility but superseded by the `obj.key` syntax.
#
# This file must be included after OI-fITS formats definitions.
#
#------------------------------------------------------------------------------

# Automatically define getters from all fields of a data-block.
for extname in EXTNAMES
    let T = get_datablock_type(extname)
        for symb in Builder.FIELDS[extname]
            eval(Meta.parse("get_$symb(db::$T) = db.$symb"))
        end
    end
end

# Define getters which rely on indirections.
get_eff_wave(db::Union{OIVis,OIVis2,OIT3}) = db.instr.eff_wave
get_eff_band(db::Union{OIVis,OIVis2,OIT3}) = db.instr.eff_band

"""
    OIFITS.get_target_id(obj)

yields the "TARGET_ID" field of OI-FITS data-block object `obj` (an instance of
`OITarget`, `OIVis`, `OIVis2`, `OIT3`, `OIFlux` or `OIPolarization`).

!!! note
    This method is to be deprecated, use `obj.target_id` instead.

""" get_target_id

"""
    OIFITS.get_target(obj)

yields the `OITarget` data-block associated with OI-FITS data-block object
`obj`, if any, `nothing` otherwise.

!!! note
    This method is to be deprecated, use `obj.target` or `obj.target[arrname]`
    instead.

"""
get_target(obj::OIMaster{T}) where {T} =
    obj.target :: Union{Nothing,OITarget{T}}
#get_target(obj::OITarget) = obj
get_target(obj::OIDataBlock) = nothing
get_target(obj::Union{OIVis,OIVis2,OIT3,OIFlux,OIPolarization}) =
    get_target(obj.owner)

"""
    OIFITS.get_array(obj)

yields the `OIArray` data-block associated with OI-FITS data-block object
`obj`, if any, `nothing` otherwise.

    OIFITS.get_array(obj, arrname)

yields the `OIArray` data-block whose name is `arrname` in OI-FITS matser
object `obj` or `nothing` if not found.

!!! note
    This method is to be deprecated, use `obj.array` or `obj.array[arrname]`
    instead.

"""
get_array(obj::OIMaster{T}, arrname::AbstractString) where {T} =
    get(get_array(obj), fix_name(arrname),
        nothing) :: Union{Nothing,OIArray{T}}
get_array(obj::OIDataBlock) = nothing
get_array(obj::OIArray) = obj
get_array(obj::Union{OIMaster,OIVis,OIVis2,OIT3,OIFlux,OIPolarization}) =
    obj.array

"""
    OIFITS.get_instrument(obj)

yields the `OIWavelength` data-block associated with OI-FITS data-block object
`obj`, if any, `nothing` otherwise.

    OIFITS.get_instrument(obj, arrname)

yields the `OIWavelength` data-block whose name is `arrname` in OI-FITS matser
object `obj` or `nothing` if not found.

!!! note
    This method is to be deprecated, use `obj.instr` or `obj.instr[arrname]`
    instead.

"""
get_instrument(obj::OIMaster{T}, insname::AbstractString) where {T} =
    get(get_instrument(obj), fix_name(insname),
        nothing) :: Union{Nothing,OIWavelength{T}}
get_instrument(obj::OIDataBlock) = nothing
get_instrument(obj::OIWavelength) = obj
get_instrument(obj::Union{OIMaster,OIVis,OIVis2,OIT3,OIFlux,OIPolarization}) =
    obj.instr

"""
    OIFITS.get_targets(obj)

yields the array of names of the targets defined in OI-FITS master object `obj`.

"""
get_targets(obj::OIMaster) =
    (tgt = obj.target) === nothing ? String[] : tgt.target

"""
    OIFITS.get_arrays(obj)

yields the names of the interferometric arrays defined in OI-FITS master object
`obj`.

"""
get_arrays(obj::OIMaster) = collect(keys(obj.array))

"""
    OIFITS.get_instruments(obj)

yields the names of the instruments defined in OI-FITS master object `obj`.

"""
get_instruments(obj::OIMaster) = collect(keys(obj.instr))
