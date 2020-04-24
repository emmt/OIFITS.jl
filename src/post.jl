#
# post.jl --
#
# Post-processing of OI-FITS data-block definitions.  Must be the last file
# included by the main OIFITS.jl source file.
#
#------------------------------------------------------------------------------

# Automatically define getters from all fields of a data-block.
for extname in EXTNAMES
    let T = get_datablock_type(extname)
        for symb in Parser.FIELDS[extname]
            eval(Meta.parse("get_$symb(db::$T) = db.$symb"))
        end
    end
end

# Define getters which rely on indirections.
get_eff_wave(db::Union{OIVis,OIVis2,OIT3}) = db.instr.eff_wave
get_eff_band(db::Union{OIVis,OIVis2,OIT3}) = db.instr.eff_band

"""

Assuming `db` is an instance of `OITarget`, `OIVis`, `OIVis2` or `OIT3`, then:

    get_target_id(db)

yields the "TARGET_ID" column of `db` which is an array of integers.

""" get_target_id
