#
# oipost.jl --
#
# Post-processing of OI-FITS data-block definitions.  Must be the last
# file included by the main OIFITS.jl source file.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015: Éric Thiébaut.
#
#------------------------------------------------------------------------------

# Automatically define getters from all fields of a data-block.
for dbname in keys(_FIELDS)
    local dbtype = _DATABLOCKS[dbname]
    for symb in _FIELDS[dbname]
        eval(parse("get_$symb(db::$dbtype) = db.contents[:$symb]"))
    end
end

# Define getters which rely on indirections.
get_eff_wave(db::Union{OIVis,OIVis2,OIT3}) = db.ins[:eff_wave]
get_eff_band(db::Union{OIVis,OIVis2,OIT3}) = db.ins[:eff_band]
