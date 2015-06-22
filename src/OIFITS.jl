#
# OIFITS.jl --
#
# Support for OI-FITS (optical interferometry data format) in Julia.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015: Éric Thiébaut.
#
#------------------------------------------------------------------------------

module OIFITS

import Base: getindex, setindex!, haskey, keys, start, done, next

include("oidata.jl")
include("fix-fitsio.jl")
include("misc.jl")
include("oifile.jl")
include("oiformat1.jl")
include("oiformat2.jl")
include("oipost.jl") # must be *after* oifile.jl and all oiformat*.jl
include("deprecations.jl")

end # module

# Local Variables:
# mode: Julia
# tab-width: 8
# indent-tabs-mode: nil
# fill-column: 79
# coding: utf-8
# ispell-local-dictionary: "american"
# End:
