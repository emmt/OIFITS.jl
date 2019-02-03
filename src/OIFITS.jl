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
# Copyright (C) 2015-2019: Éric Thiébaut.
#
#------------------------------------------------------------------------------

isdefined(Base, :__precompile__) && __precompile__(true)

module OIFITS

using Compat
using Compat: @debug, @error, @info, @warn

import Base: getindex, setindex!, haskey, keys, show

const Name = (isdefined(Core, :ASCIIString) ? ASCIIString : String)

include("oidata.jl")
include("fix-fitsio.jl")
include("misc.jl")
include("oifile.jl")
include("oiformat1.jl")
include("oiformat2.jl")
include("oipost.jl") # must be *after* oifile.jl and all oiformat*.jl
include("deprecations.jl")

end # module
