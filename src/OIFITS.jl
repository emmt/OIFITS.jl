#
# OIFITS.jl --
#
# Support for OI-FITS (optical interferometry data format) in Julia.
#
#------------------------------------------------------------------------------

isdefined(Base, :__precompile__) && __precompile__(true)

module OIFITS

import Base: getindex, setindex!, haskey, keys, show

using FITSIO
using FITSIO.Libcfitsio

import FITSIO: libcfitsio, fits_assert_ok, fits_assert_open
import FITSIO.Libcfitsio: fits_get_errstatus

include("types.jl")
include("oidata.jl")
include("bitpix.jl")
include("misc.jl")
include("oifile.jl")
include("oiformat1.jl")
include("oiformat2.jl")
include("oipost.jl") # must be *after* oifile.jl and all oiformat*.jl
include("deprecations.jl")

end # module
