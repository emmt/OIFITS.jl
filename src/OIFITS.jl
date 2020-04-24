#
# OIFITS.jl --
#
# Support for OI-FITS (optical interferometry data format) in Julia.
#
#------------------------------------------------------------------------------

module OIFITS

export
    OIArray,
    OICorrelation,
    OIDataBlock,
    OIMaster,
    OIPolarization,
    OISpectrum,
    OIT3,
    OITarget,
    OIVis,
    OIVis2,
    OIWavelength

import Base: getindex, setindex!, haskey, keys, show
using Base: IteratorSize, IteratorEltype

using FITSIO
using FITSIO.Libcfitsio

using FITSIO: libcfitsio, fits_assert_ok, fits_assert_open
using FITSIO.Libcfitsio: fits_get_errstatus

include("types.jl")
include("utils.jl")
include("parser.jl")
include("objects.jl")
include("misc.jl")
include("files.jl")
include("accessors.jl") # must be *after* oifile.jl and all oiformat*.jl
include("deprecations.jl")

end # module
