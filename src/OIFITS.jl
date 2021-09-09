#
# OIFITS.jl -
#
# Support for OI-FITS (optical interferometry data format) in Julia.
#
#------------------------------------------------------------------------------

module OIFITS

export
    OIArray,
    OICorr,
    OIData,
    OIDataBlock,
    OIFlux,
    OIInsPol,
    OIT3,
    OITarget,
    OIVis,
    OIVis2,
    OIWavelength

import Base:
    convert, copy, eltype, push!, read, haskey, keys, show,
    getindex, setindex!,
    getproperty, setproperty!, propertynames

import FITSIO
import CFITSIO
using FITSIO: FITS, HDU, TableHDU

include("types.jl")
include("formats.jl")
import .Formats: FieldDefinition, get_format
include("methods.jl")

end # module
