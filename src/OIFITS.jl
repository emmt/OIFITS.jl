#
# OIFITS.jl -
#
# Support for OI-FITS (Optical Interferometry data format) in Julia.
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
    IndexStyle,
    convert,
    copy,
    eltype,
    getindex,
    getproperty,
    haskey,
    isempty,
    keys,
    length,
    merge,
    merge!,
    propertynames,
    push!,
    read,
    setindex!,
    setproperty!,
    show,
    size,
    take!

import FITSIO
import CFITSIO
using FITSIO: FITS, HDU, TableHDU

include("types.jl")
include("formats.jl")
import .Formats: FieldDefinition, get_format
include("methods.jl")

end # module
