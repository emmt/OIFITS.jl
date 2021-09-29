#
# OIFITS.jl -
#
# Support for OI-FITS (Optical Interferometry data format) in Julia.
#
#------------------------------------------------------------------------------

module OIFITS

export
    OIDataSet,
    OIDataBlock,
    OITargetEntry,
    OI_ARRAY,
    OI_CORR,
    OI_FLUX,
    OI_INSPOL,
    OI_T3,
    OI_TARGET,
    OI_VIS,
    OI_VIS2,
    OI_WAVELENGTH

import Base:
    IndexStyle,
    axes,
    convert,
    copy,
    eachindex,
    eltype,
    empty!,
    firstindex,
    get,
    getindex,
    getproperty,
    haskey,
    isempty,
    iterate,
    keys,
    lastindex,
    length,
    merge,
    merge!,
    ndims,
    propertynames,
    push!,
    read,
    setindex!,
    setproperty!,
    show,
    size,
    take!,
    values,
    write

using Base: @propagate_inbounds

import FITSIO
import CFITSIO
using FITSIO: FITS, HDU, TableHDU, FITSHeader

include("types.jl")
include("formats.jl")
import .Formats: FieldDefinition, get_format
include("methods.jl")
include("io.jl")

end # module
