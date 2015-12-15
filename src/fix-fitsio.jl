#
# fix-fitsio.jl --
#
# Some fixes and additional routines which should be part of FITSIO.jl package.
#
# ----------------------------------------------------------------------------

using FITSIO
using FITSIO.Libcfitsio

import FITSIO: libcfitsio, fits_assert_ok, fits_assert_open
import FITSIO.Libcfitsio: fits_get_errstatus

# The exported functions cfitsio_datatype and fits_bitpix deal with conversion
# between CFITSIO type code or BITPIX value and actual Julia data types.
# They can be used as follows (assuming `T` is a Julia data type, while
# `code` and `bitpix` are integers):
#
#     cfitsio_datatype(T) --------> code (e.g., TBYTE, TFLOAT, etc.)
#     cfitsio_datatype(code) -----> T
#
#     fits_bitpix(T) -------------> bitpix (e.g., BYTE_IMG, FLOAT_IMG, etc.)
#     fits_bitpix(bitpix) --------> T
#
export cfitsio_datatype, fits_bitpix

# The following table gives the correspondances between CFITSIO "types",
# the BITPIX keyword and Julia types.
#
#     -------------------------------------------------
#     BITPIX   CFISTIO         Julia     Comments
#     -------------------------------------------------
#              int             Cint
#              long            Clong
#              LONGLONG        Int64     64-bit integer
#     -------------------------------------------------
#        8     BYTE_IMG        UInt8
#       16     SHORT_IMG       Int16
#       32     LONG_IMG        Int32
#       64     LONGLONG_IMG    Int64
#      -32     FLOAT_IMG       Float32
#      -64     DOUBLE_IMG      Float64
#     -------------------------------------------------
#              TBIT
#              TBYTE           Cuchar = UInt8
#              TSBYTE          Cchar = Int8
#              TLOGICAL        Bool
#              TSHORT          Cshort
#              TUSHORT         Cushort
#              TINT            Cint
#              TUINT           Cuint
#              TLONG           Clong
#              TLONGLONG       Int64
#              TULONG          Culong
#              TFLOAT          Cfloat
#              TDOUBLE         Cdouble
#              TCOMPLEX        Complex{Cfloat}
#              TDBLCOMPLEX     Complex{Cdouble}
#     -------------------------------------------------

# Conversion to a C `int`:
cint(x) = convert(Cint, x)
cint(x::Cint) = x

# BITPIX routines and table.
const _BITPIX = Dict{Cint, DataType}()
for (sym, val, T) in ((:BYTE_IMG,        8,       UInt8),
                      (:SHORT_IMG,      16,       Int16),
                      (:LONG_IMG,       32,       Int32),
                      (:LONGLONG_IMG,   64,       Int64),
                      (:FLOAT_IMG,     -32,       Float32),
                      (:DOUBLE_IMG,    -64,       Float64))
    val = cint(val)
    _BITPIX[val] = T
    @eval begin
        fits_bitpix(::Type{$T}) = $val
    end
end
fits_bitpix(code::Integer) = get(_BITPIX, cint(code), Void)

# Data type routines and table.
const _DATATYPE = Dict{Cint, DataType}()
for (sym, val, T) in ((:TBIT       ,   1, Void),
                      (:TBYTE      ,  11, UInt8),
                      (:TSBYTE     ,  12, Int8),
                      (:TLOGICAL   ,  14, Bool),
                      (:TSTRING    ,  16, AbstractString),
                      (:TUSHORT    ,  20, Cushort),          # Uint16
                      (:TSHORT     ,  21, Cshort),           # Int16
                      (:TUINT      ,  30, Cuint),            # Uint32
                      (:TINT       ,  31, Cint),             # Int32
                      (:TULONG     ,  40, Culong),
                      (:TLONG      ,  41, Clong),
                      (:TFLOAT     ,  42, Cfloat),           # Float32
                      (:TLONGLONG  ,  81, Int64),
                      (:TDOUBLE    ,  82, Cdouble),          # Float64
                      (:TCOMPLEX   ,  83, Complex{Cfloat}),  # Complex64
                      (:TDBLCOMPLEX, 163, Complex{Cdouble})) # Complex128
    val = cint(val)
    _DATATYPE[val] = T
    if T == AbstractString
        @eval cfitsio_datatype{S<:AbstractString}(::Type{S}) = $val
    elseif T != Void
        @eval cfitsio_datatype(::Type{$T}) = $val
    end
end
cfitsio_datatype(code::Integer) = get(_DATATYPE, cint(code), Void)
