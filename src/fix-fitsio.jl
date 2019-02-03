#
# fix-fitsio.jl --
#
# Some fixes and additional routines which should be part of FITSIO.jl package.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015-2019: Éric Thiébaut.
#
#------------------------------------------------------------------------------

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

# BITPIX routines and table.
const _BITPIX = Dict{Cint, DataType}()
for (sym, val, T) in ((:BYTE_IMG,        8,       UInt8),
                      (:SHORT_IMG,      16,       Int16),
                      (:LONG_IMG,       32,       Int32),
                      (:LONGLONG_IMG,   64,       Int64),
                      (:FLOAT_IMG,     -32,       Float32),
                      (:DOUBLE_IMG,    -64,       Float64))
    val = Cint(val)
    _BITPIX[val] = T
    @eval begin
        fits_bitpix(::Type{$T}) = $val
    end
end
fits_bitpix(code::Integer) = get(_BITPIX, Cint(code), Nothing)

# Data type routines and table.
const _DATATYPE = Dict{Cint, DataType}()
const _REVERSE_DATATYPE = Dict{DataType, Cint}()
for (sym, val, T) in ((:TBIT       , Cint(  1), Nothing),
                      (:TBYTE      , Cint( 11), UInt8),
                      (:TSBYTE     , Cint( 12), Int8),
                      (:TLOGICAL   , Cint( 14), Bool),
                      (:TSTRING    , Cint( 16), Name),
                      (:TUSHORT    , Cint( 20), Cushort),          # Uint16
                      (:TSHORT     , Cint( 21), Cshort),           # Int16
                      (:TUINT      , Cint( 30), Cuint),            # Uint32
                      (:TINT       , Cint( 31), Cint),             # Int32
                      (:TULONG     , Cint( 40), Culong),
                      (:TLONG      , Cint( 41), Clong),
                      (:TFLOAT     , Cint( 42), Cfloat),           # Float32
                      (:TLONGLONG  , Cint( 81), Int64),
                      (:TDOUBLE    , Cint( 82), Cdouble),          # Float64
                      (:TCOMPLEX   , Cint( 83), Complex{Cfloat}),  # Complex64
                      (:TDBLCOMPLEX, Cint(163), Complex{Cdouble})) # Complex128
    _DATATYPE[val] = T
    if ! haskey(_REVERSE_DATATYPE, T)
        _REVERSE_DATATYPE[T] = val
        if T == AbstractString
            @eval cfitsio_datatype{S<:AbstractString}(::Type{S}) = $val
        elseif T != Nothing
            @eval cfitsio_datatype(::Type{$T}) = $val
        end
    end
end
cfitsio_datatype(code::Integer) = get(_DATATYPE, code, Nothing)
