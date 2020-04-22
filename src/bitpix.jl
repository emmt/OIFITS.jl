#
# bitpix.jl --
#
# Convert between Julia types and FITS BITPIX or FITSIO column-type
# codes.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015-2020, Éric Thiébaut.
#
#------------------------------------------------------------------------------
#
# The functions `type_to_coltype`, `coltype_to_type`, `bitpix_to_type` and
# `type_to_bitpix` deal with conversion between CFITSIO type code or BITPIX
# value and actual Julia data types.  They can be used as follows (assuming `T`
# is a Julia data type, while `coltype` and `bitpix` are integers):
#
#     type_to_coltype(T) ------------> coltype (e.g., TBYTE, TFLOAT, etc.)
#     coltype_to_type(coltype) ------> T
#
#     type_to_bitpix(T) -------------> bitpix (e.g., BYTE_IMG, FLOAT_IMG, etc.)
#     bitpix_to_type(bitpix) --------> T
#
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

"""
    type_to_bitpix(T) -> val

yields the FITS BITPIX value `val` which corresponds to Julia type `T`.

"""
type_to_bitpix(::Type{UInt8})   =   8
type_to_bitpix(::Type{Int16})   =  16
type_to_bitpix(::Type{Int32})   =  32
type_to_bitpix(::Type{Int64})   =  64
type_to_bitpix(::Type{Float32}) = -32
type_to_bitpix(::Type{Float64}) = -64

"""
    bitpix_to_type(val) -> T

yields the Julia type `T` which corresponds to the BITPIX value `val`.

"""
bitpix_to_type(val::Integer) =
    (val ==   8 ? UInt8   :
     val ==  16 ? Int16   :
     val ==  32 ? Int32   :
     val ==  64 ? Int64   :
     val == -32 ? Float32 :
     val == -64 ? Float64 : bad_bitpix(val))

@noinline  bad_bitpix(val) =
    error(string("bad BITPIX value ", val))

"""
    coltype_to_type(val) -> T

yields the Julia type `T` which best corresponds to the FITSIO column type
`val` or `Nothing` if .

"""
coltype_to_type(val::Integer) =
    (#   ==   1 ? ...              : # TBIT
     val ==  11 ? UInt8            : # TBYTE
     val ==  12 ? Int8             : # TSBYTE
     val ==  14 ? Bool             : # TLOGICAL
     val ==  16 ? String           : # TSTRING
     val ==  20 ? Cushort          : # TUSHORT      (Uint16?)
     val ==  21 ? Cshort           : # TSHORT       (Int16?)
     val ==  30 ? Cuint            : # TUINT        (Uint32?)
     val ==  31 ? Cint             : # TINT         (Int32?)
     val ==  40 ? Culong           : # TULONG
     val ==  41 ? Clong            : # TLONG
     #   ==  42 ? ...              : # TINT32BIT
     val ==  42 ? Cfloat           : # TFLOAT      (Float32?)
     val ==  81 ? Int64            : # TLONGLONG
     val ==  82 ? Cdouble          : # TDOUBLE     (Float64?)
     val ==  83 ? Complex{Cfloat}  : # TCOMPLEX    (ComplexF32?)
     val == 163 ? Complex{Cdouble} : # TDBLCOMPLEX (ComplexF64?)
     Nothing)
