#
# fix-fitsio.jl --
#
# Some fixes and additional routines which should be part of FITSIO.jl package.
#
# ----------------------------------------------------------------------------

using FITSIO
const libcfitsio = FITSIO.libcfitsio

export fits_datatype, fits_bitpix;

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
#        8     BYTE_IMG        Uint8
#       16     SHORT_IMG       Int16
#       32     LONG_IMG        Int32
#       64     LONGLONG_IMG    Int64
#      -32     FLOAT_IMG       Float32
#      -64     DOUBLE_IMG      Float64
#     -------------------------------------------------
#              TBIT
#              TBYTE           Cuchar = Uint8
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


cint(x) = convert(Cint, x)
cint(x::Cint) = x

const _BITPIX_TABLE = Dict{Cint, DataType}()
for (sym, val, T) in ((:BYTE_IMG,        8,       Uint8),
                      (:SHORT_IMG,      16,       Int16),
                      (:LONG_IMG,       32,       Int32),
                      (:LONGLONG_IMG,   64,       Int64),
                      (:FLOAT_IMG,     -32,       Float32),
                      (:DOUBLE_IMG,    -64,       Float64))
    val = cint(val)
    _BITPIX_TABLE[val] = T
    @eval begin
        const $sym = $val
        fits_bitpix(::Type{$T}) = $val
    end
end
fits_bitpix(code::Integer) = get(_BITPIX_TABLE, cint(code), Nothing)

const _DATATYPE_TABLE = Dict{Cint, DataType}()
for (sym, val, T) in ((:TBIT       ,   1, Nothing),
                      (:TBYTE      ,  11, Uint8),
                      (:TSBYTE     ,  12, Int8),
                      (:TLOGICAL   ,  14, Bool),
                      (:TSTRING    ,  16, String),
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
    _DATATYPE_TABLE[val] = T
    @eval const $sym = $val
    if T == String
        @eval fits_datatype{S<:String}(::Type{S}) = $val
    elseif T != Nothing
        @eval fits_datatype(::Type{$T}) = $val
    end
end
fits_datatype(code::Integer) = get(_DATATYPE_TABLE, cint(code), Nothing)

# The function `fits_read_tdim()` returns the dimensions of a table column
# in a binary table. Normally this information is given by the TDIMn
# keyword, but if this keyword is not present then this routine returns [r]
# with r equals to the repeat count in the TFORM keyword.
let fn, T
    if Int == Clong
        T = Clong
        ffgtdm = "ffgtdm"
        ffgtcl = "ffgtcl"
        ffeqty = "ffeqty"
    else
        T = Clonglong
        ffgtdm = "ffgtdmll"
        ffgtcl = "ffgtclll"
        ffeqty = "ffeqtyll"
    end
    @eval begin
        function fits_get_coltype(ff::FITSFile, colnum::Integer)
            typecode = Cint[0]
            repeat = $T[0]
            width = $T[0]
            status = Cint[0]
            ccall(($ffgtcl,libcfitsio), Cint,
                  (Ptr{Void}, Cint, Ptr{$T}, Ptr{$T}, Ptr{$T}, Ptr{Cint}),
                  ff.ptr, colnum, typecode, repeat, width, status)
            fits_assert_ok(status[1])
            return (Int(typecode[1]), Int(repeat[1]), Int(width[1]))
        end
        function fits_get_eqcoltype(ff::FITSFile, colnum::Integer)
            typecode = Cint[0]
            repeat = $T[0]
            width = $T[0]
            status = Cint[0]
            ccall(($ffeqty,libcfitsio), Cint,
                  (Ptr{Void}, Cint, Ptr{$T}, Ptr{$T}, Ptr{$T}, Ptr{Cint}),
                  ff.ptr, colnum, typecode, repeat, width, status)
            fits_assert_ok(status[1])
            return (Int(typecode[1]), Int(repeat[1]), Int(width[1]))
        end
        function fits_read_tdim(ff::FITSFile, colnum::Integer)
            naxes = Array($T, 99)
            naxis = Cint[0]
            status = Cint[0]
            ccall(($ffgtdm,libcfitsio), Cint,
                  (Ptr{Void}, Cint, Cint, Ptr{Cint}, Ptr{$T}, Ptr{Cint}),
                  ff.ptr, colnum, length(naxes), naxis, naxes, status)
            fits_assert_ok(status[1])
            return naxes[1:naxis[1]]
        end
    end
end

# Local Variables:
# mode: Julia
# tab-width: 8
# indent-tabs-mode: nil
# fill-column: 79
# coding: utf-8
# ispell-local-dictionary: "american"
# End:
