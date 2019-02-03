# OIFITS.jl

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://travis-ci.org/emmt/OIFITS.jl.svg?branch=master)](https://travis-ci.org/emmt/OIFITS.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/OIFITS.jl?branch=master)](https://ci.appveyor.com/project/emmt/OIFITS-jl/branch/master)

The `OIFITS.jl` package provides support for OI-FITS data in Julia language.


## OI-FITS Summary

OI-FITS is a standard to store optical interferometry data as a collection of
data-blocks.  In the first version of the standard (see [Ref. 1](#references)),
the available data-blocks are:

* `OI_TARGET` provides a list of observed targets;
* `OI_ARRAY` describes a given array of stations;
* `OI_WAVELENGTH` describes a given instrument (notably the effective
  wavelengths and bandwidths of its spectral channels);
* `OI_VIS` contains complex visibility data;
* `OI_VIS2` contains squared visibility (powerspectrum) data;
* `OI_T3` contains triple product (bispectrum) data.

These data-blocks, are stored as binary tables in a FITS data file.  The
support for the actual FITS files is provided by the
[`FITSIO.jl`](https://github.com/JuliaAstro/FITSIO.jl) package.


## Installation

OIFITS is a [registered Julia package](http://pkg.julialang.org/), the
installation is as simple as:
```julia
Pkg.add("OIFITS")
Pkg.update()
```
The last command `Pkg.update()` may be unnecessary.


## Typical usage

Loading an OI-FITS data file:
```julia
using OIFITS
master = OIFITS.load("testdata.oifits")
```

To iterate through all data-blocks:
```julia
for db in master
    dbname = OIFITS.get_dbname(db)
    revn = OIFITS.get_revn(db)
    println("Data block is $dbname, revision $revn")
end
```

To iterate through a sub-set of the data-blocks (here the complex visibility
data, the powerspectrum data and the bispectrum data):
```julia
for db in OIFITS.select(master, "OI_VIS", "OI_VIS2", "OI_T3")
    dbname = OIFITS.get_dbname(db)
    n = length(OIFITS.get_time(db))
    println("Data block is $dbname, number of exposures is $n")
end
```


## Accessor functions

Any OI-FITS field (keyword/column) of a given data-block can be retrieved
via an accessor whose name has suffix `OIFITS.get_` followed by the name of
the field (in lower case letters and with all non-letter and all non-digit
letters replaced by the underscore character `'_'`).  A notable exception is
the revision number corresponding to the keyword "OI_REVN" which is
retrieved with the method `OIFITS.get_revn()`.  For instance:

```julia
OIFITS.get_revn(db)      # get the revison number of the format (OI_REVN)
OIFITS.get_eff_wave(db)  # get effective wavelengths (EFF_WAVE)
OIFITS.get_eff_band(db)  # get effective bandwidths (EFF_BAND)
OIFITS.get_ucoord(db)    # get the U coordinates of the data (UCOORD)
```
Of course, getting a given field must make sense.  For instance,
`OIFITS.get_eff_wave()` can be applied on any `OI_WAVELENGTH` data-blocks
but also on data-blocks which contains interferometric data such as
`OI_VIS`, `OI_VIS2`, `OI_T3`, *etc.* but not on other data-blocks like
`OI_TARGET`.


## Reading data

To load the contents of an OI-FITS file in memory, use:
```julia
master = OIFITS.load(filename)
```
where `filename` is the name of the file and the returned value, `master`,
contains all the OI-FITS data-blocks of the file.  You may have the names
of the data blocks printed as they get read with keyword `quiet=false`:
```julia
master = OIFITS.load(filename, quiet=false)
```
If you already have a `FITS` handle to the data, you can use it as the
argument to `OIFITS.load` in place of the file name.


## Constructors

It is possible to build OI-FITS data-blocks individually.  The general
syntax is:
```julia
OIFITS.new_XXX(KEY1=VAL1, KEY2=VAL2, ...)
```
where `XXX` is the type of the data-block and `KEYn=VALn` constructions
give the fields of the data-block and their values.  The names of the
fields follow the same convention as for the field accessors.

Available data-block constructors are:

* `OIFITS.new_target` => `OI_TARGET`
* `OIFITS.new_array` => `OI_ARRAY`
* `OIFITS.new_wavelength` => `OI_WAVELENGTH`
* `OIFITS.new_vis`  => `OI_VIS`
* `OIFITS.new_vis2` => `OI_VIS2`
* `OIFITS.new_t3`   => `OI_T3`

When defining a new data-block, all mandatory fields must be provided.
For instance, to create an `OI_WAVELENGTH` data-block:
```julia
µm = 1e-6  # all values are in SI units in OI-FITS
db = OIFITS.new_wavelength(insname="Amber",
                           eff_wave=[1.4µm,1.6µm,1.8µm],
                           eff_band=[0.2µm,0.2µm,0.2µm])
```
Note that the revision number (`revn=...`) can be omitted; by default, the
highest defined revision will be used.

A consistent set of OI-FITS data-blocks is made of: exactly one `OI_TARGET`
data-block, one or more `OI_WAVELENGTH` data-blocks, one or more `OI_ARRAY`
data-blocks and any number of data-blocks with interferometric data
(`OI_VIS`, `OI_VIS2` or `OI_T3`).  These data-blocks must be stored in a
container created by:
```julia
master = OIFITS.new_master()
```
Then, call:
```julia
OIFITS.attach(master, db)
```
to attach all data-block `db` to the OI-FITS container (in any order).
Finally, you must call:
```julia
OIFITS.update(master)
```
to update internal information such as links between data-blocks with
interferometric data and the related instrument (`OI_WAVELENGTH`
data-block) and array of stations (`OI_ARRAY` data-block).  If you do not
do that, then applying some accessors may not work, *e.g.*
`OIFITS.get_eff_wave()` on a data-blocks with interferometric data.

To read an OI-FITS data-block from the HDU of a FITS file:
```julia
db = OIFITS.read_datablock(hdu)
```
where `hdu` is a FITS `HDU` handle.  The result may be `nothing` if the
current HDU does not contain an OI-FITS data-block.


## Miscellaneous functions

OI-FITS implements some useful functions which can be used to deal with
FITS file (not just OI-FITS ones).  These functions could be part of `FITSIO`
package.


### Retrieving information from the header of a FITS HDU

The header of a FITS HDU can be read with the function:
```julia
fts = FITS(filename)
hdr = FITSIO.read_header(fts[1])
```
which returns an indexable and iterable object, here `hdr`.  The keys of
`hdr` are the FITS keywords of the header.  For instance:
```julia
keys(hdr)          # yield an iterator on the keys of hdr
collect(keys(hdr)) # yield all the keys of hdr
haskey(hdr, key)   # check whether key is present
hdr[key]           # retrieve the contents associated with the key
```
For commentary FITS keywords (`"HISTORY"` or `"COMMENT"`), there is no
value, just a comment but there may be any number of these *commentary*
keywords.  Other keywords must be unique and thus have a scalar value.  Use
`get_comment` to retrieve the comment of a FITS keyword:
```julia
get_comment(hdr, key)keys(hdr)          # yield an iterator on the keys of hdr
collect(keys(hdr)) # yield all the keys of hdr
haskey(hdr, key)   # check whether key is present
hdr[key]           # retrieve the contents associated with the key
```

*OIFITS* provides method `OIFITS.get_value()` and `OIFITS.get_comment()`
method to retrieve the value and comment (respectively) of a FITS keyword
with type checking and, optionaly, let you provide a default value if the
keyword is absent:
```julia
val = OIFITS.get_value(hdr, key)
val = OIFITS.get_value(hdr, key, def)
com = OIFITS.get_comment(hdr, key)
com = OIFITS.get_comment(hdr, key, def)
```
To retrieve a value and make sure it has a specific type, the following
methods are available:
```julia
OIFITS.get_logical(hdr, "SIMPLE")
OIFITS.get_integer(hdr, "BITPIX")
OIFITS.get_real(hdr, "BSCALE")
OIFITS.get_string(hdr, "XTENSION")
```
which throw an error if the keyword is not present and perform type
checking and conversion if allowed.  It is also possible to supply a
default value to return if the keyword is not present:
```julia
bscale = OIFITS.get_real(hdr, "BSCALE", 1.0)
bzero = OIFITS.get_real(hdr, "BZERO", 0.0)
xtension = OIFITS.get_string(hdr, "XTENSION", "IMAGE")
```

The function:
```julia
OIFITS.get_hdutype(hdr)
```
returns the HDU type as a Symbol, `:image_hdu` for an image, `:ascii_table`
for an ASCII table, `:binary_table` for a binary table, and `:unknown`
otherwise.

For a FITS table, the function:
```julia
OIFITS.get_dbtype(hdr)
```
returns the OI-FITS data-block type as a Symbol like `:OI_TARGET`,
`:OI_WAVELENGTH`, *etc.*


### Reading FITS tables

In addition to the method `read(tbl::TableHDU, colname::String)`
provided by FITSIO for reading a specific column of a FITS table, the
low-level function:
```julia
OIFITS.read_column(ff::FITSFile, colnum::Integer)
```
returns a Julia array with the contents of the `colnum`-th column of the
current HDU in FITS file handle `ff`.  The current HDU must be a FITS table
(an ASCII or a binary one).  The last dimension of the result corresponds
to the rows of the table.  It is also possible to read all the table:
```julia
OIFITS.read_table(ff::FITSFile)
OIFITS.read_table(hdu::Union(TableHDU,ASCIITableHDU))
```
or at high-level:
```julia
read(hdu::Union(TableHDU,ASCIITableHDU))
```
The result is a dictionary whose keys are the names of the columns (in
uppercase letters and with trailing spaces removed).  If a column has given
units, the units are stored in the dictionary with suffix `".units"`
appended to the column name.  For instance, the units for column `"TIME"`
are accessible with key `"TIME.units"`.


### FITS and Julia types conversion

The functions `cfitsio_datatype()` and `fits_bitpix()` deal with conversion
between CFITSIO type code or BITPIX value and actual Julia data types.
They can be used as follows (assuming `T` is a Julia data type, while
`code` and `bitpix` are integers):
```julia
cfitsio_datatype(T) --------> code (e.g., TBYTE, TFLOAT, etc.)
cfitsio_datatype(code) -----> T

fits_bitpix(T) -------------> bitpix (e.g., BYTE_IMG, FLOAT_IMG, etc.)
fits_bitpix(bitpix) --------> T
```

The functions `fits_get_coltype()` and `fits_get_eqcoltype()` yield the
data type, repeat count and width in bytes of a given column, their
prototypes are:
```julia
(code, repcnt, width) = fits_get_coltype(ff::FITSFile, colnum::Integer)
(code, repcnt, width) = fits_get_eqcoltype(ff::FITSFile, colnum::Integer)
```
with `colnum` the column number, `code` the CFITSIO column type (call
`cfitsio_datatype(code)` to convert it to a Julia type) of the elements in
this column, `repcnt` and `width` the repeat count and width of a cell in
this column.  The difference between `fits_get_coltype()` and
`fits_get_eqcoltype()` is that the former yields the column type as it is
stored in the file, while the latter yields the column type after automatic
scaling by the values `"TSCALn"` and `"TZEROn"` keywods if present (with
`n` the column number).  Note that reading the column data with
`fits_read_col()` or `fitsio_read_column()` automatically apply this kind
of scaling.

To retrieve the dimensions of the cells in a given column, call the
function `fits_read_tdim()`, its prototype is:
```julia
dims = fits_read_tdim(ff::FITSFile, colnum::Integer)
```
where `dims` is a vector of integer dimensions.


## Credits

The developments of this package has received funding from the European
Community's Seventh Framework Programme (FP7/2013-2016) under Grant
Agreement 312430 (OPTICON).


## References

1. Pauls, T. A., Young, J. S., Cotton, W. D., & Monnier, J. D. "A data exchange
   standard for optical (visible/IR) interferometry." Publications of the
   Astronomical Society of the Pacific, vol. 117, no 837, p. 1255 (2005).
   [[pdf]](http://arxiv.org/pdf/astro-ph/0508185)

2. Duvert, G., Young, J., & Hummel, C. "OIFITS 2: the 2nd version of the Data
   Exchange Standard for Optical (Visible/IR) Interferometry." arXiv preprint
   [[arXiv:1510.04556v2.04556]](http://arxiv.org/abs/1510.04556v2).
