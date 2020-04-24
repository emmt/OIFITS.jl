# OIFITS.jl

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |


The `OIFITS.jl` package provides support for OI-FITS data in Julia language.


## OI-FITS Summary

OI-FITS is a standard to store optical interferometry data as a collection of
data-blocks.  In the second revision of the standard (see [Ref. 1](#references)
and [Ref. 2](#references)), the available data-blocks are:

* `OITarget` provides a list of observed targets (`OI_TARGET` extension in
  OI-FITS files);
* `OIArray` describes a given array of stations (`OI_ARRAY` extension in
  OI-FITS files);
* `OIWavelength` describes a given instrument notably the effective
  wavelengths and bandwidths of its spectral channels (`OI_WAVELENGTH`
  extension in OI-FITS files);
* `OIVis` contains complex visibility data (`OI_VIS` extension in OI-FITS
  files);
* `OIVis2` contains squared visibility (powerspectrum) data (`OI_VIS2`
  extension in OI-FITS files);
* `OIT3` contains triple product (bispectrum) data (`OI_T3` extension in
  OI-FITS files);
* `OISpectrum` contains spectrum data (`OI_SPECTRUM` extension in OI-FITS
  files);
* `OIPolarization` contains instrumental polarization (`OI_INSPOL`
  extension in OI-FITS files);
* `OICorrelation` contains correlation data (`OI_CORR` extension in OI-FITS
  files).

These data-blocks, are stored as binary tables in a FITS data file.  The
support for the actual FITS files is provided by the
[`FITSIO.jl`](https://github.com/JuliaAstro/FITSIO.jl) package.


## Installation

OIFITS is a [registered Julia package](http://pkg.julialang.org/), the
installation is as simple as:

```julia
using Pkg
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

To iterate through all OI-FITS data-blocks:

```julia
for db in master
    extname = db.extname
    revn = db.revn
    println("Data-block is $extname, revision $revn")
end
```

To iterate through a sub-set of the data-blocks (here the complex visibility
data, the powerspectrum data and the bispectrum data):

```julia
for db in OIFITS.select(master, "OI_VIS", "OI_VIS2", "OI_T3")
    extname = db.extname
    n = length(db.time)
    println("Data-block is $extname, number of exposures is $n")
end
```

Any OI-FITS field (keyword/column) of a given data-block can be accessed by the
`db.key` syntax where `key` is the name of the field in lower case letters and
with non-alphanumeric letters replaced by an underscore character `_`.  For
instance `db.int_time` yields the value of `INT_TIME` the integration times of
the measurements.  The revision number corresponding to the keyword `OI_REVN`
is however accessed as `db.revn`, this is the only exception.  Other fields are
also accessible via this syntax:

- `db.extname` yields the OI-FITS name of the extension corresponding to the
  data-block `db` (for all data-block types);

- `db.array` yields the `OI_ARRAY` data-block associated with data-block `db`
  (for `OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_SPECTRUM` and `OI_POLARIZATION`
  data-block);

- `db.instr` yields the `OI_WAVELENGTH` data-block associated with data-block
  `db` (for `OI_VIS`, `OI_VIS2`, `OI_T3` and `OI_SPECTRUM` data-block);

- `db.correl` yields the `OI_CORREL` data-block associated with data-block `db`
  (for `OI_VIS`, `OI_VIS2`, `OI_T3` and `OI_SPECTRUM` data-block).

Of course, getting a given field must make sense.  For instance, `db.eff_wave`
is only possible for `OI_WAVELENGTH` data-blocks but not for an `OI_TARGET`
data-block.  The dot notation can be however be chained and:

```julia
db.instr.eff_wave
```

can be used to access the effective wavelengths of the measurements in `db` via
the instrument associated to `db`.

The dot notation can also be used on the *master* object storing all the
data-blocks:

- `master.target` is the `OI_TARGET` data-block of the OI-FITS structure;

- `master.instr` is a dictionary of `OI_WAVELENGTH` data-blocks indexed by the
  instrument names: `master.instr[insname]` yields the `OI_WAVELENGTH`
  data-block named `insname` or `nothing` if no such instrument exists;

- `master.array` is a dictionary of `OI_ARRAY` data-blocks indexed by the
  telescope array names: `master.array[arrname]` yields the `OI_ARRAY`
  data-block named `arrname` or `nothing` if no such array exists;

- `master.correl` is a dictionary of `OI_CORREL` data-blocks indexed by the
  correlation data names: `master.correl[corrname]` yields the `OI_CORREL`
  data-block named `corrname` or `nothing` if no such correlation data exists.

As show in above examples, the *master* object storing OI-FITS data-blocks can
be used as an iterator over all the stored data-blocks.  It can also be used as
a vector of data-blocks indexed by the data-block number.  For instance:

```julia
for i in eachindex(master)
    let db = master[i]
        extname = db.extname
        revn = db.revn
        println("Data-block is $extname, revisions $revn")
    end
end
```

However, beware that the *master* object storing OI-FITS data-blocks is not a
simple vector, it does many things *under the hood* to maintain the consistency
of the structure (for instance links between different data-blocks).  In order
to append data-blocks to the ones already stored, call:

```julia
push!(master, db)
```

Argument `db` above can be a variable number of data-blocks, a tuple or an
array of data-blocks.  If there many data-blocks, it is more efficient to push
them all at the same time.


## Deprecated accessor functions

Although this is considered as deprecated, OI-FITS fields of a given data-block
can also be retrieved via an accessor whose name has suffix `OIFITS.get_`
followed by the name of the field in lower case letters and with all
non-alphanumeric letters replaced by an underscore character `_`).  A notable
exception is the revision number corresponding to the keyword "OI_REVN" which
is retrieved with the method `OIFITS.get_revn()`.  For instance:

```julia
OIFITS.get_revn(db)      # get the revison number of the format (OI_REVN)
OIFITS.get_eff_wave(db)  # get effective wavelengths (EFF_WAVE)
OIFITS.get_eff_band(db)  # get effective bandwidths (EFF_BAND)
OIFITS.get_ucoord(db)    # get the U coordinates of the data (UCOORD)
```



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
push!(master, db)
```

to attach all data-block `db` to the OI-FITS container (in any order).

To read OI-FITS data from a Header Data Units (HDU) of a FITS file, call:

```julia
dat = OIFITS.read_datablock(hdu)
```

where `hdu` is a FITS HDU.  The result may be `nothing` if the current HDU does
not contain an OI-FITS data-block; otherwise the result is a 3-tuple `(extname,
revn, dict)` with the name of the FITS extension, the OI-FITS revision number
and a dictionary of the OI-FITS keywords and columns. These can be directly
provided to `OIFITS.build_datablock` to build an instance of `OIDataBlock`:

```julia
db = OIFITS.build_datablock(extname, revn, dict)
```


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
OIFITS.get_float(hdr, "BSCALE")
OIFITS.get_string(hdr, "XTENSION")
```

which throw an error if the keyword is not present and perform type
checking and conversion if allowed.  It is also possible to supply a
default value to return if the keyword is not present:

```julia
bscale = OIFITS.get_float(hdr, "BSCALE", 1.0)
bzero = OIFITS.get_float(hdr, "BZERO", 0.0)
xtension = OIFITS.get_string(hdr, "XTENSION", "IMAGE")
```

The function:

```julia
OIFITS.get_hdu_type(hdr)
```

returns the HDU type as a Symbol, `:image_hdu` for an image, `:ascii_table` for
an ASCII table, `:binary_table` for a binary table, and `:unknown` otherwise.
The returned symbol should match the result of the low level method
`FITSIO.Libcfitsio.fits_get_hdu_type`.


For a FITS table, the function:

```julia
OIFITS.get_datablock_type(hdr)
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

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/OIFITS.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/OIFITS.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.org/emmt/OIFITS.jl.svg?branch=master
[travis-url]: https://travis-ci.org/emmt/OIFITS.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/OIFITS.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/OIFITS-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/OIFITS.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/OIFITS.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/OIFITS.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/OIFITS.jl?branch=master
