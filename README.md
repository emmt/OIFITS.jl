# OIFITS.jl
Support for OI-FITS (optical interferometry data format) in Julia.

## OI-FITS Summary

OI-FITS is a standard to store optical interferometry data as a collection
of data-blocks.  In the first version of the standard, the available
data-blocks are:

* `OI_TARGET` provides a list of observed targets;
* `OI_ARRAY` describes a given array of stations;
* `OI_WAVELENGTH` describes a given instrument (notably the effective
  wavelengths and bandwidths of its spectral channels);
* `OI_VIS` contains complex visibility data;
* `OI_VIS2` contains squared visibility (powerspectrum) data;
* `OI_T3` contains triple product (bispectrum) data.

These data-blocks, are stored as binary tables in a FITS data file.

The objective of the `OIFITS.jl` package is to provide support of
OI-FITS data in Julia language.


## Prerequisites

You need the [`FITSIO.jl`](https://github.com/JuliaAstro/FITSIO.jl) package
installed to use `OIFITS.jl`:
```julia
Pkg.add("FITSIO")
Pkg.update()
```


## Typical usage

Loading an OI-FITS data file:
```julia
using OIFITS
master = oifits_load("testdata.oifits")
```

To iterate through all data-blocks:
```julia
for db in master
    dbname = oifits_dbname(db)
    revn = oifits_get_revn(db)
    println("Data block is $dbname, revision $revn")
end
```

To iterate through a sub-set of the data-blocks (here the complex visibility
data, the powerspectrum data and the bispectrum data):
```julia
for db in oifits_select(master, "OI_VIS", "OI_VIS2", "OI_T3")
    dbname = oifits_dbname(db)
    n = length(oifits_get_time(db))
    println("Data block is $dbname, number of exposures is $n")
end
```


## Accessor functions

Any OI-FITS field (keyword/column) of a given data-block can be retrieved
via an accessor whose name has suffix `oifits_get_` followed by the name of
the field (in lower case letters and with all non-letter and all non-digit
letters replaced by the underscore character `'_'`).  A notable exception is
the revision number corresponding to the keyword "OI_REVN" which is
retrieved with the method `oifits_get_revn()`.  For instance:

```julia
oifits_get_revn(db)      # get the revison number of the format (OI_REVN)
oifits_get_eff_wave(db)  # get effective wavelengths (EFF_WAVE)
oifits_get_eff_band(db)  # get effective bandwidths (EFF_BAND)
oifits_get_ucoord(db)    # get the U coordinates of the data (UCOORD)
```
Of course, getting a given field must make sense.  For instance,
`oifits_get_eff_wave()` can be applied on any `OI_WAVELENGTH` data-blocks
but also on data-blocks which contains interferometric data such as
`OI_VIS`, `OI_VIS2`, `OI_T3`, *etc.* but not on other data-blocks like
`OI_TARGET`.


## Reading data

To load the contents of an OI-FITS file in memory, use:
```julia
master = oifits_load(filename)
```
where `filename` is the name of the file and the returned value, `master`,
contains all the OI-FITS data-blocks of the file.


## Constructors

It is possible to build OI-FITS data-blocks individually.  The general
syntax is:
```julia
oifits_new_XXX(KEY1=VAL1, KEY2=VAL2, ...)
```
where `XXX` is the type of the data-block and `KEYn=VALn` constructions
give the fields of the data-block and their values.  The names of the
fields follow the same convention as for the field accessors.

Available data-block constructors are:

* `oifits_new_target` => `OI_TARGET`
* `oifits_new_array` => `OI_ARRAY`
* `oifits_new_wavelength` => `OI_WAVELENGTH`
* `oifits_new_vis`  => `OI_VIS`
* `oifits_new_vis2` => `OI_VIS2`
* `oifits_new_t3`   => `OI_T3`

When defining a new data-block, all mandatory fields must be provided.
For instance, to create an `OI_WAVELENGTH` data-block:
```julia
µm = 1e-6  # all values are in SI units in OI-FITS
db = oifits_new_wavelength(insname="Amber",
                           eff_wave=[1.4µm,1.6µm,1.8µm],
                           eff_band=[0.2µm,0.2µm,0.2µm])
```
Note that the revision number (`revn=...`) can be omitted; by default, the
highest defined revision will be used.

A consistent set of OI-FITS data-blocks is made of: exactly one `OI_TARGET`
data-block, one or more `OI_WAVELENGTH` data-block, one or more `OI_ARRAY`
data-block and any number of data-blocks with interferometric data
(`OI_VIS`, `OI_VIS2` or `OI_T3`).  These data-blocks must be stored in a
container created by:
```julia
master = oifits_new_master()
```
Then, call:
```julia
oifits_attach(master, db)
```
to attach all data-block `db` to the OI-FITS container (in any order).
Finally, you must call:

```julia
oifits_update(master)
```
to update internal information such as links between data-blocks with
interferometric data and the related instrument (`OI_WAVELENGTH`
data-block) and array of stations (`OI_ARRAY` data-block).  If you do not
do that, then applying some accessors may not work, *e.g.*
`oifits_get_eff_wave()` on a data-blocks with interferometric data.

To read an OI-FITS data-block from the current HDU of a FITS file:
```julia
db = oifits_read_datablock(ff)
```
where `ff` is a `FITSFile` handle.  The result may be `nothing` if the
current HDU does not contain an OI-FITS data-block.


## Miscellaneous functions

OI-FITS implements some useful functions which can be used to deal with
FITS file (not just OI-FITS ones).  These functions could be part of `FITSIO`
package.


### Reading the header of a FITS HDU

The header of a FITS HDU can be read with the function:
```julia
hdr = oifits_read_header(ff::FITSFile)
```
which returns an indexable and iterable object, here `hdr`.  The keys of
`hdr` are the FITS keywords of the header.  For instance:
```julia
hdr = oifits_read_header(ff)
keys(hdr)          # yield an iterator on the keys of hdr
collect(keys(hdr)) # yield all the keys of hdr
haskey(hdr, key)   # check whether key is present
hdr[key]           # retrieve the contents associated with the key
```
The contents associated with a keyword are `(val,com)` tuples with `val`
the value of the keyword (in a suitable Julia type) and `com` the
associated comment.   An exemple to iterate over all keys is:
```julia
for (key, (value, comment)) in hdr
    println("$key = $value / $comment")
end
```
For commentary FITS keywords (`"HISTORY"` or `"COMMENT"`), the value and
the comment are identical (in the sense that they are references to the
same object) and are a vector of strings (one for each commentary card of
the given keyword).  Other keywords must be unique and thus have a scalar
value.

Retrieving only the value or the comment part is done with the
`oifits_get_value()` or the `oifits_get_comment()` methods which also
accept a default value provided if the keyword is not present:
```julia
val = oifits_get_value(hdr, key)
val = oifits_get_value(hdr, key, def)
com = oifits_get_comment(hdr, key)
com = oifits_get_comment(hdr, key, def)
```

In addition to the indexation, *e.g.,* `hdr[key]`, to the
`oifits_get_value()` or the `oifits_get_comment()` methods, the specific
value of a keyword can be retrieved by one of the following four different
methods:
```julia
oifits_get_logical(hdr, "SIMPLE")
oifits_get_integer(hdr, "BITPIX")
oifits_get_real(hdr, "BSCALE")
oifits_get_string(hdr, "XTENSION")
```
which throw an error if the keyword is not present and perform type
checking and conversion if allowed.  It is also possible to supply a
default value to return if the keyword is not present:
```julia
bscale = oifits_get_real(hdr, "BSCALE", 1.0)
bzero = oifits_get_real(hdr, "BZERO", 0.0)
xtension = oifits_get_string(hdr, "XTENSION", "IMAGE")
```

The function:
```julia
oifits_get_hdutype(hdr)
```
returns the HDU type as a symbol, `:image_hdu` for an image, `:ascii_table`
for an ASCII table, `:binary_table` for a binary table, and `:unknown`
otherwise.

For a FITS table, the function:
```julia
oifits_get_colnum(hdr, colname)
```
returns the index of the column whose name is `colname`; -1 is returned if
the column is not present or if the HDU is not a table.  Column names
should be given in uppercase letters and without trailing spaces (according
to FITS convention, letter case and trailing spaces are insignificant).

For a FITS table, the function:
```julia
oifits_get_dbtype(hdr)
```
returns the OI-FITS data-block type as a symbol like `:OI_TARGET`,
`:OI_WAVELENGTH`, *etc.*


### Reading a given column from a FITS table

The function
```julia
oifits_read_column(ff::FITSFile, colnum::Integer)
```
returns a Julia array with the contents of the `colnum`-th column of the
current HDU in FITS file handle `ff`.  The current HDU must be a FITS table
(an ASCII or a binary one).  The last dimension of the result corresponds
to the rows of the table.


### FITS and Julia types conversion

The functions `fits_datatype()` and `fits_bitpix()` deal with conversion
between CFITSIO type code or BITPIX value and actual Julia data types.
They can be used as follows (assuming `T` is a Julia data type, while
`code` and `bitpix` are integers):
```julia
fits_datatype(T) --------> code (e.g., TBYTE, TFLOAT, etc.)
fits_datatype(code) -----> T

fits_bitpix(T) ----------> bitpix (e.g., BYTE_IMG, FLOAT_IMG, etc.)
fits_datatype(bitpix) ---> T
```

The functions `fits_get_coltype()` and `fits_get_eqcoltype()` yield the
data type, repeat count and width in bytes of a given column, their
prototypes are:
```julia
(code, repcnt, width) = fits_get_coltype(ff::FITSFile, colnum::Integer)
(code, repcnt, width) = fits_get_eqcoltype(ff::FITSFile, colnum::Integer)
```
with `colnum` the column number, `code` the CFITSIO column type (call
`fits_datatype(code)` to convert it to a Julia type) of the elements in
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
