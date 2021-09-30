# OIFITS.jl

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |


The `OIFITS.jl` package provides support for OI-FITS data in Julia language.


## OI-FITS types

OI-FITS is a standard to store optical interferometry data as a collection of
data-blocks.  In the second revision of the standard (see [Ref. 1](#references)
and [Ref. 2](#references)), the following data-blocks are available:

* an `OI_TARGET` data-block provides a list of observed targets;
* each `OI_ARRAY` data-block describes a given array of stations;
* each `OI_WAVELENGTH` data-block describes a given instrument notably the effective
  wavelengths and bandwidths of its spectral channels;
* `OI_CORR` data-blocks store correlation data;
* `OI_VIS` data-blocks store complex visibility dat;
* `OI_VIS2` data-blocks store squared visibility (powerspectrum)
  data;
* `OI_T3` data-blocks store triple product (bispectrum) data;
* `OI_FLUX` data-blocks store spectral flux data;
* `OI_INSPOL` data-blocks contain instrumental polarization data.

These data-blocks, are stored as binary tables in a FITS data file.  The
support for FITS files is provided by the
[`FITSIO.jl`](https://github.com/JuliaAstro/FITSIO.jl) package.

The julia type of an OI-FITS data-block is named as the
corresponding OI-FITS extension. In addition to these types for individual OI-FITS data-blocks, the `OIFITS.jl`
package provides data-sets (of type `OIDataSet`) that contain several OI-FITS
data-blocks. Each data-set is an efficient representation of the contents of a
compliant OI-FITS file.


## Reading and writing OI-FITS files

Reading an OI-FITS data file yields a data-set and is done by:

```julia
using OIFITS
ds = read(OIDataSet, input)
```

where `input` it the name of the OI-FITS file or an instance of `FITSIO.FITS` which
represents an open FITS file.   The above `read` call is equivalent
to the shortcut:

```julia
ds = OIDataSet(input)
```

It is possible to merge the contents of several OI-FITS file, say `inp1`,
`inp2`, etc., by one of:

```julia
ds = read(OIDataSet, inp1, inp2, ...)
ds = OIDataSet(inp1, inp2, ...)
```

or to merge them into an existing data-set `ds`:

```julia
read!(ds, inp1, inp2, ...)
```

To write an OI-FITS file, just write the data-set:

```julia
write(filename, ds)
```

Overwriting is forbidden by default, but the keyword `overwrite=true` may be
specified to allow for silently overwriting an existing file.


## Accessing the contents of data-blocks and data-sets

The contents of OI-FITS data-blocks and data-sets may be accessed by the dot
notation but also by indexation.


### Contents of data-sets

The dot notation can be used on a *data-set* object, say `ds`, storing a
consistent set of OI-FITS data-blocks.  The following properties are available:

- `ds.target` is the `OI_TARGET` data-block of the OI-FITS structure.

- `ds.instr` is a list of `OI_WAVELENGTH` data-blocks indexed by a regular
  integer index or by the instrument name:

  ```julia
  ds.instr[i]       # yields the i-th OI_WAVELENGTH data-block
  ds.instr[insname] # yields the OI_WAVELENGTH data-block whose name
                    # matches insname
  ```

  Matching of names follows FITS conventions that case of letters and trailing
  spaces are ignored.  An error is thrown if the index (integer or name) is not
  valid.  The `get` method can be used to provide a default value, for
  instance:

  ```julia
  get(ds.instr, insname, nothing)
  ```

  would yield `nothing` is `insname` is not found in `ds.instr`.

- `ds.array` is a list of `OI_ARRAY` data-blocks indexed like `ds.instr` except
  that interferometric array names are assumed.

- `ds.correl` is a list of `OI_CORR` data-blocks indexed like `ds.instr` except
  that correlation data array names are assumed.

- `ds.vis` is a vector of `OI_VIS` data-blocks.

- `ds.vis2` is a vector of `OI_VIS2` data-blocks.

- `ds.t3` is a vector of `OI_T3` data-blocks.

- `ds.flux` is a vector of `OI_FLUX` data-blocks.

- `ds.inspol` is a vector of `OI_INSPOL` data-blocks.

Other fields of data-sets shall be considered as **private** and not accessed
directly.

Using the dot notaion, it is easy to access to the different data-blocks containing measurements.  For
instance:

```julia
for db in ds.vis2
    ...
end
```

is convenient to loop across all `OI_VIS2` instances stored by `ds`.


### Contents of data-blocks

The contents of a data-block, say `db`, may also be accessed by the dot
notation.  As a general rule, `db.key` or `db.col` yield the value of the
keyword `key` or the contents of the column `col` of the OI-FITS table
corresponding to the data-block `db`.  In order to follow Julia conventions and
to accommodate for some restrictions, `key` and `col` are the FITS keyword of
column name converted to lower case letters and with non-alphanumeric letters
replaced by underscores.  For instance `db.int_time` yields the values of the
column `INT_TIME`, that is the integration times of the measurements.  The
revision number corresponding to the keyword `OI_REVN` is however accessed as
`db.revn`, this is the only exception.  Other properties are also accessible
via this syntax:

- `db.extname` yields the OI-FITS name of the extension corresponding to the
  data-block `db` (for all data-block types);

- `db.array` yields the `OI_ARRAY` data-block associated with data-block `db`
  (only for `OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, and `OI_INSPOL`
  data-block);

- `db.instr` yields the `OI_WAVELENGTH` data-block associated with data-block
  `db` (only for `OI_VIS`, `OI_VIS2`, `OI_T3`, and `OI_FLUX` data-block);

- `db.correl` yields the `OI_CORR` data-block associated with data-block `db`
  (only for `OI_VIS`, `OI_VIS2`, `OI_T3`, and `OI_FLUX` data-block).

- `db.name` is an alias for `db.arrname` for `OI_ARRAY` instances, for
  `db.insname` for `OI_WAVELENGTH` instances, and for `db.corrname` for
  `OI_CORR` instances,

Of course, getting a given property must make sense.  For example,
`db.sta_name` is only possible for `OI_ARRAY` data-blocks but not for an
`OI_WAVELENGTH` data-block.  The dot notation can be however be chained and:

```julia
db.instr.eff_wave
```

can be used to access the effective wavelengths of the measurements in `db` via
the instrument associated to `db`.  Shortcuts are provided:

```julia
λ  = db.eff_wave # get effective wavelength
Δλ = db.eff_band # get effective bandwidth
```

for `OI_VIS`, `OI_VIS2`, `OI_T3`, and `OI_FLUX` data-blocks.

Some fields of a data-block `db` may however be undefined because:

- the field is not yet defined (the data-block is being constructed);

- the field is optional in the revision `db.revn` of the data-block;

- the field (for example `db.instr` for an `OI_VIS` data-block) involves links
  with other data-blocks (the *dependencies)* and these links are only defined
  when a data-block is pushed into a data-set.


## Building of data-sets

Reading an OI-FITS file is the easiest way to define a data-set but a new
OI-FITS data-set may be built by creating an empty data-set with `OIDataSet()`,
and then pushing OI-FITS data-blocks **in order** vith `push!(...)`.  Indeed,
in order to ensure the consistency of a data-set, it is required to push into a
data-set the dependencies (`OI_TARGET`, `OI_ARRAY`, `OI_WAVELENGTH`, and
`OI_CORR` data-blocks) **before** the data-blocks containing measurements
(`OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, and `OI_INSPOL`) that may refer to
them.

For example, building a new data-set, say `ds`, looks like:

```julia
ds = OIDataSet() # create empty data-set
push!(ds, arr)   # push OI_ARRAY data-block(s)
push!(ds, ins)   # push OI_WAVELENGTH data-block(s)
push!(ds, cor)   # push OI_CORR data-block(s)
push!(ds, tgt)   # push OI_TARGET data-block
push!(ds, db1)   # push data
push!(ds, db2)   # push more data
...
```

with the dependencies:

- `arr` an `OI_ARRAY` instance defining the interferometric array (several such
   instances may be pushed),

- `ins` an `OI_WAVELENGTH` instance defining the instrument (several such
   instances can be pushed),

- `cor` an `OI_COORREL` instance defining the correlations (zero or any number
   of such instances can be pushed),

- `tgt` an `OI_TARGET` instance defining the list of observed targets (at least
  one such instance is required, if more such instances are pushed in the same
  data-set, they are merged in a single one);

and `db1`, `db2`, ... instances of `OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, or
`OI_INSPOL` that provide measurements.

You may push all data-blocks in a single `push!` call:

```julia
ds = push!(OIDataSet(), arr, ins, cor, tgt, d1, db2, ...)
```

and the following shortcut is implemented:

```julia
ds = OIDataSet(arr, ins, cor, tgt, d1, db2, ...)
```

These two are equivalent to the long example, but remember that pushing
data-blocks in order (i.e., dependencies before they may be referenced) is
required to have a consistent data-set.  Apart from this constraint,
dependencies may be pushed in any order **before** the data-blocks with
measurements and data-blocks with measurements can be be pushed in any order
**after** dependencies.

As a benefit of the constraint of pushing data-blocks in order, data-blocks
with dependencies are automatically linked to these dependencies when pushed on
the data-set (which implies that the dependencies already exist in the
data-set).  This allows for syntaxic sugar like:

```julia
ds.vis2[i].eff_wave # the wavelengths of the i-th OI_VIS2 data-block in ds
ds.t3[i].array      # the interferometric array for the i-th OI_T3 data-block in ds
ds.vis[i].instr     # the instrument used for the i-th OI_VIS data-block in ds
```

Without linked dependencies, the first above example would require to (1) find
in the data-set `ds` the `OI_WAVELENGTH` instance, say `ins`, whose name is
matching `ds.vi2[i].insname` and (2) extract the field `eff_wave` of `ins`.
The latter setp is as simple as `ins.eff_wave` but the former one has some
overheads and scales as `O(n)` with `n` the number of `OI_WAVELENGTH` instances
in the data-set.

Since an OI-FITS data-set has a single list of targets (an `OI_TARGET` instance
accessible via `ds.target`), a mean to merge list of targets had to de defined.
The adopted rule is pretty simple:

> The `target_id` field of any data-block that is part of a data-set correspond
> to the index of the target entry in the list of targets stored by the
> data-set.

As a consequence, whenever a data-block is pushed into a data-set, the target
identifiers of the data-block have to be rewritten according to this rule.  Of
course this does not apply for data-blocks with no `target_id` field such as
`OI_ARRAY`, `OI_WAVELENGTH`, and `OI_CORR`.

To summarize, here is what happens under the hood when a data-block `db` is
pushed into a data-set `ds`:

- When an `OI_ARRAY`, `OI_WAVELENGTH`, or `OI_CORR` instance `db` is pushed in
  a data-set `ds`, it is appended to the corresponding list (`ds.array`,
  `ds.instr`, or `ds.correl`) unless this list already has an entry with a name
  matching `db.name`.  In this latter case, an assertion exception is thrown if
  the two data-blocks whose names are matching do not have the same contents
  (otherwise inconsistent data-sets could be built).

- When an `OI_TARGET` instance is pushed in a data-set, the new targets
  (according to their names) are appended to the list of targets in the
  data-set and their identifiers set to their index in this list.  This also
  re-initialize an internal dictionary used to perform the conversion from all
  the target identifiers of the `OI_TARGET` instance that has been pushed to
  the target identifiers in the data-set.  Until it is reinitialized (by
  pushing another `OI_TARGET` instance), this mapping is used to rewrite the
  target identifiers of subsequent data-blocks pushed in the data-set.

- When an `OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, or `OI_INSPOL` instance
  `db` a data-set `ds`, it is appended to the corresponding list (`ds.vis`,
  `ds.vis2`, `db.t3`, `db.flux`, or `ds.inspol`), after it has been linked to
  its dependencies (`OI_ARRAY`, `OI_WAVELENGTH`, etc., which must already exist
  in the data-set), and its target identifiers are rewritten according to the
  mapping defined by the last `OI_TARGET` instance that has been pushed to the
  data-set.  Rewriting of the target identifiers may be avoided by using the
  keyword `rewrite_target_id=false`, this assumes that the target identifiers
  in the newly pushed data-blocks are equal to the index in the list of targets
  `ds.target`.

Pushing a data-block in a data-set does check the consistency of the
data-block.  This is to allow for building the data-blocks step by step so that
they not need to be consistent at all times (just when pushed into a data-set).

Pushing a data-block in a data-set left the data-block unchanged.  A swallow
copy of it is added to the data-blocks stored by the data-set.  Most members of
the pushed data-blocks are shared by the one stored by the data-set whith the
notable exception of the target identifiers which are rewritten and the links
to the dependencies which are updated.

While it sounds complicated, the default rule of rewriting the target
identifiers just amounts to assuming that the target identifiers of `OI_VIS`,
`OI_VIS2`, `OI_T3`, `OI_FLUX`, or `OI_INSPOL` instances pushed in a data-set
refer to the last `OI_TARGET` instance pushed on the same data-set.

Pushing several groups of data-blocks, each group making a consistent data-set,
in the same data-set is easy.  Typically:

```julia
# Push dependencies for group #1.
push!(ds, group1_arr)
push!(ds, group1_ins)
push!(ds, group1_cor)
# Push targets for group #1 (reinitializing target_id mapping).
push!(ds, group1_tgt)
# Push data for group #1 (using current target_id mapping).
push!(ds, group1_db1)
push!(ds, group1_db2)
...
# Push dependencies for group #2.
push!(ds, group2_arr)
push!(ds, group2_ins)
push!(ds, group2_cor)
# Push targets for group #2 (reinitializing target_id mapping).
push!(ds, group2_tgt)
# Push data for group #2 (using current target_id mapping).
push!(ds, group2_db1)
push!(ds, group2_db2)
...
```

Since they are referenced by their names, it is not needed to push `OI_ARRAY`,
`OI_WAVELENGTH`, and `OI_COORREL` dependencies if they already exist in the
data-set (according to their name).  It is however mandatory to push an
`OI_TARGET` instance with all targets and their identifiers as assumed by the
subsequent data-blocks.


# Merging of data-sets

Two OI-FITS data-sets (or more), say `A` and `B`, can be consistently merged
together by:

```julia
C = merge(A, B)
```

As much as possible, the resulting data-set `C` will share its contents with
`A` and/or `B` but without affecting `A` and `B` which are guaranteed to remain
unchanged.  As for pushing data-blocks, the target identifiers (the `target_id`
field) may be rewritten in the result.

Merging of data-sets assumes that the two merged data-sets are *consistent* and
*compatible*.  Here *compatible* means that targets and dependencies with
matching names must have the same contents.  This is checked during the merge
operation.

It is also allowed to merge several data-sets and/or merge data-sets *in-place*:

```julia
ds = merge(ds1, ds2, ds3, ...)  # merge ds1, ds2, ... in new data-set ds
merge!(ds, ds1, ds2, ds3, ...)  # merge ds1, ds2, ... in existing data-set ds
```

Note that `merge!(ds,...)` yields the destination `ds`.

Also note that, after merging, the internal dictionary used for rewriting
target identifiers is left with the mapping built from the targets of the last
merged data-set.


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
