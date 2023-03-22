# OIFITS.jl

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |


The `OIFITS.jl` package provides support for OI-FITS data in Julia language.


## OI-FITS types

OI-FITS is a standard to store optical interferometry data as a collection of
data-blocks. In the second revision of the standard (see [Ref. 1](#references)
and [Ref. 2](#references)), an OI-FITS file may contain the following
data-blocks:

* an `OI_TARGET` data-block stores a list of observed targets;
* each `OI_ARRAY` data-block describes a given array of telescope stations;
* each `OI_WAVELENGTH` data-block describes a given instrument notably the
  effective wavelengths and bandwidths of its spectral channels;
* `OI_CORR` data-blocks store correlation data;
* `OI_VIS` data-blocks store complex visibility data;
* `OI_VIS2` data-blocks store squared visibility (powerspectrum) data;
* `OI_T3` data-blocks store triple product (bispectrum) data;
* `OI_FLUX` data-blocks store spectral flux data;
* `OI_INSPOL` data-blocks store instrumental polarization data.

These data-blocks are stored as binary tables in a FITS data file. The support
for FITS files is provided by the
[`FITSIO.jl`](https://github.com/JuliaAstro/FITSIO.jl) package.

The Julia type of an OI-FITS data-block is named as the corresponding OI-FITS
extension. In addition to these types for individual OI-FITS data-blocks, the
`OIFITS.jl` package provides data-sets (of type `OIDataSet`) that contain
several OI-FITS data-blocks. Each data-set is an efficient representation of
the contents of a compliant OI-FITS file.


## Reading and writing OI-FITS files

### Reading and writing OI-FITS data-sets

Reading an OI-FITS data file in Julia yields a data-set and is done by:

```julia
using OIFITS
ds = read(OIDataSet, input)
```

where `input` it the name of the OI-FITS file or an instance of `FITSIO.FITS`
which represents an open FITS file. The above `read` call is equivalent to the
shortcut:

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

Creating an OI-FITS file is as simple as writing the data-set `ds`:

```julia
write(filename, ds)
```

Overwriting is forbidden by default, but the keyword `overwrite=true` may be
specified to allow for silently overwriting an existing file.


### Reading individual OI-FITS data-blocks

It may be useful to read individual OI-FITS data-blocks, to debug or to fix the
contents of a non-compliant OI-FITS file. To that end, you must open the FITS
file and can then read a given HDU as an OI-FITS data-block:

```julia
using FITSIO, OIFITS
f = FITS(filename, "r")     # open FITS file for reading
tgt = OI_TARGET(f[i])       # read OI_TARGET extension in i-th HDU
tgt = read(OI_TARGET, f[i]) # idem
db = OI_VIS2(f[j])          # read OI_VIS2 extension in j-th HDU
db = read(OI_VIS2, f[j])    # idem
...
```

any OI-FITS data-block type can be used in that way.  If the type of the `i`-th
extension is not known, `OIDataBlock` can be used instead but the result is not
type-stable:

```julia
db = OIDataBlock(f[i])       # read OI-FITS extension extension in i-th HDU
db = read(OIDataBlock, f[i]) # idem
```

Writing individual OI-FITS data-blocks is also possible:

```julia
using FITSIO, OIFITS
f = FITS(filename, "w") # open FITS file for writing
write(f, db)            # write db in the next HDU of f
```

To fix a non-compliant OI-FITS file (usually duplicate target or instrument
names), you can read all the data-blocks, fix those which are wrong and push
them in **order** in an `OIDataSet` to have a consistent data-set which you can
then directly use or write in an OI-FITS file for later. Thanks to the
automatic rewriting of target identifiers and of the fact that targets (and
other dependencies) are identified by their name and consistently merged, it is
possible to push an `OI_TARGET` with multiply defined identical targets (apart
maybe their identifiers).


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
  ds.instr[insname] # yields the OI_WAVELENGTH data-block whose name matches insname
  ```

  Matching of names follows FITS conventions that case of letters and trailing
  spaces are ignored. An exception is thrown if the index (integer or name) is
  not valid. The `get` method can be used to provide a default value, for
  example:

  ```julia
  get(ds.instr, insname, nothing)
  ```

  would yield `nothing` if `insname` is not found in `ds.instr` instead of
  throwing an exception.

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

Using the dot notation, it is easy to access the different data-blocks
containing measurements. For instance:

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
to accommodate for a number of restrictions, `key` or `col` are the FITS
keyword or column name converted to lower case letters and with
non-alphanumeric letters replaced by underscores.  For instance `db.date_obs`
yields the value of the keyword `DATE-OBS`, that is the UTC start date of
observations.  The revision number corresponding to the keyword `OI_REVN` is
however accessed as `db.revn`, this is the only exception.  Other properties
are also accessible via this syntax:

- `db.extname` yields the OI-FITS name of the extension corresponding to the
  data-block `db` (for all data-block types);

- `db.array` yields the `OI_ARRAY` data-block associated with data-block `db`
  (only for `OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, and `OI_INSPOL`
  data-block).  Beware that the association with an `OI_ARRAY` is optional, so
  `db.array` may be actually undefined; this can be checked by
  `isdefined(db,:array)`.

- `db.instr` yields the `OI_WAVELENGTH` data-block associated with data-block
  `db` (only for `OI_VIS`, `OI_VIS2`, `OI_T3`, and `OI_FLUX` data-block).

- `db.correl` yields the `OI_CORR` data-block associated with data-block `db`
  (only for `OI_VIS`, `OI_VIS2`, `OI_T3`, and `OI_FLUX` data-block).

- `db.name` is an alias for `db.arrname` for `OI_ARRAY` instances, for
  `db.insname` for `OI_WAVELENGTH` instances, and for `db.corrname` for
  `OI_CORR` instances.

Of course, getting a given property must make sense.  For example,
`db.sta_name` is only possible for an `OI_ARRAY` data-block but not for an
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

for `OI_WAVELENGTH` data-blocks but also for `OI_VIS`, `OI_VIS2`, `OI_T3`, and
`OI_FLUX` data-blocks.

Some fields of a data-block `db` may however be undefined because:

- the field is not yet defined (the data-block is being constructed);

- the field is optional in the revision `db.revn` of the data-block;

- the field (for example `db.instr` for an `OI_VIS` data-block) involves links
  with other data-blocks (the *dependencies)* and these links are only defined
  when a data-block is part of a data-set (see [Building of
  data-sets](#building-of-data-sets) below).


### `OI_TARGET` data-blocks

For efficiency, instances of `OI_TARGET` data-blocks do not follow the same
rules as other types of OI-FITS data-blocks whose properties are the columns of
the corresponding OI-FITS table: in an `OI_TARGET` instance, all parameters
describing a target are repesented by an `OITargetEntry` structure and all
targets are stored as a vector of `OITargetEntry`.  An `OI_TARGET` instance,
say `db`, has the 3 following properties:

```julia
db.extname # yields "OI_TARGET"
db.list    # yields a vector of OITargetEntry instances
db.revn    # yields the revision number
```

The list of targets `db.list` can be indexed by an integer (as any Julia
vector) or by the target name (case of letters and trailing spaces are
irrelevant).

As an `OI_TARGET` data-blocks is essentially a vector of target entries, it can
be used as an iterable and it can indexed by an integer index or by a target
name:

```julia
length(db) # the number of targets, shortcut for `length(db.list)`
db[i]      # the i-th target, shortcut for `db.list[i]`
db[key]    # the target whose name matches string `key`, shortcut for `db.list[key]`
```

Standard methods `get` and `haskey`, applied to `db.list` or directly to `db`,
work as expected and according to the type (integer or string) of the key.  For
the `keys` method, the default is to return an iterator over the target names,
but the type of the expected keys can be specified:

```julia
get(db,key,def)   # yields `db[key]` or `def` if `key` not found
keys(db)          # iterator over target names
keys(String, db)  # idem
keys(Integer, db) # iterator over target indices
keys(Int, db)     # idem
```

The method `OIFITS.get_column` is a helper to recover a single target field as
a vector:

```julia
OIFITS.get_column([T,] db, col)
```

yields the column `col` of an OI-FITS data-block `db`.  Column is identified by
`col` which is either `sym` or `Val(sym)` where `sym` is the symbolic name of
the corresponding field in `OITargetEntry`.  Optional argument `T` is to
specify the element type of the returned array.

To build an `OI_TARGET` instance, you may provide the list of targets and the
revision number:

```julia
OI_TARGET(lst=OITargetEntry[]; revn=0)
```

yields an `OI_TARGET` data-block.  Optional argument `lst` is a vector of
`OITargetEntry` specifying the targets (none by default).  Keyword `revn`
specifies the revision number.

A target entry may be constructed by specifying all its fields (there are many)
by keywords, all of which but `category` are mandatory:

```julia
x = OITargetEntry(;
        target_id ::Integer,
        target    ::AbstractString,
        raep0     ::AbstractFloat,
        decep0    ::AbstractFloat,
        equinox   ::AbstractFloat,
        ra_err    ::AbstractFloat,
        dec_err   ::AbstractFloat,
        sysvel    ::AbstractFloat,
        veltyp    ::AbstractString,
        veldef    ::AbstractString,
        pmra      ::AbstractFloat,
        pmdec     ::AbstractFloat,
        pmra_err  ::AbstractFloat,
        pmdec_err ::AbstractFloat,
        parallax  ::AbstractFloat,
        para_err  ::AbstractFloat,
        spectyp   ::AbstractString,
        category  ::AbstractString = "")
```

It is also possible to specify another target entry, say `ref`, which is used
as a template: any unspecified keyword is assume to have the same value as in
`ref`:

```julia
x = OITargetEntry(ref;
        target_id = ref.target_id,
        target    = ref.target,
        ...)
```

Note that, when an `OI_TARGET` instance is pushed in a data-set, target
identifiers (field `target_id`) are automatically rewritten to be identical to
the index in the list of targets of the data-set.


## Building of data-sets

### Pushing data-blocks to data-sets

Reading an OI-FITS file is the easiest way to define a data-set but a new
OI-FITS data-set may be built by creating an empty data-set with `OIDataSet()`,
and then pushing OI-FITS data-blocks **in order** with `push!(...)`. Indeed, in
order to ensure the consistency of a data-set, it is required to push the
dependencies (`OI_TARGET`, `OI_ARRAY`, `OI_WAVELENGTH`, and `OI_CORR`
data-blocks) **before** the data-blocks containing measurements (`OI_VIS`,
`OI_VIS2`, `OI_T3`, `OI_FLUX`, and `OI_INSPOL`) that may refer to them.

For example, building a new data-set, say `ds`, looks like:

```julia
ds = OIDataSet() # create empty data-set
push!(ds, arr)   # push OI_ARRAY data-block(s)
push!(ds, ins)   # push OI_WAVELENGTH data-block(s)
push!(ds, cor)   # push OI_CORR data-block(s)
push!(ds, tgt)   # push OI_TARGET data-block
push!(ds, db1)   # push data
push!(ds, db2)   # push more data
push!(ds, db3)   # push even more data
...
```

with the dependencies:

- `arr` an `OI_ARRAY` instance defining the interferometric array (zero or any
   number of such instances may be pushed),

- `ins` an `OI_WAVELENGTH` instance defining the instrument (several such
   instances can be pushed),

- `cor` an `OI_COORREL` instance defining the correlations (zero or any number
   of such instances can be pushed),

- `tgt` an `OI_TARGET` instance defining the list of observed targets (at least
  one such instance is required, if more such instances are pushed in the same
  data-set, they are merged in a single one);

and where `db1`, `db2`, `db3`, etc., are instances of `OI_VIS`, `OI_VIS2`,
`OI_T3`, `OI_FLUX`, or `OI_INSPOL` that provide measurements.

You may push all data-blocks in a single `push!` call:

```julia
ds = push!(OIDataSet(), arr, ins, cor, tgt, d1, db2, ...)
```

and the following shortcut is implemented:

```julia
ds = OIDataSet(arr, ins, cor, tgt, d1, db2, ...)
```

These two are equivalent to the multi-line example above, but remember that
pushing data-blocks in order (i.e., dependencies before they may be referenced)
is required to have a consistent data-set. Apart from this constraint,
dependencies may be pushed in any order **before** the data-blocks with
measurements and data-blocks with measurements can be pushed in any order
**after** dependencies.

As a benefit of the constraint of pushing data-blocks in order, data-blocks
with dependencies are automatically linked to these dependencies when pushed on
the data-set (which implies that the dependencies already exist in the
data-set). This allows for syntactic sugar like:

```julia
ds.vis2[i].eff_wave # the wavelengths of the i-th OI_VIS2 data-block in ds
ds.t3[i].array      # the interferometric array for the i-th OI_T3 data-block in ds
ds.vis[i].instr     # the instrument used for the i-th OI_VIS data-block in ds
```

Without linked dependencies, the first above example would require to (1) find
in the data-set `ds` the `OI_WAVELENGTH` instance, say `instr`, whose name is
matching `ds.vi2[i].insname` and (2) extract the field `eff_wave` of `instr`.
The latter step is as simple as `instr.eff_wave` but the former one has some
overheads and scales as `O(n)` with `n` the number of `OI_WAVELENGTH` instances
in the data-set.

Since an OI-FITS data-set has a single list of targets (an `OI_TARGET` instance
accessible via `ds.target`), a mean to merge list of targets had to be defined.
The adopted rule is pretty simple:

> The `target_id` field of any data-block that is part of a data-set
> corresponds to the index of the target entry in the list of targets stored by
> the data-set.

As a consequence, whenever a data-block is pushed into a data-set, the target
identifiers of the data-block have to be rewritten according to this rule. Of
course, this does not apply for data-blocks with no `target_id` field such as
`OI_ARRAY`, `OI_WAVELENGTH`, and `OI_CORR`.

To summarize, here is what happens under the hood when a data-block `db` is
pushed into a data-set `ds`:

- When an `OI_ARRAY`, `OI_WAVELENGTH`, or `OI_CORR` instance `db` is pushed in
  a data-set `ds`, it is appended to the corresponding list (`ds.array`,
  `ds.instr`, or `ds.correl`) unless this list already has an entry with a name
  matching `db.name`. In this latter case, nothing is done unless that an
  assertion exception is thrown if the two data-blocks whose names are matching
  do not have the same contents (to prevent building inconsistent data-sets).

- When an `OI_TARGET` instance is pushed in a data-set, the new targets
  (according to their names) are appended to the list of targets in the
  data-set and their identifiers set to their index in this list. This also
  re-initializes an internal dictionary used to perform the conversion from all
  the target identifiers of the `OI_TARGET` instance that has been pushed to
  the target identifiers in the data-set. Until it is reinitialized (by pushing
  another `OI_TARGET` instance), this mapping is used to rewrite the target
  identifiers of subsequent data-blocks pushed in the data-set.

- When an `OI_VIS`, `OI_VIS2`, `OI_T3`, `OI_FLUX`, or `OI_INSPOL` instance `db`
  is pushed in a data-set `ds`, it is appended to the corresponding list
  (`ds.vis`, `ds.vis2`, `db.t3`, `db.flux`, or `ds.inspol`), after it has been
  linked to its dependencies (`OI_ARRAY`, `OI_WAVELENGTH`, etc., which must
  already exist in the data-set), and its target identifiers have been
  rewritten according to the mapping defined by the last `OI_TARGET` instance
  previously pushed to the data-set. Rewriting of the target identifiers may be
  avoided by using the keyword `rewrite_target_id=false`, this assumes that the
  target identifiers in the pushed data-block are already set according to the
  index in the list of targets `ds.target`.

Pushing a data-block in a data-set does check the consistency of the
data-block. This is to allow for building the data-blocks step by step so that
they not need to be consistent at all times (just when pushed into a data-set).

Pushing a data-block in a data-set lefts the data-block unchanged. A swallow
copy of it is added to the data-blocks stored by the data-set. Most members of
the pushed data-blocks are shared by the one stored by the data-set with the
notable exception of the target identifiers which are rewritten and the links
to the dependencies which are updated.

While it sounds complicated, the default rule of rewriting the target
identifiers just amounts to assuming that the target identifiers of `OI_VIS`,
`OI_VIS2`, `OI_T3`, `OI_FLUX`, or `OI_INSPOL` instances pushed in a data-set
refer to the last `OI_TARGET` instance previously pushed on the same data-set.

Pushing several groups of data-blocks, each group making a consistent data-set,
in the same data-set is easy. Typically:

```julia
# First push dependencies for group 1.
push!(ds, group1_arr) # push OI_ARRAY
push!(ds, group1_ins) # push OI_INS
push!(ds, group1_cor) # push OI_CORR
push!(ds, group1_tgt) # push OI_TARGET (reinitializing target_id mapping)
# Then push data for group 1 (using current target_id mapping).
push!(ds, group1_db1)
push!(ds, group1_db2)
...
# First push dependencies for group 2.
push!(ds, group2_arr) # push OI_ARRAY
push!(ds, group2_ins) # push OI_INS
push!(ds, group2_cor) # push OI_CORR
push!(ds, group2_tgt) # push OI_TARGET (reinitializing target_id mapping)
# Then push data for group 2 (using current target_id mapping).
push!(ds, group2_db1)
push!(ds, group2_db2)
...
```

Since they are referenced by their names, it is not necessary to push
`OI_ARRAY`, `OI_WAVELENGTH`, and `OI_COORREL` dependencies if they already
exist in the data-set (according to their name), but it doesn't hurt. It is
however mandatory to push an `OI_TARGET` instance with all targets and their
identifiers as assumed by the subsequent data-blocks.


## Merging data-sets

Two OI-FITS data-sets (or more), say `A` and `B`, can be consistently merged
together by:

```julia
C = merge(A, B)
```

As much as possible, the resulting data-set `C` will share its contents with
`A` and/or `B` but without affecting `A` and `B` which are guaranteed to remain
unchanged. As for pushing data-blocks, the target identifiers (the `target_id`
field) may be rewritten in the result.

Merging of data-sets assumes that the two merged data-sets are *consistent* and
*compatible*. Here *compatible* means that targets and dependencies with
matching names must have the same contents. This is checked during the merge
operation.

It is also allowed to merge several data-sets and/or merge data-sets
*in-place*:

```julia
ds = merge(ds1, ds2, ds3, ...) # merge ds1, ds2, ... in new data-set ds
merge!(ds, ds1, ds2, ds3, ...) # merge ds1, ds2, ... in existing data-set ds
```

Note that `merge!(ds,...)` yields the destination `ds`.

Also note that, after merging, the internal dictionary used for rewriting
target identifiers is left with the mapping built from the targets of the last
merged data-set.


## Credits

The development of this package has received funding from the European
Community's Seventh Framework Programme (FP7/2013-2016) under Grant Agreement
312430 (OPTICON).


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

[travis-img]: https://travis-ci.com/emmt/OIFITS.jl.svg?branch=master
[travis-url]: https://travis-ci.com/emmt/OIFITS.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/OIFITS.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/OIFITS-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/OIFITS.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/OIFITS.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/OIFITS.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/OIFITS.jl?branch=master
