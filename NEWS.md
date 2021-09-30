# User visible changes in Julia interface to OI-FITS


## Things to be done

`OIFITS` package is *work-in-progress*, it can currently read any compliant
OI-FITS file and provides easy and fast access to the contents of a set of OI
data.  In a near future more capabilities will be implemented to meet the
requirements of data processing softwares:

- Provide filter method to extract sub-sets of OI-FITS data: to select a chosen
  wavelength range, a given target, etc.

- Enforce column types in written OI-FITS file.

- Automatically rewrite `sta_index`, and indices in correlation matrices so
  that they match the indices in the related arrays.

- Automatically deal with revision number of `OI_TARGET` instances.

- More fully check data-blocks on `push!`.

- If a data-block has linked dependencies, push them before pushing the
  data-block itself.

- Hide private fields.


## Versions 1.0

This version introduces major changes.

- Data-blocks are now structured objects and no longer dictionaries.  As a
  result, type-stability of structures and methods has been much improved.
  Accessing the contents of data-blocks is considerably faster.

- An OI-FITS file can be read as a data-set (that is an instance of
  `OIDataSet`), and a data-set can be written to an OI-FITS file.

- The `obj.key` syntax is encouraged for any OI-FITS object `obj`.  Accessors
  functions are no longer needed (and have been discarded).  Provided
  dependencies have been correctly set, shortcuts are provided so that
  `obj.eff_wave` (or `obj.eff_band`) can be used instead of
  `obj.instr.eff_wave` (or `obj.instr.eff_band`) for any object `obj` storing
  OI-FITS data.

- Building a data-set (that is an instance of `OIDataSet`) by calling `push!`
  automatically takes care of linking dependencies and of verifying the
  consistencey of the data-set.  A *copy-on-write* policy is applied to avoid
  side-effects when pushing data-blocks.  Merging data-sets is easily done by
  `merge` or `merge!`.

- Calling `using OIFITS` only exports OI-FITS types (all prefixed with `OI*`)
  and no methods other than type constructors.

- `Base.copy` has been extended and replaces `OIFITS.clone`, it yields a
  swallow copy of a data-block.

- `OIMaster` has been renamed as `OIDataSet` and all types renamed to follow
  the names of the extensions in the OI-FITS specifications.  Thus `OITarget`,
  `OIArray`, `OIWavelength`, `OICorr` (or `OICorrelation`), `OIVis`, `OIVis`,
  `OIT`, `OIFlux`, and `OIInsPol` (or `OIPolarization`) renamed as `OI_TARGET`,
  `OI_ARRAY`, `OI_WAVELENGTH`, `OI_CORR`, `OI_VIS`, `OI_VIS2`, `OI_T3`,
  `OI_FLUX`, and `OI_INSPOL` respectively.

- Macros `@header` and `@column` are provided to define OI-FITS formats with a
  syntax very close to the tables in OI-FITS specifications.

- The package no longer hacks `FITSIO` and `CFITSIO` packages to handle FITS
  files.  As a result, `OIFITS` should be much less sensitive to the evolution
  of these dependencies.


## Versions 0.4

 - All `oifits_` prefixes have been removed and no methods are exported.

 - In most cases, `oifits_*` prefix has to be replaced by `OIFITS.`; for
   instance, use `OIFITS.load` instead of `oifits_load`.  The same for
   `OIFITS.read_datablock`, `OIFITS.read_column`, `OIFITS.new_target`,
   `OIFITS.new_array`, `OIFITS.new_wavelength`, `OIFITS.new_vis`,
   `OIFITS.new_vis2`, `OIFITS.new_t3`, `OIFITS.new_master`,
   `OIFITS.attach!`, `OIFITS.update`, `OIFITS.select`,
   `OIFITS.get_hdutype`, `OIFITS.get_dbtype`, `OIFITS.get_time`,
   `OIFITS.get_mjd`, `OIFITS.get_int_time`, `OIFITS.get_sta_index`,
   `OIFITS.get_flag`, `OIFITS.get_visamp`, `OIFITS.get_visamperr`,
   `OIFITS.get_visphi`, `OIFITS.get_visphierr`, `OIFITS.get_vis2data`,
   `OIFITS.get_vis2err`, `OIFITS.get_t3amp`, `OIFITS.get_t3amperr`,
   `OIFITS.get_t3phi`, `OIFITS.get_t3phierr`, `OIFITS.get_ucoord`,
   `OIFITS.get_vcoord`, `OIFITS.get_u1coord`, `OIFITS.get_v1coord`,
   `OIFITS.get_u2coord`, `OIFITS.get_v2coord`, `OIFITS.get_date_obs`,
   `OIFITS.get_arrname`, `OIFITS.get_insname`, `OIFITS.get_revn`,
   `OIFITS.get_frame`, `OIFITS.get_arrayx`, `OIFITS.get_arrayy`,
   `OIFITS.get_arrayz`, `OIFITS.get_tel_name`, `OIFITS.get_sta_name`,
   `OIFITS.get_sta_index`, `OIFITS.get_diameter`, `OIFITS.get_staxyz`,
   `OIFITS.get_eff_wave`, `OIFITS.get_eff_band`, `OIFITS.get_target_id`,
   `OIFITS.get_target`, `OIFITS.get_raep0`, `OIFITS.get_decep0`,
   `OIFITS.get_equinox`, `OIFITS.get_ra_err`, `OIFITS.get_dec_err`,
   `OIFITS.get_sysvel`, `OIFITS.get_veltyp`, `OIFITS.get_veldef`,
   `OIFITS.get_pmra`, `OIFITS.get_pmdec`, `OIFITS.get_pmra_err`,
   `OIFITS.get_pmdec_err`, `OIFITS.get_parallax`, `OIFITS.get_para_err`,
   `OIFITS.get_spectyp`.

 - Removed methods: `oifits_get_colnum`.

 - Method `oifits_read_header` replaced by `readheader` in FITSIO package.

 - Method `oifits_dbname` available as `OIFITS.get_dbname`.


# Internals

 - Lots of methods moved to [FITSIO](https://github.com/JuliaAstro/FITSIO.jl)

 - No longer exported methods: `oifits_get_hdutype` (available as `OIFITS.get_hdutype`).
