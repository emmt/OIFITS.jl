# v0.2.0 (unreleased)

## User interface

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
