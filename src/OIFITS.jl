#
# OIFITS.jl --
#
# Support for OI-FITS (optical interferometry data format) in Julia.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015, Éric Thiébaut.
#
#------------------------------------------------------------------------------

module OIFITS

import Base: getindex, setindex!, haskey, keys, start, done, next

export oifits_get_time, oifits_get_mjd, oifits_get_int_time,
       oifits_get_sta_index, oifits_get_flag, oifits_get_visamp,
       oifits_get_visamperr, oifits_get_visphi, oifits_get_visphierr,
       oifits_get_vis2data, oifits_get_vis2err, oifits_get_t3amp,
       oifits_get_t3amperr, oifits_get_t3phi, oifits_get_t3phierr,
       oifits_get_ucoord, oifits_get_vcoord, oifits_get_u1coord,
       oifits_get_v1coord, oifits_get_u2coord, oifits_get_v2coord,
       oifits_get_date_obs, oifits_get_arrname, oifits_get_insname,
       oifits_get_revn, oifits_get_frame, oifits_get_arrayx,
       oifits_get_arrayy, oifits_get_arrayz, oifits_get_tel_name,
       oifits_get_sta_name, oifits_get_sta_index, oifits_get_diameter,
       oifits_get_staxyz, oifits_get_eff_wave, oifits_get_eff_band,
       oifits_get_target_id, oifits_get_target, oifits_get_raep0,
       oifits_get_decep0, oifits_get_equinox, oifits_get_ra_err,
       oifits_get_dec_err, oifits_get_sysvel, oifits_get_veltyp,
       oifits_get_veldef, oifits_get_pmra, oifits_get_pmdec,
       oifits_get_pmra_err, oifits_get_pmdec_err, oifits_get_parallax,
       oifits_get_para_err, oifits_get_spectyp

export oifits_new_target, oifits_new_array, oifits_new_wavelength,
       oifits_new_vis, oifits_new_vis2, oifits_new_t3

export oifits_read_header, oifits_get_hdutype, oifits_get_colnum,
       oifits_get_dbtype, oifits_get_value, oifits_get_comment,
       oifits_get_logical, oifits_get_integer,
       oifits_get_real, oifits_get_string,
       oifits_read_column, oifits_read_datablock

include("oidata.jl")
include("oifile.jl")
include("oiformat1.jl")
include("oipost.jl") # must be the last one

end # module

# Local Variables:
# mode: Julia
# tab-width: 8
# indent-tabs-mode: nil
# fill-column: 79
# coding: utf-8
# ispell-local-dictionary: "american"
# End:
