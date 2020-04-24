# Deprecations in OIFITS module.

import Base: @deprecate

# Deprecated in v0.5
@deprecate symbolicname(name::AbstractString) to_fieldname(name)
@deprecate get_dbname get_extname
@deprecate attach!(master::OIMaster, db::OIDataBlock) push!(master, db)
@deprecate get_real(hdr::FITSHeader, key::AbstractString) get_float(hdr, key)
@deprecate get_real(hdr::FITSHeader, key::AbstractString, def) get_float(hdr, key, def)
@deprecate fixname fix_name
@deprecate get_dbtype get_datablock_type

# Deprecated in v0.2
@deprecate readtable(ff::FITSFile) read_table(ff::FITSFile)

# Deprecated in v0.3.2
@deprecate name2Symbol(name::AbstractString) symbolicname(name::AbstractString)
