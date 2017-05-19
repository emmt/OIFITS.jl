# Deprecations in OIFITS module.

import Base: @deprecate

# Deprecated in v0.2
@deprecate readtable(ff::FITSFile) read_table(ff::FITSFile)

# Deprecated in v0.3.2
@deprecate name2Symbol(name::AbstractString) symbolicname(name::AbstractString)
