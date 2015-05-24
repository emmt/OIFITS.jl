# Deprecations in OIFITS module.

import Base: @deprecate

# Deprecated in v0.2
@deprecate readtable(ff::FITSFile) read_table(ff::FITSFile)
