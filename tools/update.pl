#! /usr/bin/perl -pi.orig
#
# Notes:
#  - use flags '-pi.orig' to backup original files with extension .orig
#  - use flags '-pi' to just do in-place replacement (unsafe)
#
s/\bOIData\b/OIDataSet/g;
s/\bOITarget\b/OI_TARGET/g;
s/\bOIArray\b/OI_ARRAY/g;
s/\bOIWavelength\b/OI_WAVELENGTH/g;
s/\bOICorr\b/OI_CORR/g;
s/\bOIVis\b/OI_VIS/g;
s/\bOIVis2\b/OI_VIS2/g;
s/\bOIT3\b/OI_T3/g;
s/\bOIFlux\b/OI_FLUX/g;
s/\bOIInsPol\b/OI_INSPOL/g;
