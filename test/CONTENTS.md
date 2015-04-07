## Contents

* Files `contest-2004*.oifits` are those for the 2004 Opt/IR imaging Beauty
  Contest.

* Files `contest-2008*.oifits` are those for the 2008 Opt/IR imaging Beauty
  Contest.  There are two test sources, a model representing an AGB star
  (`contest-2008-obj1-*.oifits`) and the other, a model representing an AGN
  (`contest-2008-obj2-*.oifits`).  The total field of view of each data set
  is tapered with a 15 mas FWHM Gaussian.  Data are in the Optical/IR
  interferometry FITS format and are simulations of data obtained with the
  CHARA array.  Simulated data for each model are given in each H, J and K
  band.  The models for each band were generated from a set of components
  of differing spectra and thus differ in detail from band to band.  Each
  data set consists of data at 8 sub-bands with a constant brightness of
  the model with sub-band.

  In addition, there is data for a binary star (`contest-2008-binary.oifits`)
  to verify that the data has been read correctly, the parameters of
  this binary are:

   - Separation: 5.0 mas
   - PA (bright to faint): 30 degs PA east of north
   - Brightness ratio: 8.9
   - Uniform disk size of primary: 1.2 mas
   - Uniform disk size of secondary: 0.75 mas
