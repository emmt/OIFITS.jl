# Things to be done

`OIFITS` package is *work-in-progress*, the following changes/improvements may
be implemented in a near future:

- Provide filter method to extract sub-sets of OI-FITS data: to select a chosen
  wavelength range, a given target, etc.

- Automatically rewrite `sta_index`, and indices in correlation matrices so
  that they match the indices in the related arrays.

- Automatically deal with revision number of `OI_TARGET` instances.

- If a data-block has linked dependencies, push them before pushing the
  data-block itself.

- Optionally, do not copy linked dependencies when copying a data-block.

- Use `EasyFITS` package and extend the `write!` method.
