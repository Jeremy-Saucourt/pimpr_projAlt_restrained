# pimpr_projAlt
PIM-PR method for phasing laser arrays.

This archive is restrained due to proprietary code.

### Subfolders contents
  * **pimpr_projAlt_aleat**: PIM-PR method with random complex transfer matrices.
      
  * **pimpr_projAlt_gaussian_net_1D**: PIM-PR method with a propagating gaussian beam network. No diffractive element can be created. Propagation is modeled by gaussian beam equations. For speed purposes, calculations are made on a single dimension.
      
  * **pimpr_projAlt_objDiffrac_1D**: PIM-PR method with a diffractive element. Propagation is modeled by free-space transfer function and uses FFTs. For speed purposes, calculations are made on a single dimension.
      
  * **pimpr_projAlt_objDiffrac_2D**: PIM-PR method with a diffractive element. Propagation is modeled by free-space transfer function and uses FFTs. Here, calculations are operated in two dimensions.

  * **pimpr_projAlt_speckle_1D**: PIM-PR method with a suface diffuser. The diffuser is modeled as a random phase mask. Propagation is modeled by free-space transfer function and uses FFTs. For speed purposes, calculations are made on a single dimension.
