# VNP43NRTAlbedo.jl

Near-Real-Time Implementation of the VNP43 VIIRS BRDF Correction Algorithm for VNP09GA Surface Reflectance in Julia

The `VNP43NRTAlbedo.jl` package provides a Julia implementation of the Near-Real-Time (NRT) VIIRS Bidirectional Reflectance Distribution Function (BRDF) correction algorithm, specifically designed for VNP09GA Surface Reflectance data. This package allows for the calculation of key land surface albedo products, including **White-Sky Albedo (WSA)**, **Black-Sky Albedo (BSA)**, and **Nadir BRDF Adjusted Reflectance (NBAR)**, by inverting a BRDF model from observed reflectances and associated solar-view geometry.

## Contact

Margaret Johnson (she/her)<br>
[maggie.johnson@jpl.nasa.gov](mailto:maggie.johnson@jpl.nasa.gov)<br>
Algorithm implementation<br>
NASA Jet Propulsion Laboratory

Gregory H. Halverson (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
Algorithm re-implementation<br>
NASA Jet Propulsion Laboratory

Crystal Schaaf (she/her)<br>
[schaaf@bu.edu](mailto:schaaf@bu.edu)<br>
Algorithm inventor<br>
Boston University

Kerry Cawse-Nicholson (she/her)<br>
[kerry-anne.cawse-nicholson@jpl.nasa.gov](mailto:kerry-anne.cawse-nicholson@jpl.nasa.gov)<br>
Concept development and project management<br>
NASA Jet Propulsion Laboratory

## Features

* **BRDF Model Inversion**: Implements a linear least squares approach to invert the kernel-driven BRDF model (RossThick/LiSparseR).
* **Albedo and Nadir Reflectance Calculation**: Derives WSA, BSA, and NBAR from the fitted BRDF parameters.
* **Uncertainty Propagation**: Provides standard errors for the calculated albedo and nadir products based on the BRDF fit uncertainty.
* **Weighted Least Squares**: Supports weighted least squares regression with an exponential decay scheme, prioritizing more recent observations for NRT applications.
* **Solar Zenith Angle Calculation**: Includes functions to compute solar zenith angles from time and location.
* **Sinusoidal Tile Dimension Utility**: Provides a helper function to calculate spatial dimensions for MODIS/VIIRS sinusoidal tiles.
* **Compatibility**: Designed to work with `Rasters.jl` and `DimensionalData.jl` for robust geospatial data handling.

## Installation

To install `VNP43NRTAlbedo.jl`, open the Julia Pkg REPL by pressing `]` and run:

```julia
pkg> add VNP43NRTAlbedo
```

## Usage

The core functionality of the package is exposed through functions like `NRT_BRDF_all`, which takes observed reflectances and angular information to produce albedo and nadir products.

## Module Structure

The package is organized into several key components:

* **Constants**: `SINUSOIDAL_CRS` defines the Well-Known Text (WKT) for the Sinusoidal projection.
* **Coordinate Utilities**:
    * `sinusoidal_tile_dims`: Calculates `X` and `Y` dimensions for a specific sinusoidal tile, useful for `Rasters.jl` integration.
* **BRDF Model Kernels**:
    * `Kvol`, `Kvol_vec`, `Kvol_sc`: Implement the RossThick volumetric scattering kernel for matrix, vector, and scalar inputs.
    * `Kgeo`, `Kgeo_vec`, `Kgeo_sc`: Implement the LiSparseR geometric-optical scattering kernel for matrix, vector, and scalar inputs.
* **Solar Geometry**:
    * `zenith_from_solarnoon_time`, `zenith_from_solarnoon_vec`: Functions to calculate solar zenith angle from time and location.
    * `calculate_SZA`: A simplified SZA calculation based on day of year, hour, and latitude.
* **BRDF Inversion and Product Generation**:
    * `BRDFParameters`: A struct to hold the results of BRDF model inversion (coefficients, standard errors, R-squared, covariance).
    * `NRT_BRDF`: Performs the core BRDF model inversion for a given set of observations and kernels.
    * `NRT_BRDF_albedo`: Calculates White-Sky and Black-Sky Albedo.
    * `NRT_BRDF_nadir`: Calculates Nadir BRDF Adjusted Reflectance.
    * `NRT_BRDF_all`: A consolidated function to compute WSA, BSA, and Nadir reflectance along with their uncertainties and fit statistics.
* **Date Utilities**:
    * `date_range`: Generates a vector of `Date` objects within a specified range.
* **Internal Modules**:
    * `MODLAND.jl` and `VIIRS.jl`: These are included to encapsulate functionalities specific to MODIS and VIIRS satellite sensors, respectively (e.g., data structures, constants).

## Contributing

Contributions to `VNP43NRTAlbedo.jl` are welcome! Please refer to the guidelines for contributing to Julia packages.

## References

1. **Schaaf, C. B., Gao, F., Strahler, A. H., Lucht, W., Li, X., Tsang, T., Strugnell, N. C., Zhang, X., Jin, Y., Muller, J.-P., Lewis, P., Barnsley, M. J., Hobson, P., Disney, M., Roberts, G., Dunderdale, M., Doll, C., d’Entremont, R. P., Hu, B., Liang, S., Privette, J. L., & Roy, D. (2002).**  
   First operational BRDF, albedo nadir reflectance products from MODIS. *Remote Sensing of Environment*, 83(1-2), 135-148.  
   [https://doi.org/10.1016/S0034-4257(02)00091-3](https://doi.org/10.1016/S0034-4257(02)00091-3)

2. **Lucht, W., Schaaf, C. B., & Strahler, A. H. (2000).**  
   An algorithm for the retrieval of albedo from space using semiempirical BRDF models. *IEEE Transactions on Geoscience and Remote Sensing*, 38(2), 977-998.  
   [https://doi.org/10.1109/36.841980](https://doi.org/10.1109/36.841980)

3. **Schaaf, C. B., Liu, J., Gao, F., Strahler, A. H., & Jin, Y. (2011).**  
   VIIRS Bidirectional Reflectance Distribution Function and Albedo Climate Data Record (BRDF/Albedo CDR) Algorithm Theoretical Basis Document (ATBD).  
   [NOAA/NCEI CDR Documentation](https://www.ncei.noaa.gov/pub/data/sds/cdr/CDRs/VIIRS_Land_Surface_Reflectance/VNP43MA3/DOC/ATBD_VNP43MA3_CDR_v1.0.pdf)

4. **USGS MODIS & VIIRS Land Products:**  
   [MODIS BRDF/Albedo Product Overview](https://modis-land.gsfc.nasa.gov/brdf.html)  
   [VIIRS BRDF/Albedo Product Overview (VNP43)](https://lpdaac.usgs.gov/products/vnp43ma3v001/)

