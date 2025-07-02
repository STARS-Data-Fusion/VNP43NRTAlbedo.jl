module VNP43NRTAlbedo

using Dates
using LinearAlgebra
using Statistics
using Rasters
using DimensionalData.Dimensions.LookupArrays
using ProgressMeter

# Include modules defining MODIS and VIIRS specific functionalities.
# These modules likely contain data structures, helper functions, or constants
# relevant to processing data from these satellite sensors.
include("MODLAND.jl")
using .MODLAND

include("VIIRS.jl")
using .VIIRS


"""
    SINUSOIDAL_CRS

A constant representing the Well-Known Text (WKT) definition of the Sinusoidal projection
Coordinate Reference System (CRS). This projection is widely used for global land products,
including MODIS and VIIRS, due to its equal-area property.
"""
const SINUSOIDAL_CRS = WellKnownText("PROJCS[\"unknown\",GEOGCS[\"unknown\",DATUM[\"unknown\",SPHEROID[\"unknown\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"longitude_of_center\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]")


"""
    sinusoidal_tile_dims(h::Int, v::Int, tile_width_cells::Int)::Tuple{X,Y}

Calculates the `X` and `Y` dimensions (spatial coordinates) for a specific tile
within the global Sinusoidal projection grid.

# Arguments
- `h::Int`: The horizontal tile index (column) of the desired tile.
- `v::Int`: The vertical tile index (row) of the desired tile.
- `tile_width_cells::Int`: The number of cells (pixels) along the width (and height)
  of the tile. Assumes square cells.

# Returns
- `Tuple{X,Y}`: A tuple containing two `DimensionalData.jl` dimension objects, `X` and `Y`,
  which define the spatial extent and resolution of the tile in the Sinusoidal CRS.

# Details
This function uses hardcoded global projection boundaries and tile size information
(typical for MODIS/VIIRS sinusoidal tiles) to calculate the precise geographic
coordinates for the given tile. The `X` and `Y` dimensions are created as `Projected`
dimensions with `LinRange` for coordinates, `Regular` span, and `Start` sampling,
all referenced to the `SINUSOIDAL_CRS`.
"""
function sinusoidal_tile_dims(h::Int, v::Int, tile_width_cells::Int)::Tuple{X,Y}
    # Boundaries of the global sinusoidal projection in meters.
    # These are standard constants for the MODIS/VIIRS sinusoidal grid.
    GLOBE_UPPER_LEFT_X = -20015109.355798
    GLOBE_UPPER_LEFT_Y = 10007554.677899
    GLOBE_LOWER_RIGHT_X = 20015109.355798
    GLOBE_LOWER_RIGHT_Y = -10007554.677899

    # Size across (width or height) of any equal-area sinusoidal tile in meters.
    # This is also a standard constant.
    TILE_SIZE = 1111950.5197665554

    # Total number of rows and columns in the global sinusoidal tile grid.
    TOTAL_ROWS = 18
    TOTAL_COLUMNS = 36

    # Calculate the size of an individual cell (pixel) within the tile.
    cell_size = TILE_SIZE / tile_width_cells

    # Calculate the exact geographic (projected) coordinates for the given tile.
    # X coordinates: from left to right.
    tile_left_x = GLOBE_UPPER_LEFT_X + h * TILE_SIZE
    tile_right_x = GLOBE_UPPER_LEFT_X + (h + 1) * TILE_SIZE - cell_size # Subtract cell_size to get the end of the last cell, not beyond the tile.
    # Y coordinates: from upper to lower. Note the inversion for v index.
    # The 'v' index counts from top to bottom (0 to 17), but the Y-coordinates decrease upwards.
    # So, (TOTAL_ROWS - v) gives the number of tiles from the bottom.
    tile_upper_y = GLOBE_LOWER_RIGHT_Y + (TOTAL_ROWS - v) * TILE_SIZE - cell_size
    tile_lower_y = GLOBE_LOWER_RIGHT_Y + (TOTAL_ROWS - 1 - v) * TILE_SIZE

    # Define the sampling strategy for the dimensions (intervals start at the coordinate value).
    sampling = Intervals(Start())

    # Create the X dimension.
    # `LinRange` creates a range of evenly spaced points.
    # `ForwardOrdered` indicates coordinates increase along the dimension.
    # `Regular` span indicates a constant step size.
    x_dim = X(Projected(LinRange(tile_left_x, tile_right_x, tile_width_cells),
                        order=ForwardOrdered(), span=Regular(cell_size),
                        sampling=sampling, crs=SINUSOIDAL_CRS))

    # Create the Y dimension.
    # `ReverseOrdered` indicates coordinates decrease along the dimension (from top to bottom).
    # Note the negative `cell_size` for `Regular` span to reflect the decreasing order.
    y_dim = Y(Projected(LinRange(tile_upper_y, tile_lower_y, tile_width_cells),
                        order=ReverseOrdered(), span=Regular(-cell_size),
                        sampling=sampling, crs=SINUSOIDAL_CRS))

    dims = (x_dim, y_dim)

    return dims
end

# Export the function to make it accessible when the module is used.
export sinusoidal_tile_dims

"""
    BRDFParameters

A struct to encapsulate the results of a BRDF (Bidirectional Reflectance Distribution Function)
model inversion.

# Fields
- `brdf::AbstractMatrix`: A matrix where each row contains the fitted BRDF model
  coefficients (e.g., [f0, f1, f2] for isotropic, volumetric, and geometric kernels).
- `se::AbstractVector`: A vector of standard errors of the BRDF model fit for each location.
- `R2::AbstractVector`: A vector of R-squared values, indicating the coefficient of
  determination for each BRDF fit.
- `uncert::AbstractArray{Float64,3}`: A 3D array where `uncert[:, :, i]` represents the
  covariance matrix of the BRDF coefficients for the i-th location. This matrix is
  typically the inverse of (X'X) or (X'SX) in a least squares regression.
"""
struct BRDFParameters
    brdf::AbstractMatrix
    se::AbstractVector
    R2::AbstractVector
    uncert::AbstractArray{Float64,3}
end

"""
    zenith_from_solarnoon_time(tm::Float64, lat::Float64, lon::Float64)::Float64

Calculates the solar zenith angle (SZA) at a given time and location using
astronomical formulas. This function is designed for a single time point.

# Arguments
- `tm::Float64`: Time in seconds (likely since a specific epoch, e.g., J2000).
- `lat::Float64`: Latitude in degrees.
- `lon::Float64`: Longitude in degrees.

# Returns
- `Float64`: The solar zenith angle in degrees.

# Details
This function implements a standard algorithm for calculating solar position,
involving conversions to Julian Date, calculation of mean anomaly, eccentricity,
equation of time, solar declination, and finally the hour angle to derive the SZA.
It ensures the cosine of zenith angle is clamped between -1 and 1 to avoid `NaN` from `acos`.
"""
function zenith_from_solarnoon_time(tm::Float64, lat::Float64, lon::Float64)::Float64
    rad = π/180.0 # Degrees to radians conversion factor

    # Calculate Julian Date (Jd) from input time (assuming tm is seconds since J2000 for 2440587.5 base)
    Jd = tm/86400.0 + 2440587.5
    # Calculate Julian Century (Jc)
    Jc = (Jd - 2451545.0)/36525.0

    # Calculate Geocentric Mean Longitude (L0) of the Sun
    L0 = mod((280.46646 + Jc * (36000.76983 + 0.0003032 * Jc)),360)
    # Calculate Mean Anomaly (M) of the Sun
    M = 357.52911 + Jc * (35999.05029 - 0.0001537 * Jc)
    # Calculate Eccentricity (e) of Earth's orbit
    e = 0.016708634 - Jc * (4.2037e-05 + 1.267e-07 * Jc)

    # Calculate Equation of Center (eqctr)
    eqctr = sin(rad * M) * (1.914602 - Jc * (0.004817 + 1.4e-05 *
        Jc)) + sin(rad * 2.0 * M) * (0.019993 - 0.000101 * Jc) +
        sin(rad * 3.0 * M) * 0.000289

    # Calculate True Longitude (lambda0) of the Sun
    lambda0 = L0 + eqctr
    # Calculate Longitude of Ascending Node (omega)
    omega = 125.04 - 1934.136 * Jc
    # Calculate Apparent Longitude (lambda) of the Sun
    lambda = lambda0 - 0.00569 - 0.00478 * sin(rad * omega)

    # Calculate Obliquity of the Ecliptic (obliq)
    seconds = 21.448 - Jc * (46.815 + Jc * (0.00059 - Jc * (0.001813)))
    obliq0 = 23.0 + (26.0 + (seconds/60.0))/60.0
    obliq = obliq0 + 0.00256 * cos(rad * omega)

    # Calculate 'y' term for Equation of Time
    y = tan(rad * obliq/2)^2
    # Calculate Equation of Time (eqnTime) in minutes
    eqnTime = 4/rad * (y*sin(rad*2*L0) - 2*e*sin(rad*M) +
        4*e*y*sin(rad*M)*cos(rad*2*L0) - y^2/2*sin(rad*4*L0) -
        e^2*1.25*sin(rad*2*M))

    # Calculate Solar Declination (solarDec)
    solarDec = asin(sin(rad*obliq)*sin(rad*lambda))
    sinSolarDec = sin(solarDec)
    cosSolarDec = cos(solarDec)

    # Calculate Local Solar Time in minutes
    solarTime = (mod(Jd-1/2,1)*1440+eqnTime)/4
    # Calculate Hour Angle (hourAngle)
    hourAngle = solarTime+lon-180 # Note: Solar time usually in 0-360 range for calculations.

    # Calculate cosine of Solar Zenith Angle (cosZenith)
    cosZenith = sin(rad*lat)*sinSolarDec+cos(rad*lat)*cosSolarDec*cos(rad*hourAngle)

    # Clamp cosZenith to the valid range [-1, 1] to prevent domain errors in acos
    if cosZenith < -1
        cosZenith=-1
    end
    if cosZenith > 1
        cosZenith=1
    end

    # Return Solar Zenith Angle in degrees
    return acos(cosZenith)/rad
end

"""
    zenith_from_solarnoon_vec(tm::Vector{Float64}, lat::Vector{Float64}, lon::Vector{Float64})::Vector{Float64}

Calculates the solar zenith angle (SZA) for multiple time points and locations
by vectorizing the `zenith_from_solarnoon_time` function.

# Arguments
- `tm::Vector{Float64}`: A vector of times in seconds.
- `lat::Vector{Float64}`: A vector of latitudes in degrees.
- `lon::Vector{Float64}`: A vector of longitudes in degrees.

# Returns
- `Vector{Float64}`: A vector of solar zenith angles in degrees, corresponding
  to each input time and location.
"""
function zenith_from_solarnoon_vec(tm::Vector{Float64}, lat::Vector{Float64}, lon::Vector{Float64})::Vector{Float64}
    n = length(tm)
    szn = Vector{Float64}(undef, n) # Preallocate output vector

    # Loop through each element and calculate SZA using the scalar function
    for i in 1:n
        szn[i] = zenith_from_solarnoon_time(tm[i], lat[i], lon[i])
    end

    return szn
end

"""
    Kvol(sz::Matrix, vz::Matrix, rz::Matrix)::Matrix

Calculates the RossThick volumetric scattering kernel (`Kvol`) for given solar zenith (sz),
view zenith (vz), and relative azimuth (rz) angles. This version operates on matrices
of angles.

# Arguments
- `sz::Matrix`: Matrix of solar zenith angles (degrees).
- `vz::Matrix`: Matrix of view zenith angles (degrees).
- `rz::Matrix`: Matrix of relative azimuth angles (degrees).

# Returns
- `Matrix`: A matrix of RossThick volumetric kernel values.

# Details
The RossThick kernel describes the scattering of light due to a homogeneous,
randomly distributed canopy of leaves. The formula is based on geometric optics
and radiative transfer principles. `eps` represents the phase angle between
the incident and observed light rays.
"""
function Kvol(sz::Matrix, vz::Matrix, rz::Matrix)::Matrix
    sc = π / 180 # Conversion factor from degrees to radians
    # Calculate the phase angle (eps) in radians
    eps = acos.(cos.(sc .* sz) .* cos.(sc .* vz) .+ sin.(sc .* sz) .* sin.(sc .* vz) .* cos.(sc .* rz))
    # Calculate the RossThick volumetric kernel
    Kvol = 1 ./ (cos.(sc .* sz) .+ cos.(sc .* vz)) .* ((π / 2 .- eps) .* cos.(eps) .+ sin.(eps)) .- π / 4

    return Kvol
end

"""
    Kvol_vec(sz::AbstractVector, vz::AbstractVector, rz::AbstractVector)::AbstractVector

Calculates the RossThick volumetric scattering kernel (`Kvol`) for given solar zenith (sz),
view zenith (vz), and relative azimuth (rz) angles. This version operates on vectors
of angles.

# Arguments
- `sz::AbstractVector`: Vector of solar zenith angles (degrees).
- `vz::AbstractVector`: Vector of view zenith angles (degrees).
- `rz::AbstractVector`: Vector of relative azimuth angles (degrees).

# Returns
- `AbstractVector`: A vector of RossThick volumetric kernel values.
"""
function Kvol_vec(sz::AbstractVector, vz::AbstractVector, rz::AbstractVector)::AbstractVector
    sc = pi/180 # Conversion factor from degrees to radians
    # Calculate the phase angle (eps) in radians
    eps = acos.(cos.(sc .* sz) .* cos.(sc .* vz) .+ sin.(sc .* sz) .* sin.(sc .* vz) .* cos.(sc .* rz))
    # Calculate the RossThick volumetric kernel
    Kvol_vec = 1 ./ (cos.(sc .* sz) .+ cos.(sc .* vz)) .* ((pi/2 .- eps) .* cos.(eps) .+ sin.(eps)) .- pi/4

    return Kvol_vec
end

"""
    Kvol_sc(sz::Real, vz::Real, rz::Real)::Real

Calculates the RossThick volumetric scattering kernel (`Kvol`) for scalar solar zenith (sz),
view zenith (vz), and relative azimuth (rz) angles.

# Arguments
- `sz::Real`: Solar zenith angle (degrees).
- `vz::Real`: View zenith angle (degrees).
- `rz::Real`: Relative azimuth angle (degrees).

# Returns
- `Real`: The RossThick volumetric kernel value.
"""
function Kvol_sc(sz::Real, vz::Real, rz::Real)::Real
    sc = pi/180 # Conversion factor from degrees to radians
    # Calculate the phase angle (eps) in radians
    eps = acos(cos(sc*sz)*cos(sc*vz) + sin(sc*sz)*sin(sc*vz)*cos(sc*rz))
    # Calculate the RossThick volumetric kernel
    Kvol_sc = 1/(cos(sc*sz) + cos(sc*vz))*((pi/2-eps)*cos(eps) + sin(eps)) - pi/4

    return Kvol_sc
end

"""
    Kgeo_sc(sz::Real, vz::Real, rz::Real)::Real

Calculates the LiSparseR geometric-optical scattering kernel (`Kgeo`) for scalar solar zenith (sz),
view zenith (vz), and relative azimuth (rz) angles.

# Arguments
- `sz::Real`: Solar zenith angle (degrees).
- `vz::Real`: View zenith angle (degrees).
- `rz::Real`: Relative azimuth angle (degrees).

# Returns
- `Real`: The LiSparseR geometric-optical kernel value.

# Details
The LiSparseR kernel models the geometric effects of canopy structure, such as
shadowing and obscuring. It is based on a sparse canopy approximation.
`D` and `t` are intermediate variables related to the geometry, and `O` is the
overlap function.
"""
function Kgeo_sc(sz::Real, vz::Real, rz::Real)::Real
    sc = pi/180 # Conversion factor from degrees to radians
    # Calculate the phase angle (eps) in radians
    eps = acos(cos(sc*sz)*cos(sc*vz) + sin(sc*sz)*sin(sc*vz)*cos(sc*rz))

    # Calculate D, a geometric factor related to the distance between sun and sensor directions.
    D = sqrt(tan(sc*sz)^2 + tan(sc*vz)^2 - 2*tan(sc*sz)*tan(sc*vz)*cos(sc*rz))

    # Calculate 'cost', related to the cosine of the angle between sun and sensor.
    cost = 2*sqrt(D^2 + (tan(sc*sz)*tan(sc*vz)*sin(sc*rz))^2)/(1/cos(sc*sz) + 1/cos(sc*vz))
    cost = clamp(cost, -1.0, 1.0) # Clamp to valid range for acos.

    # Calculate 't', the angle related to cost.
    t = acos(cost)

    # Calculate O, the overlap function.
    O = 1/pi*(t-sin(t)*cos(t))*(1/cos(sc*sz) + 1/cos(sc*vz))
    O = max(O, 0) # Ensure O is non-negative.

    # Calculate the LiSparseR geometric-optical kernel.
    Kgeo_sc = O - 1/cos(sc*sz) - 1/cos(sc*vz) + 0.5*(1+cos(eps))/cos(sc*sz)/cos(sc*vz)

    return Kgeo_sc
end

"""
    Kgeo(sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix)::Matrix

Calculates the LiSparseR geometric-optical scattering kernel (`Kgeo`) for given solar zenith (sz),
view zenith (vz), and relative azimuth (rz) angles. This version operates on matrices
of angles.

# Arguments
- `sz::AbstractMatrix`: Matrix of solar zenith angles (degrees).
- `vz::AbstractMatrix`: Matrix of view zenith angles (degrees).
- `rz::AbstractMatrix`: Matrix of relative azimuth angles (degrees).

# Returns
- `Matrix`: A matrix of LiSparseR geometric-optical kernel values.
"""
function Kgeo(sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix)::Matrix
    n, t = size(sz) # Get dimensions (n: number of locations, t: number of observations/times)
    Kg = zeros(n, t) # Preallocate output matrix

    # Loop through each element and calculate Kgeo using the scalar function
    for i in 1:n
        for j in 1:t
            Kg[i,j] = Kgeo_sc(sz[i,j], vz[i,j], rz[i,j])
        end
    end
    return Kg
end

"""
    Kgeo_vec(sz::AbstractVector, vz::AbstractVector, rz::AbstractVector)::AbstractVector

Calculates the LiSparseR geometric-optical scattering kernel (`Kgeo`) for given solar zenith (sz),
view zenith (vz), and relative azimuth (rz) angles. This version operates on vectors
of angles.

# Arguments
- `sz::AbstractVector`: Vector of solar zenith angles (degrees).
- `vz::AbstractVector`: Vector of view zenith angles (degrees).
- `rz::AbstractVector`: Vector of relative azimuth angles (degrees).

# Returns
- `AbstractVector`: A vector of LiSparseR geometric-optical kernel values.
"""
function Kgeo_vec(sz::AbstractVector, vz::AbstractVector, rz::AbstractVector)::AbstractVector
    t = length(vz) # Get number of observations/times
    Kg = zeros(t) # Preallocate output vector

    # Loop through each element and calculate Kgeo using the scalar function
    for j in 1:t
        Kg[j] = Kgeo_sc(sz[j], vz[j], rz[j])
    end

    return Kg
end

"""
    NRT_BRDF(Y::AbstractMatrix, kv::AbstractMatrix, kg::AbstractMatrix, weighted::Bool, scale::Real)::BRDFParameters

Performs Near Real-Time (NRT) BRDF model inversion using a linear least squares approach.
This function estimates the BRDF model parameters (coefficients) for each location.

# Arguments
- `Y::AbstractMatrix`: A matrix of observed reflectances. Each row represents a spatial
  location (e.g., pixel), and each column represents a different observation (e.g., time
  or acquisition).
- `kv::AbstractMatrix`: A matrix of RossThick volumetric kernel values, corresponding to `Y`.
- `kg::AbstractMatrix`: A matrix of LiSparseR geometric-optical kernel values, corresponding to `Y`.
- `weighted::Bool`: If `true`, a weighted least squares regression is performed.
  If `false`, an ordinary least squares regression is performed.
- `scale::Real`: A parameter used in the exponential weighting scheme if `weighted` is `true`.
  It controls the decay rate of weights for older observations.

# Returns
- `BRDFParameters`: A struct containing the fitted BRDF coefficients (`brdf`),
  standard errors of the fit (`se`), R-squared values (`R2`), and the covariance
  matrix of the coefficients (`uncert`) for each location.

# Details
The BRDF model assumed is `Y = f0 + f1 * Kvol + f2 * Kgeo`, where `f0`, `f1`, `f2`
are the BRDF parameters to be estimated. The function iterates through each row
(location) of the input data:
1. Identifies non-missing observations.
2. Requires a minimum of 7 valid observations for a stable fit.
3. Constructs the design matrix `x = [1, Kvol, Kgeo]`.
4. Performs (weighted) least squares regression to find `brdf` coefficients.
5. Calculates `se` (residual standard error) and `R2`.
6. Calculates `uncert` (covariance matrix of the coefficients).
"""
function NRT_BRDF(Y::AbstractMatrix, kv::AbstractMatrix, kg::AbstractMatrix, weighted::Bool, scale::Real)::BRDFParameters
    n, p = size(Y) # n: number of locations, p: number of observations/times per location
    brdf = fill(NaN, n, 3) # Preallocate matrix for BRDF coefficients (3 per location: f0, f1, f2)
    R2 = fill(NaN, n)      # Preallocate vector for R-squared values
    se = fill(NaN, n)      # Preallocate vector for standard errors of the fit
    uncert = fill(NaN, 3, 3, n) # Preallocate 3D array for covariance matrices (3x3 per location)

    # Initialize the design matrix 'x' for the linear regression.
    # It will have 3 rows (for intercept, Kvol, Kgeo) and 'p' columns (for observations).
    x = ones(3, p)

    # Iterate over each spatial location (row in Y)
    for i in 1:n
        yt = Y[i,:] # Reflectance observations for the current location
        non_missing = findall(isfinite.(yt)) # Find indices of finite (non-NaN) observations
        nt = length(non_missing) # Number of valid observations

        # Skip locations with insufficient valid observations for a stable fit
        if nt < 7
            continue
        else
            # Populate the design matrix with Kvol and Kgeo for the current location
            x[2,:] = kv[i,:] # Kvol values
            x[3,:] = kg[i,:] # Kgeo values

            # Select only the valid observations and corresponding design matrix columns
            xt = x[:,non_missing]
            ytt = yt[non_missing]

            # Perform (weighted) least squares regression
            if weighted
                # Calculate weights using an exponential decay based on observation order.
                # Older observations get lower weights.
                xx = exp.(-0.5 .* range(p-1, stop=0; length=p) ./ scale)
                xxt = xx[non_missing] # Weights for non-missing observations

                SSi = Diagonal(xxt) # Create a diagonal matrix from weights
                Si = inv(xt * SSi * xt') # Weighted inverse of (X'WX) for covariance
                brdf[i,:] = (Si * xt * SSi * ytt)' # Weighted least squares solution
                yp = (brdf[i,:] * xt)' # Predicted reflectances
                # Weighted standard error calculation
                se[i] = sqrt(sum(((ytt - yp) .* xxt).^2)/(nt-3)) # (nt-3) because 3 parameters are fitted
            else
                Si = inv(xt * xt') # Inverse of (X'X) for covariance (Ordinary Least Squares)
                brdf[i,:] = (Si * xt * ytt)' # Ordinary least squares solution
                yp = (brdf[i,:] * xt)' # Predicted reflectances
                # Ordinary standard error calculation
                se[i] = sqrt(sum((ytt - yp).^2)/(nt-3))
            end

            # Store the estimated covariance matrix (Si) of the BRDF coefficients
            uncert[:,:,i] = Si
            # Calculate R-squared value
            # R2 = 1 - (residual sum of squares / total sum of squares)
            # The original formula uses (nt-3) in numerator and denominator,
            # which is an adjusted R2-like measure for error variance estimation.
            # A more common R2 is 1 - SS_res / SS_tot.
            # Here it seems to calculate R2 = 1 - (estimated error variance / sample variance).
            R2[i] = 1 - se[i]^2 * (nt-3) / var(ytt; corrected=true) / nt
        end
    end

    return BRDFParameters(brdf, se, R2, uncert)
end

"""
    NRT_BRDF_albedo(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool, scale::Real)::AbstractMatrix

Calculates White-Sky Albedo (WSA) and Black-Sky Albedo (BSA) from observed reflectances
and view/solar geometry using a BRDF model inversion.

# Arguments
- `Y::AbstractMatrix`: Matrix of observed reflectances (rows: locations, cols: observations).
- `sz::AbstractMatrix`: Matrix of solar zenith angles (degrees) corresponding to `Y`.
- `vz::AbstractMatrix`: Matrix of view zenith angles (degrees) corresponding to `Y`.
- `rz::AbstractMatrix`: Matrix of relative azimuth angles (degrees) corresponding to `Y`.
- `soz_noon::AbstractVector`: Vector of solar zenith angles at local solar noon (degrees) for each location.
  This is used for calculating Black-Sky Albedo.
- `weighted::Bool`: If `true`, weighted least squares is used for BRDF inversion.
- `scale::Real`: Weighting scale parameter if `weighted` is `true`.

# Returns
- `AbstractMatrix`: A matrix where each row corresponds to a location and columns contain:
  [WSA, BSA, WSA_SE, BSA_SE, BRDF_RMSE, BRDF_R2, NumberOfObservations].

# Details
This function extends `NRT_BRDF` by calculating specific albedo products:
- **White-Sky Albedo (WSA)**: Directional-hemispherical reflectance, representing albedo
  under isotropic (diffuse) illumination. It's calculated using fixed coefficients `gwsa`.
- **Black-Sky Albedo (BSA)**: Bi-hemispherical reflectance, representing albedo
  under direct illumination at a specific solar zenith angle (here, `soz_noon`).
  Its coefficients `gbsa` are dependent on `soz_noon`.
The standard errors for WSA and BSA are propagated from the BRDF parameter uncertainties.
"""
function NRT_BRDF_albedo(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool, scale::Real)::AbstractMatrix
    # RossThick kernel constants for albedo calculation
    g0vol = -0.007574
    g1vol = -0.070987
    g2vol = 0.307588

    # LiSparseR kernel constants for albedo calculation
    g0geo = -1.284909
    g1geo = -0.166314
    g2geo = 0.041840

    # Coefficients for White-Sky Albedo (WSA)
    # [isotropic, volumetric, geometric]
    gwsa = [1.0, 0.189184, -1.377622]

    # Coefficients for Black-Sky Albedo (BSA). The second and third elements
    # are calculated based on solar zenith angle at noon.
    gbsa = [1.0, 0.0, 0.0]

    n, p = size(Y) # n: number of locations, p: number of observations/times
    @info "reflectance n: $(n) p: $(p)"
    # Preallocate results matrix: (wsa, bsa, wsa_se, bsa_se, rmse, R2, nt)
    results = fill(NaN, n, 7)
    @info "results shape rows: $(size(results)[1]) cols: $(size(results)[2])"

    # Initialize the design matrix 'x' for the linear regression.
    x = ones(3, p)

    # Iterate over each spatial location
    for i in 1:n
        yt = Y[i,:] # Reflectance observations for the current location
        non_missing = findall(isfinite.(yt)) # Find indices of finite observations
        nt = length(non_missing) # Number of valid observations

        # Skip locations with insufficient valid observations
        if nt < 7
            continue
        else
            # Calculate BSA kernel coefficients based on solar zenith angle at local solar noon
            sznrad = soz_noon[i] * pi/180 # Convert to radians
            gbsa[2] = g0vol + g1vol * sznrad^2 + g2vol * sznrad^3 # Volumetric component
            gbsa[3] = g0geo + g1geo * sznrad^2 + g2geo * sznrad^3 # Geometric component

            # Populate the design matrix with Kvol and Kgeo for the current location
            x[2,:] = Kvol_vec(sz[i,:], vz[i,:], rz[i,:])
            x[3,:] = Kgeo_vec(sz[i,:], vz[i,:], rz[i,:])

            # Select only valid observations and corresponding design matrix columns
            xt = x[:,non_missing]
            ytt = yt[non_missing]

            # Perform BRDF inversion (weighted or unweighted least squares)
            if weighted
                xx = exp.(-0.5 .* range(p-1, stop=0; length=p) ./ scale)
                xxt = xx[non_missing]
                SSi = Diagonal(xxt)
                Si = inv(xt * SSi * xt') # Covariance matrix
                brdf = Si * xt * SSi * ytt # Fitted BRDF coefficients
                yp = (brdf' * xt)' # Predicted reflectances
                se = sqrt(sum(((ytt - yp) .* xxt).^2)/(nt-3)) # Standard error of fit
            else
                Si = inv(xt * xt') # Covariance matrix
                brdf = Si * xt * ytt # Fitted BRDF coefficients
                yp = (brdf' * xt)' # Predicted reflectances
                se = sqrt(sum((ytt - yp).^2)/(nt-3)) # Standard error of fit
            end

            # Calculate WSA and BSA using the fitted BRDF coefficients and their respective 'g' vectors
            results[i,1] = dot(gwsa, brdf) # White-Sky Albedo (WSA)
            results[i,2] = dot(gbsa, brdf) # Black-Sky Albedo (BSA)

            # Propagate uncertainty to WSA and BSA (se * sqrt(g' * Si * g))
            results[i,3] = se * sqrt(dot(gwsa, Si * gwsa)) # WSA standard error
            results[i,4] = se * sqrt(dot(gbsa, Si * gbsa)) # BSA standard error

            # Store fit statistics
            results[i,5] = se # BRDF Root Mean Square Error (RMSE)
            results[i,6] = 1 - se^2 * (nt-3) / var(ytt; corrected=true) / nt # BRDF R-squared
            results[i,7] = nt # Number of observations used for BRDF estimation
        end
    end

    return results
end

"""
    NRT_BRDF_nadir(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool, scale::Real)::AbstractMatrix

Calculates Nadir (0-degree view zenith, 0-degree view azimuth) reflectance from
observed reflectances and view/solar geometry using a BRDF model inversion.

# Arguments
- `Y::AbstractMatrix`: Matrix of observed reflectances (rows: locations, cols: observations).
- `sz::AbstractMatrix`: Matrix of solar zenith angles (degrees) corresponding to `Y`.
- `vz::AbstractMatrix`: Matrix of view zenith angles (degrees) corresponding to `Y`.
- `rz::AbstractMatrix`: Matrix of relative azimuth angles (degrees) corresponding to `Y`.
- `soz_noon::AbstractVector`: Vector of solar zenith angles at local solar noon (degrees) for each location.
  This is used to calculate the kernel values for nadir view at solar noon.
- `weighted::Bool`: If `true`, weighted least squares is used for BRDF inversion.
- `scale::Real`: Weighting scale parameter if `weighted` is `true`.

# Returns
- `AbstractMatrix`: A matrix where each row corresponds to a location and columns contain:
  [Nadir, Nadir_SE, BRDF_RMSE, BRDF_R2, NumberOfObservations].

# Details
Similar to `NRT_BRDF_albedo`, this function performs a BRDF inversion and then
uses specific kernel values to calculate nadir reflectance. The `Knadir` coefficients
are computed for a nadir view (view zenith = 0, relative azimuth = 0) at the solar
zenith angle corresponding to local solar noon.
"""
function NRT_BRDF_nadir(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool, scale::Real)::AbstractMatrix
    n, p = size(Y) # n: number of locations, p: number of observations/times
    # Preallocate results matrix: (nadir, nadir_se, rmse, R2, nt)
    results = fill(NaN, n, 5)
    # Initialize the design matrix 'x' for the linear regression.
    x = ones(3, p)
    # Coefficients for Nadir reflectance calculation.
    # The volumetric and geometric components will be calculated dynamically.
    Knadir = [1.0, 0.0, 0.0]

    # Iterate over each spatial location
    for i in 1:n
        yt = Y[i,:] # Reflectance observations for the current location
        non_missing = findall(isfinite.(yt)) # Find indices of finite observations
        nt = length(non_missing) # Number of valid observations

        # Skip locations with insufficient valid observations
        if nt < 7
            continue
        else
            szn = soz_noon[i] # Solar zenith angle at local solar noon
            # Populate the design matrix with Kvol and Kgeo for the current location
            x[2,:] = Kvol_vec(sz[i,:], vz[i,:], rz[i,:])
            x[3,:] = Kgeo_vec(sz[i,:], vz[i,:], rz[i,:])

            # Select only valid observations and corresponding design matrix columns
            xt = x[:,non_missing]
            ytt = yt[non_missing]

            # Perform BRDF inversion (weighted or unweighted least squares)
            if weighted
                xx = exp.(-0.5 .* range(p-1, stop=0; length=p) ./ scale)
                xxt = xx[non_missing]
                SSi = Diagonal(xxt)
                Si = inv(xt * SSi * xt') # Covariance matrix
                brdf = Si * xt * SSi * ytt # Fitted BRDF coefficients
                yp = (brdf' * xt)' # Predicted reflectances
                se = sqrt(sum(((ytt - yp) .* xxt).^2)/(nt-3)) # Standard error of fit
            else
                Si = inv(xt * xt') # Covariance matrix
                brdf = Si * xt * ytt # Fitted BRDF coefficients
                yp = (brdf' * xt)' # Predicted reflectances
                se = sqrt(sum((ytt - yp).^2)/(nt-3)) # Standard error of fit
            end

            # Calculate Nadir kernel coefficients at solar noon for nadir view (vz=0, rz=0)
            Knadir[2] = Kvol_sc(szn, 0.0, 0.0)
            Knadir[3] = Kgeo_sc(szn, 0.0, 0.0)

            # Calculate Nadir reflectance and its standard error
            results[i,1] = dot(Knadir, brdf) # Nadir reflectance
            results[i,2] = se * sqrt(dot(Knadir, Si * Knadir)) # Nadir standard error

            # Store fit statistics
            results[i,3] = se # BRDF RMSE
            results[i,4] = 1 - se^2 * (nt-3) / var(ytt; corrected=true) / nt # BRDF R-squared
            results[i,5] = nt # Number of observations used
        end
    end

    return results
end

"""
    NRT_BRDF_all(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool = true, scale::Real = 1.87)::AbstractMatrix

Performs comprehensive Near Real-Time (NRT) BRDF model inversion and calculates
White-Sky Albedo (WSA), Black-Sky Albedo (BSA), and Nadir reflectance for each
spatial location.

# Arguments
- `Y::AbstractMatrix`: A matrix of observed reflectances. Rows are spatial locations,
  columns are observations over time.
- `sz::AbstractMatrix`: Solar zenith angles (degrees) corresponding to `Y`.
- `vz::AbstractMatrix`: View zenith angles (degrees) corresponding to `Y`.
- `rz::AbstractMatrix`: Relative azimuth angles (degrees) corresponding to `Y`.
- `soz_noon::AbstractVector`: Solar zenith angles at local solar noon (degrees) for each location.
- `weighted::Bool`: (Optional) If `true` (default), weighted least squares is used.
- `scale::Real`: (Optional) Weighting scale parameter (default: 1.87) if `weighted` is `true`.

# Returns
- `AbstractMatrix`: A matrix where each row corresponds to a location.
  The columns contain:
  [WSA, BSA, Nadir, WSA_SE, BSA_SE, Nadir_SE, BRDF_RMSE, BRDF_R2, NumberOfObservations].

# Details
This function combines the logic from `NRT_BRDF_albedo` and `NRT_BRDF_nadir`.
It iterates through each location, fits the BRDF model using (optionally weighted)
least squares, and then derives the three main BRDF-derived products along
with their propagated standard errors and fit statistics. A progress bar is
displayed during the computation for user feedback.
"""
function NRT_BRDF_all(Y::AbstractMatrix, sz::AbstractMatrix, vz::AbstractMatrix, rz::AbstractMatrix, soz_noon::AbstractVector, weighted::Bool = true, scale::Real = 1.87)::AbstractMatrix
    @info "processing BRDF"
    @info "reflectance rows: $(size(Y)[1]) cols: $(size(Y)[2])"
    @info "solar zenith rows: $(size(sz)[1]) cols: $(size(sz)[2])"
    @info "sensor zenith rows: $(size(vz)[1]) cols: $(size(vz)[2])"
    @info "relative azimuth rows: $(size(rz)[1]) cols: $(size(rz)[2])"
    @info "solar zenith noon size: $(size(soz_noon)[1])"

    # RossThick kernel constants for albedo/nadir calculation
    g0vol = -0.007574
    g1vol = -0.070987
    g2vol = 0.307588

    # LiSparseR kernel constants for albedo/nadir calculation
    g0geo = -1.284909
    g1geo = -0.166314
    g2geo = 0.041840

    # Coefficients for White-Sky Albedo (WSA)
    gwsa = [1.0, 0.189184, -1.377622]

    # Coefficients for Black-Sky Albedo (BSA) - second and third will be calculated dynamically
    gbsa = [1.0, 0.0, 0.0]
    # Coefficients for Nadir reflectance - second and third will be calculated dynamically
    gnbar = [1.0, 0.0, 0.0]

    n, p = size(Y) # n: number of locations, p: number of observations/times
    @info "n: $(n) p: $(p)"
    # Preallocate results matrix: (wsa, bsa, nadir, wsa_se, bsa_se, nadir_se, rmse, R2, nt)
    results = fill(NaN, n, 9)
    @info "results rows: $(size(results)[1]) cols: $(size(results)[2])"

    # Initialize the design matrix 'x' for the linear regression.
    x = ones(3, p)

    # Use ProgressMeter to display a progress bar during the loop
    @showprogress for i in 1:n
        yt = Y[i,:] # Reflectance observations for the current location
        non_missing = findall(isfinite.(yt)) # Find indices of finite observations
        nt = length(non_missing) # Number of valid observations

        # Skip locations with insufficient valid observations
        if nt < 7
            continue
        else
            sznrad = soz_noon[i] * pi/180 # Solar zenith angle at local solar noon in radians
            # Calculate BSA kernel coefficients based on solar zenith angle at local solar noon
            gbsa[2] = g0vol + g1vol * sznrad^2 + g2vol * sznrad^3
            gbsa[3] = g0geo + g1geo * sznrad^2 + g2geo * sznrad^3

            # Calculate Nadir kernel coefficients at solar noon for nadir view (vz=0, rz=0)
            gnbar[2] = Kvol_sc(soz_noon[i], 0.0, 0.0)
            gnbar[3] = Kgeo_sc(soz_noon[i], 0.0, 0.0)

            # Populate the design matrix with Kvol and Kgeo for the current location
            x[2,:] = Kvol_vec(sz[i,:], vz[i,:], rz[i,:])
            x[3,:] = Kgeo_vec(sz[i,:], vz[i,:], rz[i,:])

            # Select only valid observations and corresponding design matrix columns
            xt = x[:,non_missing]
            ytt = yt[non_missing]

            # Perform BRDF inversion (weighted or unweighted least squares)
            if weighted
                # Calculate weights using an exponential decay.
                xx = exp.(-0.5 .* range(p-1, stop=0; length=p) ./ scale)
                xxt = xx[non_missing]
                SSi = Diagonal(xxt) # Diagonal matrix of weights
                Si = inv(xt * SSi * xt') # Weighted inverse of (X'WX)
                brdf = Si * xt * SSi * ytt # Weighted least squares solution
                yp = (brdf' * xt)' # Predicted reflectances
                se = sqrt(sum(((ytt - yp) .* xxt).^2)/(nt-3)) # Weighted standard error of fit
            else
                Si = inv(xt * xt') # Inverse of (X'X)
                brdf = Si * xt * ytt # Ordinary least squares solution
                yp = (brdf' * xt)' # Predicted reflectances
                se = sqrt(sum((ytt - yp).^2)/(nt-3)) # Ordinary standard error of fit
            end

            # Calculate and store WSA, BSA, Nadir and their standard errors
            results[i,1] = dot(gwsa, brdf) # White-Sky Albedo (WSA)
            results[i,2] = dot(gbsa, brdf) # Black-Sky Albedo (BSA)
            results[i,3] = dot(gnbar, brdf) # Nadir reflectance

            results[i,4] = se * sqrt(dot(gwsa, Si * gwsa)) # WSA standard error
            results[i,5] = se * sqrt(dot(gbsa, Si * gbsa)) # BSA standard error
            results[i,6] = se * sqrt(dot(gnbar, Si * gnbar)) # Nadir standard error

            # Store fit statistics
            results[i,7] = se # BRDF RMSE
            results[i,8] = 1 - se^2 * (nt-3) / var(ytt; corrected=true) / nt # BRDF R-squared
            results[i,9] = nt # Number of observations used for BRDF estimation

            # Replace any `missing` values (if they somehow appeared, though fill(NaN) should prevent) with NaN.
            # This line might be redundant if NaN is always used instead of Missing.
            replace!(results, missing=>NaN)
        end
    end

    return results
end

# Export the main function for external use.
export NRT_BRDF_all

"""
    date_range(start::Date, stop::Date)::Vector{Date}

Generates a vector of `Date` objects starting from `start` and ending at `stop`,
inclusive, with a daily step.

# Arguments
- `start::Date`: The starting date.
- `stop::Date`: The ending date.

# Returns
- `Vector{Date}`: A vector of dates.
"""
function date_range(start::Date, stop::Date)::Vector{Date}
    return collect(start:Day(1):stop) # Use `collect` to turn the Date range into a vector
end

"""
    calculate_SZA(day_of_year::Array{Int64,1}, hour_of_day::Union{Array{Float64,1}, Float64}, lat::Array{Float64,1})::Array{Float64,1}

Calculates the Solar Zenith Angle (SZA) for given day(s) of the year, hour(s) of the day,
and latitude(s). This is a simplified astronomical calculation.

# Arguments
- `day_of_year::Array{Int64,1}`: A vector of day of the year (1-365/366).
- `hour_of_day::Union{Array{Float64,1}, Float64}`: A vector or scalar representing the
  hour of the day (0.0 - 24.0).
- `lat::Array{Float64,1}`: A vector of latitudes in degrees.

# Returns
- `Array{Float64,1}`: A vector of calculated Solar Zenith Angles in degrees.

# Details
This function uses approximate formulas to calculate solar declination and hour angle
to derive the SZA.
- `day_angle`: Angular position of Earth in its orbit.
- `dec`: Solar declination angle (angle of the sun north or south of the equator).
- `hour_angle`: Angle of the sun relative to the local meridian.
"""
function calculate_SZA(day_of_year::Array{Int64,1}, hour_of_day::Union{Array{Float64,1}, Float64}, lat::Array{Float64,1})::Array{Float64,1}
    # Convert day of year to an angle in radians
    day_angle = (2 * pi * (day_of_year .- 1)) / 365
    # Convert latitude from degrees to radians
    lat = deg2rad.(lat)
    # Calculate solar declination angle in radians, then convert to degrees
    dec = deg2rad.((0.006918 .- 0.399912 .* cos.(day_angle) .+ 0.070257 .* sin.(day_angle) .- 0.006758 .* cos.(2 .* day_angle) .+ 0.000907 .* sin.(2 .* day_angle) .- 0.002697 .* cos.(3 .* day_angle) .+ 0.00148 .* sin.(3 .* day_angle)) .* (180 / pi))
    # Calculate hour angle (15 degrees per hour, centered at solar noon 12:00)
    hour_angle = deg2rad.(hour_of_day * 15.0 .- 180.0) # Assuming hour_of_day is 0-24 and 12 is solar noon.

    # Calculate Solar Zenith Angle using the cosine formula
    # cos(SZA) = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hour_angle)
    SZA = rad2deg.(acos.(sin.(lat) .* sin.(dec) .+ cos.(lat) .* cos.(dec) .* cos.(hour_angle)))
    return SZA
end

end # module VNP43NRTAlbedo