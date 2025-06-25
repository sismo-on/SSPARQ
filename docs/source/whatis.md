# Project Overview

<span style="display:block;text-align:center">![image](_static/images/logo_mini.png)

[Python](https://www.python.org/)-based solution for processing seismological data to estimate:

- Sensor misorientation (deviation from true north)
- Clock synchronization errors (timing irregularities)
- Instrumental gain variations (component sensitivity differences)

through P-wave polarization analysis in three-component seismic recordings.

## Methodological Approach

### Theoretical Basis

The method assumes a homogeneous, isotropic medium beneath the station, where:

- P-wave energy should appear primarily on the vertical (Z) and radial (R) components
- Transverse (T) component energy should be minimal

### Processing Workflow

**1. Coordinate Transformation**

Raw seismic signal components (Z, N, E) are transformed into the ZRT coordinate system (Vertical, Radial, Transverse) using earthquake-to-station backazimuth information.

**2. Optimal Orientation Estimation**

Test rotations from 0° to 180° to:
  - Minimize transverse-component energy (P-wave leakage)
  - Maximize vertical-radial correlation (Pearson coefficient positive)

Resolve 180° ambiguity using cross-component phase analysis.

**3. Quality Control**

Apply strict validation criteria:
- Signal-to-noise ratio (SNR) on vertical component
- Transverse-to-radial (T/R) energy ratio
- Radial-to-vertical (R/Z) energy ratio

**4. Statistical Robustness**

Aggregate results from multiple events across different backazimuths to:
- Reduce bias from local structures
- Improve measurement reliability

**5. Trend Analysis**

Apply DBSCAN (Density-Based Spatial Clustering of Applications with Noise) clustering to identify:
- Consistent misorientation trends
- Outliers (potentially indicating timing/gain issues)

## Key Advantages

- Python Implementation: Leverages Python's scientific computing ecosystem
- Multi-parameter analysis: Simultaneously evaluates multiple instrumentation parameters
- Empirical Validation: Uses real-world seismic events for calibration