# Useful functions

## Find sensor misorientation

```{eval-rst}
.. py:function:: find_orientation(baz, SS, CCRZ, ETER, EREZ):

    This function calculates the best back azimuth (φ) and sensor misorientation (θ) based on the 
    given quality criteria: 
    
    - Transverse signal strength (SS),
    - Similarity of radial and vertical components (CCRZ), 
    - Transverse-to-radial energy ratio (ETER), and 
    - Radial-to-vertical energy ratio (EREZ).

    The cost function combines these criteria in such a way that minimazing the cost function helps to
    find the optimal back azimuth and corresponding orientation. The function outputs the best back azimuth, orientation, and the values of the quality criteria at the best azimuth index.

    :type baz: float
    :param baz: Initial back-azimuth from taup model (degrees)
    :type SS: numpy.ndarray
    :param SS: Transverse signal strength array (0 to 1 normalized)
    :type CCRZ: numpy.ndarray
    :param CCRZ: Vertical-radial similarity array (-1 to 1 normalized)
    :type ETER: numpy.ndarray
    :param ETER: Transverse-to-radial energy ratios
    :type EREZ: numpy.ndarray
    :param EREZ: Radial-to-vertical energy ratios

    :returns: φ, θ, SS_best, SZR_best, ETERR_best, EREZ_best
    :rtype: tuple of (float, float, float, float, float, float)

```

```{tip}
For optimal results:

- Use 0.1° azimuth search increments
- Normalize before optimization
```

```{admonition} Formula
The optimization minimazes:

$$
\mathcal{Cost}(\phi) = SS_{T} - CC_{RZ}
$$

Where:
- $SS_{T}$ = Normalized transverse component energy
- $CC_{RZ}$ = Radial-to-vertical Pearson correlation
```

-------------

## Applying the algorithm of Braunmiller & Pornsopin


```{eval-rst}
.. py:function:: Braunmiller_Pornsopin_algorithm(tr1, tr2, trZ, noise, baz, time_ins, CCVR_MIN=0.45, SNR_MIN=10, TRR_MIN=0.45, RVR_MIN=-1):

    Implements the P-wave particle motion analysis for:

    - Back-azimuth estimation
    - Sensor misorientation detection
    - Time instability
    - Instrument gain calculation
    - Data quality assessment

    Applies quality criteria with configurable thresholds.

    :type tr1: numpy.ndarray
    :param tr1: First horizontal component (typically North)
    :type tr2: numpy.ndarray
    :param tr2: Second horizontal component (typically East)
    :type trZ: numpy.ndarray
    :param trZ: Vertical component
    :type noise: numpy.ndarray
    :param noise: Noise window for SNR calculation
    :type baz: float
    :param baz: Theoretical back-azimuth (degrees)
    :type time_ins: float
    :param time_ins: Observed-predicted time difference (seconds)
    :returns: Dictionary containing 15 result metrics
    :rtype: dict
```

```{important}
**Quality Criteria Thresholds**:

1. SNR ≥ 10 dB
2. Vertical-radial correlation ≥ 0.5
3. Transverse/radial energy ratio ≤ 0.2
4. -90 sec < Time residual |Δt| < 90 sec
5. Radial/vertical energy ratio ≤ 2
```

```{admonition} Component rotation implemented as:

$$
\begin{pmatrix}
R \\ 
T 
\end{pmatrix} = 
\begin{pmatrix}
\cos\phi & \sin\phi \\ 
-\sin\phi & \cos\phi  
\end{pmatrix}
\begin{pmatrix}
H_{N} \\ 
H_{E} 
\end{pmatrix}
$$

Where energy ratios are calculated as:

$$
ETER = \frac{\int T^2(t)dt}{\int R^2(t)dt}
$$

$$
EREZ = \frac{\int R^2(t)dt}{\int Z^2(t)dt}
$$
```


```{seealso}
- Braunmiller, J.,J. Nabelek, and A. Ghods (2020). Sensor Orientation of Iranian Broadband Seismic Stations from P-Wave Particle Motion, Seismol. Res. Lett. 91, 1660–1671, doi: 10.1785/0220200019.

- Pornsopin, P., Pananont, P., Furlong, K.P. et al. Sensor orientation of the TMD seismic network (Thailand) from P-wave particle motions. Geosci. Lett. 10, 24 (2023). https://doi.org/10.1186/s40562-023-00278-7.

- Zhu, G., H. Yang,J. Lin, and Q. You (2020). Determining the Orientation of Ocean-Bottom Seismometers on the Seafloor and Correcting for Polarity Flipping via Polarization Analysis and Waveform Modeling, Seismol. Res. Lett. XX, 1–12, doi: 10.1785/0220190239.
```