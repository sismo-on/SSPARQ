import numpy as np
import pandas as pd
import math
from obspy import read_events

from parameters_py.config import (
					TIME_WINDOW,PERIOD_BANDS_MAX,PERIOD_BANDS_MIN,REMOVE_RESPONSE,OUTPUT_UNIT
				   )

def quakeml_to_dataframe(quakeml_file):
    # Read the QuakeML file using ObsPy
    catalog = read_events(quakeml_file)
    
    entries = []
    
    for event in catalog:
        # Extract basic event information
        origin = event.preferred_origin() or event.origins[0]
        magnitude = event.preferred_magnitude() or event.magnitudes[0]
                
        if event.focal_mechanisms:
            fm = event.focal_mechanisms[0]
            if fm.moment_tensor:
                mt = fm.moment_tensor.tensor
                moment_tensor = [
                    mt.m_rr,
                    mt.m_tt,
                    -mt.m_pp,  # Mff = -Mpp
                    mt.m_rt,
                    -mt.m_rp,  # Mrf = -Mrp
                    -mt.m_tp   # Mtf = -Mtp
                ]

            # Append event data to entries list
            entries.append({
                'time': origin.time.datetime,
                'latitude': origin.latitude,
                'longitude': origin.longitude,
                'depth': origin.depth / 1000,  # Convert from m to km
                'mag': magnitude.mag,
                'magType': magnitude.magnitude_type,
                'moment tensor': moment_tensor})
        else:
            # Append event data to entries list
            entries.append({
                'time': origin.time.datetime,
                'latitude': origin.latitude,
                'longitude': origin.longitude,
                'depth': origin.depth / 1000,  # Convert from m to km
                'mag': magnitude.mag,
                'magType': magnitude.magnitude_type})
        
    # Convert to DataFrame
    df = pd.DataFrame(entries)
    return df


def moment_tensor_to_nodal_planes(input_mt):

    mrr, mtt, mff, mrt, mrf, mtf = input_mt 
    
    """
    Function Name: moment
    Description: Computes scalar seismic moment, compensated linear vector dipole (CLVD) ratio, deviatoric components, isotropic component and its ratio, eigenvectors, and position on the Hudson diagram.
   
    Extracted from: https://github.com/Jose-Alvarez/FMC/blob/master/FMC.py
    Original Author: Jose A. Alvarez-Gomez
    Year: 2015
    """
    
    # Construct the seismic moment tensor (M)
    M = np.array([
        [mrr, mrt, mrf],
        [mrt, mtt, mtf],
        [mrf, mtf, mff]
    ])

    # Remove the isotropic part (mean trace)
    trace = np.trace(M) / 3.0
    M_dev = M - np.eye(3) * trace  # Deviatoric moment tensor

    # Compute eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(M_dev)

    # Sort eigenvalues in ascending order (λ1 < λ2 < λ3)
    idx = np.argsort(eigvals)
    lambda1, lambda2, lambda3 = eigvals[idx]
    v1, v2, v3 = eigvecs[:, idx[0]], eigvecs[:, idx[1]], eigvecs[:, idx[2]]

    # Define principal axes
    P_axis = v1  # Maximum compression (smallest eigenvalue)
    B_axis = v2  # Neutral axis (intermediate eigenvalue)
    T_axis = v3  # Maximum tension (largest eigenvalue)

    # Compute plunges of principal axes   
    plungP = calculate_plunge(P_axis)  # Plunge of P-axis
    plungB = calculate_plunge(B_axis)  # Plunge of B-axis
    plungT = calculate_plunge(T_axis)  # Plunge of T-axis

    # Return plunges
    return (plungP, plungB, plungT)
       
def calculate_plunge(v):
    """Returns the plunge of vector v in degrees"""
    return math.degrees(math.asin(abs(v[2])))

def mecclass(plunges):
    """
    Function Name: mecclass
    Description: Classifies the rupture type of an earthquake based on the plunges of the P, B, and T axes.
    
    Extracted from: https://github.com/Jose-Alvarez/FMC/blob/master/FMC.py
    Original Author: Jose A. Alvarez-Gomez
    Year: 2015
    """

    plunges = np.asarray(plunges)

    P = plunges[0]
    B = plunges[1]
    T = plunges[2]

    maxplung, axis = plunges.max(0), plunges.argmax(0)
    
    if maxplung >= 67.5:
        if axis == 0:  # P max
            clase = 'N'  # normal faulting
        elif axis == 1:  # B max
            clase = 'SS'  # strike-slip faulting
        elif axis == 2:  # T max
            clase = 'R'  # reverse faulting
    else:
        if axis == 0:  # P max
            if B > T:
                clase = 'N-SS'  # normal - strike-slip faulting
            else:
                clase = 'N'  # normal faulting
        if axis == 1:  # B max
            if P > T:
                clase = 'SS-N'  # strike-slip - normal faulting
            else:
                clase = 'SS-R'  # strike-slip - reverse faulting
        if axis == 2:  # T max
            if B > P:
                clase = 'R-SS'  # reverse - strike-slip faulting
            else:
                clase = 'R'  # reverse faulting
    return clase

def adjust_baz_for_ZEN(baz_original):
    """
    Ajusta o BAZ (Back-Azimute) para o sistema ZEN (troca de N e E).
    
    Parâmetros:
    -----------
    baz_original : float
        Back-azimute no sistema ZNE (em graus, 0° a 360°).
    
    Retorna:
    --------
    float
        Novo BAZ no sistema ZEN (em graus, 0° a 360°).
    """
    baz_ZEN = baz_original - 90
    # Ajusta para o intervalo [0°, 360°)
    baz_ZEN = baz_ZEN % 360
    return baz_ZEN

def rms(x):
    """
    Function to calculate root-mean-square of array

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Input array

    Returns
    -------
    rms : float
        Root-Mean-Square value of `x`
    """

    return np.sqrt(np.mean(x**2))

def energy(x):
    """
    Function to calculate energy of array

    Parameters
    ----------
    x : :class:`~numpy.ndarray`
        Input array

    Returns
    -------
    energy : float
        Square value of `x`
    """

    return np.sum(x**2)

def calculate_quartis_mask(orientations):
    # Estimating quartis
    Q1 = np.percentile(orientations, 25)
    Q3 = np.percentile(orientations, 75)
    IQR = Q3 - Q1
                            
    # Defining limits
    lower_bound = Q1 - 1 * IQR
    upper_bound = Q3 + 1 * IQR
    
    # Filter mask
    mask_good = (orientations >= lower_bound) & (orientations <= upper_bound)
    mask_outliers = ~mask_good  # Inverse

    return mask_good,mask_outliers

def format_y_ticks(value, _):
    '''
    Format y-label to degrees
    '''
    return f"{value:.0f}°"
#-------------------------------------------------------------------------------

def preprocess_trace(trace, evtime, inventory, remove_response=REMOVE_RESPONSE,
                     output=OUTPUT_UNIT, pre_filt=(0.005, 0.01, 45.0, 50.0)):
    
    """
    Preprocess a single seismic trace for further analysis.

    This function performs a sequence of standard preprocessing steps on a
    seismic trace, including trimming around the event origin time, removal
    of linear trend and mean (demean), tapering, optional removal of the
    instrumental response using StationXML metadata, and bandpass filtering.

    The instrumental response removal is optional and can be controlled by
    the `remove_response` parameter. When enabled, the trace is converted to
    the desired physical unit (displacement, velocity, or acceleration).

    Parameters
    ----------
    trace : obspy.Trace
        The seismic trace to be processed.
    evtime : obspy.UTCDateTime
        Origin time of the seismic event used to define the trimming window.
    inventory : obspy.Inventory
        Station metadata (StationXML) used to remove the instrument response.
    remove_response : bool, optional
        If True, the instrumental response is removed from the trace.
        Default is False.
    output : str, optional
        Physical unit of the output signal after response removal.
        Accepted values are "DISP", "VEL", or "ACC". Default is "VEL".
    pre_filt : tuple of float, optional
        Pre-filter frequencies (f1, f2, f3, f4) in Hz used during response
        removal to stabilize the deconvolution. Default is (0.005, 0.01, 45, 50).

    Returns
    -------
    trace : obspy.Trace
        The processed seismic trace, ready for analysis (e.g., component
        rotation and orientation optimization).
    """

    # Trim
    trace.trim(evtime - TIME_WINDOW, evtime + TIME_WINDOW)

    # Detrend + demean (ObsPy)
    trace.detrend("linear")
    trace.detrend("demean")

    # Taper
    trace.taper(type="cosine", max_percentage=0.1)

    # Remove instrument response (optional)
    if remove_response:
        trace.remove_response(
            inventory=inventory,
            output=output,
            pre_filt=pre_filt,
            zero_mean=False,
            taper=False,
            plot=False
        )

    # Bandpass
    trace.filter(
        'bandpass',
        freqmin=PERIOD_BANDS_MIN,
        freqmax=PERIOD_BANDS_MAX,
        zerophase=True,
        corners=4
    )

    return trace