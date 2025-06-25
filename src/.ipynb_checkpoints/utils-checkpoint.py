import numpy as np
import pandas as pd
import math
from obspy import read_events

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
#-------------------------------------------------------------------------------