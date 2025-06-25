import os
import glob
from datetime import datetime

import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path

import obspy as op
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel
from obspy.signal.headers import clibsignal
from obspy.signal.rotate import rotate_ne_rt
import pyarrow.feather as feather

from parameters_py.config import (
					WAVEFORM_DIR,CATALOG_FILE,XML_DIR,SSPARQ_OUTPUT,num_processes,TAUPY_MODEL,TIME_WINDOW,PERIOD_BANDS_MAX,PERIOD_BANDS_MIN,TIME_FINAL_P
				   )

from src.utils import (
					moment_tensor_to_nodal_planes,calculate_plunge,mecclass,adjust_baz_for_ZEN,rms,energy
				   )
#-------------------------------------------------------------------------------

def aic_simple(a):
    r"""
    Simple Akaike Information Criterion [Maeda1985]_.

    It's computed directly from input data :math:`a` and defined as

    .. math::
        \text{AIC}(k) = k\log(\text{Var}(a_{1..k})) +
                        (N-k-1)\log(\text{Var}(a_{k+1..N}))

    which variance denoted as :math:`\text{Var}`.

    The true output is one data sample less. To make it convenient with other
    metrics in this module, where the output length is preserved, the last
    element is appended to the output: ``aic[-2] == aic[-1]``.

    :type a: :class:`numpy.ndarray` or :class:`list`
    :param a: Input time series
    :rtype: :class:`numpy.ndarray`
    :return: aic - Akaike Information Criterion array

    Extracted from: https://docs.obspy.org/_modules/obspy/signal/trigger.html#aic_simple

    """
    n = len(a)
    if n <= 2:
        return np.zeros(n, dtype=np.float64)
    a = np.ascontiguousarray(a, np.float64)
    aic_res = np.empty(n, dtype=np.float64)
    clibsignal.aic_simple(aic_res, a, n)
    return aic_res
    
#-------------------------------------------------------------------------------

def find_orientation(baz,SS,SZR,ERTR,ERRZ):

    """
    This function calculates the best back azimuth (phi) and sensor misorientation (theta) based on the 
    given quality criteria: signal strength (SS), similarity of vertical and radial components (SZR), 
    transverse-to-radial energy ratio (ERTR), and radial-to-vertical energy ratio (ERRZ).

    The cost function combines these criteria in such a way that minimazing the cost function helps to
    find the optimal back azimuth and corresponding orientation. The function outputs the best back azimuth, orientation,
    and the values of the quality criteria at the best azimuth index.

    Parameters:
    ----------
    baz : float
        Initial back azimuth value from the taup model (degrees).
    SS : np.array
        Array of signal strength values for each azimuth angle.
    SZR : np.array
        Array of similarity values between vertical and radial components for each azimuth angle.
    ERTR : np.array
        Array of transverse-to-radial energy ratios for each azimuth angle.
    ERRZ : np.array
        Array of radial-to-vertical energy ratios for each azimuth angle.

    Returns:
    -------
    phi : float
        The best back azimuth angle (degrees) that minimizes the cost function.
    theta : float
        The sensor misorientation angle (degrees), defined as the difference between the true back azimuth 
        and the estimated back azimuth.
    SS_best : float
        The signal strength value at the best azimuth.
    SZR_best : float
        The similarity between vertical and radial components at the best azimuth.
    ERTR_best : float
        The transverse-to-radial energy ratio at the best azimuth.
    ERRZ_best : float
        The radial-to-vertical energy ratio at the best azimuth.
    """
    
    # Find best index
    cost_function = (
                SS -  # Minimizing energy
                SZR ) # Maximizing similarity
    
    # Best index will minimize the cost function                    
    best_index = np.argmin(cost_function)
    
    # --------------------
    # Search Space of BAZ

    # Step size for the azimuth search (in degrees).
    dphi = 0.1

    # Array of azimuth angles to search through (in degrees).
    ang = np.arange(0., 360., dphi)
                            
    # Get azimuth and correct for angles above 360
    phi = round(ang[best_index])
    theta = round(baz - ang[best_index])

    # Expressed as a deviation from North
    theta = theta % 360          # Convert to (0°, 360°)
    if theta > 180:              
        theta -= 360             # Convert to (-180°, 180°)
                            
    # Get argument of maximum coherence:
    SS_best = SS[best_index]
    SZR_best = SZR[best_index]
    ERTR_best = ERTR[best_index]
    ERRZ_best = ERRZ[best_index]

    return phi,theta,SS_best,SZR_best,ERTR_best,ERRZ_best

# --------------------------------------------------------------------------

def Braunmiller_Pornsopin_algorithm(tr1,tr2,trZ,noise,baz,time_ins,CCVR_MIN=0.5,SNR_MIN=10,TRR_MIN=0.2,RVR_MIN=2):

    """
    Estimate back azimuth using P-wave particle motion and apply quality criteria.

    This algorithm estimates the back azimuth by analyzing P-wave particle motion in an isotropic, 
    homogeneous layered medium. In such a medium, the P-wave energy propagates along a great circle 
    path between the source and receiver, with horizontal components defining the radial direction.
    The angle between the radial direction and true north gives the back azimuth. The P-wave energy 
    is confined to the vertical and radial components, with no energy in the transverse component.

    The sensor 'misorientation angle' is the difference between the true back azimuth (from the taup model) 
    and the empirically estimated back azimuth, with positive values representing a clockwise misorientation.

    This method applies several quality criteria to filter out unreliable results from component malfunctions or 
    missing horizontal components:
    
    == Quality criteria for automatic processing ==

    To select reliable back azimuths in automatic processing, the following five quality criteria are applied:
     - (1) Overall signal strength of the radial component.
     - (2) Similarity between the vertical and radial components.
     - (3) Transverse-to-radial energy ratio.
     - (4) Radial-to-vertical energy ratio.
     - (5) Signal-to-noise ratio (SNR) on the vertical component.

    The function uses these criteria to assess the quality of the estimated back azimuth, and classifies the result 
    as 'good' or 'bad' based on the user-defined thresholds.

    -----------  
    Parameters:
    ----------
    tr1 : np.array
        The first horizontal component of the seismogram.
    tr2 : np.array
        The second horizontal component of the seismogram.
    trZ : np.array
        The vertical component of the seismogram.
    noise : np.array
        Noise window to calculate the signal-to-noise ratio (SNR).
    baz : float
        Back azimuth from the taup model (in degrees).
    time_ins : float
        Difference between observed and predicted travel time
        using available data of distant earthquake (in seconds).
    CCVR_MIN : float, optional
        Minimum required similarity of vertical and radial components (default is 0.45).
    SNR_MIN : float, optional
        Minimum required signal-to-noise ratio (default is 10).
    TRR_MIN : float, optional
        Minimum required transverse-to-radial energy ratio (default is 0.45).
    RVR_MIN : float, optional
        Minimum allowed radial-to-vertical energy ratio (default is -1).
        
    
    Returns:
    -------
    dict
        A dictionary containing the following calculated quality criteria, estimated azimuth, and additional results:
        
        - 'phi' : float
            The estimated back azimuth angle (in degrees) based on the best-fit azimuth search.
        
        - 'baz' : float
            The true back azimuth (in degrees) from the taup model, which serves as a reference for comparison with the estimated azimuth.
        
        - 'SNR' : float
            The signal-to-noise ratio (SNR) of the vertical component, expressed in decibels (dB), representing the strength of the signal relative to noise.
        
        - 'quality' : str
            A classification of the estimated azimuth quality ('good' or 'bad'), based on the comparison of various quality criteria and thresholds.
        
        - 'theta' : float
            The sensor misorientation angle (in degrees), representing the difference between the true back azimuth (from the taup model) and the empirically estimated azimuth.
        
        - 'SS_best' : float
            The best signal strength value for the optimal azimuth, quantifying the energy of the radial component for the best-fit azimuth.
        
        - 'SZR_best' : float
            The best similarity score between the vertical and radial components for the optimal azimuth, indicating how well the vertical and radial components align.
        
        - 'ERTR_best' : float
            The best transverse-to-radial energy ratio for the optimal azimuth, assessing the degree to which the transverse component contaminates the radial component.
        
        - 'ERRZ_best' : float
            The best radial-to-vertical energy ratio for the optimal azimuth, showing the relative strength of the radial component compared to the vertical component.
        
        - 'signal_strength' : np.array
            A NumPy array containing the signal strength values for each azimuth tested in the search range, reflecting the overall energy of the radial component.
        
        - 'similarity_ZR' : np.array
            A NumPy array containing the similarity (correlation) coefficients between the vertical and radial components for each azimuth tested.
        
        - 'energy_ratio_TR' : np.array
            A NumPy array containing the transverse-to-radial energy ratios for each azimuth tested, evaluating the amount of transverse energy relative to radial energy.
        
        - 'max_value_HHR_N' : float
            The gain (amplification factor) of radial maximum amplitude of the North-South (HHN) component.
        
        - 'max_value_HHR_E' : float
            The gain (amplification factor) of radial maximum amplitude of the East-West (HHE) component.
        
    Notes:
    ------
    The algorithm assumes that the back azimuth is between 0 and 360 degrees and that the sensor 
    misorientation is defined as the difference between the true back azimuth and the empirically estimated back azimuth.
    """
    
    # --------------------
    # Search Space of BAZ

    # Step size for the azimuth search (in degrees).
    dphi = 0.1

    # Array of azimuth angles to search through (in degrees).
    ang = np.arange(0., 360., dphi)
    
    # ---------------------------
    # Initialize quality criteria
    
    signal_strength = np.zeros(len(ang))
    similarity_ZR = np.zeros(len(ang))
    energy_ratio_TR = np.zeros(len(ang))
    energy_ratio_RZ = np.zeros(len(ang))
  
    # Search through azimuths and find best-fit azimuth
    for k, an in enumerate(ang):
        R, T = rotate_ne_rt(tr1, tr2, an)

        # (1) Overall signal strength of the transversal component
        signal_strength[k] = energy(T)

        # (2) Similarity of vertical and radial components
        similarity_ZR[k] = np.corrcoef(trZ, R)[0, 1]
        
        # (3) Transverse-to-radial energy ratio
        energy_ratio_TR[k] = energy(T) / energy(R)
        
        # (4) Radial-to-vertical energy ratio
        energy_ratio_RZ[k] = energy(R) / energy(trZ)  
    
    # (5) Signal-to-noise ratio on vertical component
    SNR = round(10.0 * np.log10(rms(trZ)**2 / rms(noise)**2), 1)

    # Normalizing the signal strength of the transversal component
    signal_strength = (signal_strength - np.min(signal_strength)) / (np.max(signal_strength) - np.min(signal_strength)) 

    phi,theta,SS_best,SZR_best,ERTR_best,ERRZ_best = find_orientation(baz,signal_strength,similarity_ZR,energy_ratio_TR,energy_ratio_RZ)

    # Estimating: instrument gain HHN
    new_R_N, new_T_N = rotate_ne_rt(tr1, tr2, phi)

    max_value_HHR_N = np.max(abs(new_R_N))
    
    # Estimating: instrument gain HHE
    new_R_E, new_T_E = rotate_ne_rt(tr2, tr1, adjust_baz_for_ZEN(phi))

    max_value_HHR_E = np.max(abs(new_R_E))

    # Estimating: instrument gain HHZ
    
    max_value_HHZ = np.max(abs(trZ))


    if (SZR_best >= CCVR_MIN) & (SNR >= SNR_MIN) & (ERTR_best < TRR_MIN) & (ERRZ_best >= RVR_MIN) &  (-90 < time_ins < 90):
        quality = 'good'
    else:
        quality = 'bad'

    # Collect results
    results = {
        'phi': phi,
        'baz': baz,
        'SNR': SNR,
        'quality': quality,
        'theta': theta,
        'SS_best': SS_best,
        'SZR_best': SZR_best,
        'ERTR_best': ERTR_best,
        'ERRZ_best': ERRZ_best,
        'signal_strength': signal_strength,
        'similarity_ZR': similarity_ZR,
        'energy_ratio_TR': energy_ratio_TR,
        'energy_ratio_RZ': energy_ratio_RZ,
        'gain_HHN': max_value_HHR_N,
        'gain_HHE': max_value_HHR_E,        
        'gain_HHZ': max_value_HHZ,        
    }
    
    return results

# ---------------------------------------------------------------------------------------------------

def calculate_metrics(input_lst):

    """
    This function calculates the optimal orientation of horizontal seismic components (phi, theta) for a given station, 
    using a set of qualifying seismic events. It processes waveform data for P-wave arrivals and evaluates their quality 
    based on signal-to-noise ratio (SNR), signal strength, and component energy ratios.

    The function reads metadata from a StationXML file and searches for event-specific waveform files in a given directory. 
    For each eligible event, the function estimates the arrival time of the P-wave using the TauP model and applies a 
    modified version of the Braunmiller and Pornsopin algorithm to find the best azimuthal orientation of the sensor.

    The results—including orientation angles, waveform data, energy metrics, and event classification—are stored in 
    a `.feather` file per event for later analysis.

    Parameters:
    ----------
    input_lst : list
        A list containing two elements:
            - XML_FILE (str): Path to the station XML metadata file (StationXML format).
            - WAVE_DIR (str): Path to the directory containing waveform files associated with the events.
            - CAT (list): Catalog with events.
            - SSPARQ_OUTPUT (str): Path to the output folder.

    Returns:
    -------
    None
        This function does not return a value. Instead, it saves a `.feather` file for each valid event to:
        SSPARQ_OUTPUT/FEATHER_FILES/METRICS/{NETWORK}.{STATION}/
        The file contains a DataFrame with waveform segments, results (phi, theta), SNR, and event metadata.
    """

    XML_FILE = input_lst[0]
    WAVE_DIR = input_lst[1]
    CAT = input_lst[2]
    SSPARQ_OUTPUT = input_lst[3]
    
    # ---------------
    # Read XML file
                        
    station_xml = op.read_inventory(XML_FILE)
    network = station_xml[0].code
    station = station_xml[0][0].code    
    
    # ---------------------------
    # Retrieving events waveforms
    # ---------------------------
    
    for evid,event in tqdm(CAT.iterrows(), total=len(CAT),desc=station+' estimation'):
        # ------------------------------
        # Check if the event is eligible

        evdp = event['depth']
        evmag = event['mag']
        evtype = event['magType']
        evtime = UTCDateTime(event['time'])
        evname = UTCDateTime(event['time']).strftime('%Y.%j.%H.%M.%S')

        year = UTCDateTime(event['time']).strftime('%Y')
        julian_day = UTCDateTime(event['time']).strftime('%j')

        # -------------------------------
        # Epicentral distance estimation:
            
        stlo = station_xml[-1][-1][-1].longitude
        stla = station_xml[-1][-1][-1].latitude
            
        evla = event['latitude']
        evlo = event['longitude']
            
        dist,az,baz = gps2dist_azimuth(evla, evlo,stla, stlo)
        gcarc = kilometers2degrees(dist/1000)

        # -------------------------------
        # Taup: theoretical travel times 

        model = TauPyModel(model=TAUPY_MODEL)
        arrivals = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=gcarc,phase_list=['P','PKP','PKIKP'])
        
        if len(arrivals) > 0 and (gcarc < 100 or 140 < gcarc <= 180):
            
            # Event time + first phase arrival time
            evtime = evtime+arrivals[0].time
            
            # ----------------------------
            # Check if feather file exists
            
            output_FEATHER_FILES_METRICS = SSPARQ_OUTPUT+'FEATHER_FILES/METRICS/'+network+'.'+station+'/'
            
            file_feather_name = output_FEATHER_FILES_METRICS+network+'.'+station+'.'+evname+'.metrics.feather'
    
            station_pwd = list(Path(WAVE_DIR).rglob(f'*{evname}*'))
  
            if os.path.isfile(file_feather_name):
                pass
        
            else:
                # -------------------------------
                # Check if components file exists
                        
                files = [str(x) for x in station_pwd if str(x).endswith(evname)]

                if len(files) >= 3:
                        
                    file_HHE = [x for x in files if "HE." in x or "H2." in x][0]
                    file_HHN = [x for x in files if "HN." in x or "H1." in x][0]
                    file_HHZ = [x for x in files if "HZ." in x][0]
            
                    # --------
                    # Data HHE
                                
                    tr2_data_file = op.read(file_HHE)
                    tr2_data_file.trim(evtime-TIME_WINDOW,evtime+TIME_WINDOW)
                    tr2_data_file.taper(type='cosine',max_percentage=0.1)
                    tr2_data_file.filter('bandpass',freqmin=PERIOD_BANDS_MIN,freqmax=PERIOD_BANDS_MAX,zerophase=True, corners=4)
                    
                    # --------
                    # Data HHN
                    
                    tr1_data_file = op.read(file_HHN)
                    tr1_data_file.trim(evtime-TIME_WINDOW,evtime+TIME_WINDOW)
                    tr1_data_file.taper(type='cosine',max_percentage=0.1)
                    tr1_data_file.filter('bandpass',freqmin=PERIOD_BANDS_MIN,freqmax=PERIOD_BANDS_MAX,zerophase=True, corners=4)
                            
                    # --------
                    # Data HHZ
                                
                    trZ_data_file = op.read(file_HHZ)
                    trZ_data_file.trim(evtime-TIME_WINDOW,evtime+TIME_WINDOW)
                    trZ_data_file.taper(type='cosine',max_percentage=0.1)
                    trZ_data_file.filter('bandpass',freqmin=PERIOD_BANDS_MIN,freqmax=PERIOD_BANDS_MAX,zerophase=True, corners=4)

                    if len(tr2_data_file) > 0 and len(tr1_data_file) > 0 and len(trZ_data_file) > 0:
                               
                            if (tr2_data_file[0].stats.npts == tr1_data_file[0].stats.npts == trZ_data_file[0].stats.npts) and (trZ_data_file[0].stats.npts > TIME_WINDOW * tr2_data_file[0].stats.sampling_rate):

                                # -------------------------------------------------------------------------------------
                                # Remove 5 seconds from the beginning and end of the waveform to eliminate edge effects
                                
                                # HHE
                                tr2_data_filtered = tr2_data_file[0].data[500:-500]
                                
                                # HHN
                                tr1_data_filtered = tr1_data_file[0].data[500:-500]
                                    
                                # HHZ
                                trZ_data_filtered = trZ_data_file[0].data[500:-500]
                                trZ_time = trZ_data_file[0].times()[500:-500]-TIME_WINDOW

                                # -------------------------------------------
                                # Function estimates the Akaike Information directly from data 
                                # The Summed Log Likelihood section implies that a natural 
                                # changepoint estimate is the sample index that minimizes 
                                # the AIC in equation
                                
                                aic_curve = aic_simple(trZ_data_filtered)
                                
                                k_min_index = np.argmin(aic_curve)
                            
                                time_P_arr = trZ_time[k_min_index]
                                                        
                                # Time error:
                                time_ins = round(time_P_arr,1)
                            
                                if time_ins > 0:
                                    signal_window_start = time_ins
                                    signal_window_final = time_ins+TIME_FINAL_P
                                    noise_window_start = time_ins
                                    noise_window_final = time_ins-TIME_FINAL_P
                                else:
                                    signal_window_start = time_ins
                                    signal_window_final = TIME_FINAL_P+time_ins
                                    noise_window_start = time_ins
                                    noise_window_final = -(abs(time_ins)+TIME_FINAL_P) 
                                # -------------------------------------------------------------------------------------------------------------------------------
                                # Signal and noise windows
                                        
                                signal_window = (trZ_time >= signal_window_start) & (trZ_time <= signal_window_final)
                                noise_window = (trZ_time >= noise_window_final) & (trZ_time <= noise_window_start)
                                
                                noise = trZ_data_filtered[noise_window]
                                trZ_noise_time = trZ_time[noise_window]
            
                                tr2 = tr2_data_filtered[signal_window]
                                tr1 = tr1_data_filtered[signal_window]
                                trZ = trZ_data_filtered[signal_window]
                                trZ_signal_time = trZ_time[signal_window]
                                    
                                # -------------------------------------------------------------------------------------------------------------------------------
                                # Calculating the optimal orientation
                                                    
                                results = Braunmiller_Pornsopin_algorithm(tr1,tr2,trZ,noise,baz,time_ins,CCVR_MIN=0.45,SNR_MIN=10,TRR_MIN=0.45,RVR_MIN=-1)

                                # -------------------------------------------------------------------------------------------------------------------------------
                                # Calculating the Plunge of: P, B, and T axis

                                if not event.get('moment tensor'):
                                    # ----------------------------------------------------------------------------------------------------                               
                                    # Creating a Pandas DataFrame:
                                    column_info = [network,station,stla,stlo,evname,evla,evlo,evtime,evmag,evtype,evdp,dist,gcarc,baz,tr1_data_filtered,tr2_data_filtered,trZ_data_filtered,trZ_time,results['SS_best'],results['signal_strength'],results['SZR_best'],results['similarity_ZR'],results['ERTR_best'],results['energy_ratio_TR'],results['ERRZ_best'],results['energy_ratio_RZ'],results['SNR'],results['phi'],results['theta'],aic_curve,time_ins,results['quality'],results['gain_HHN'],results['gain_HHE'],results['gain_HHZ']]
                                    columns_header = ['network','station','stla','stlo','evname','evla','evlo','evtime','evmag','evtype','evdp','distance','gcarc','baz','tr1_data','tr2_data','trZ_data','trZ_time','SS_best','signal_strength','SZR_best','similarity_vertical_radial','ERTR_best','energy_transverse_radial','ERRZ_best','energy_radial_vertical','SNR','phi','theta','aic_curve','clock_error','quality','gain_HHN','gain_HHE','gain_HHZ']
                                    
                                else:                                
                                    nodal_planes = moment_tensor_to_nodal_planes(event['moment tensor'])
                                    event_class = mecclass(nodal_planes)

                                    # ----------------------------------------------------------------------------------------------------                               
                                    # Creating a Pandas DataFrame:
                                    column_info = [network,station,stla,stlo,evname,evla,evlo,evtime,evmag,evtype,evdp,dist,gcarc,baz,tr1_data_filtered,tr2_data_filtered,trZ_data_filtered,trZ_time,results['SS_best'],results['signal_strength'],results['SZR_best'],results['similarity_ZR'],results['ERTR_best'],results['energy_ratio_TR'],results['ERRZ_best'],results['energy_ratio_RZ'],results['SNR'],results['phi'],results['theta'],aic_curve,time_ins,results['quality'],results['gain_HHN'],results['gain_HHE'],results['gain_HHZ'],event['moment tensor'],nodal_planes,event_class]
                                    columns_header = ['network','station','stla','stlo','evname','evla','evlo','evtime','evmag','evtype','evdp','distance','gcarc','baz','tr1_data','tr2_data','trZ_data','trZ_time','SS_best','signal_strength','SZR_best','similarity_vertical_radial','ERTR_best','energy_transverse_radial','ERRZ_best','energy_radial_vertical','SNR','phi','theta','aic_curve','clock_error','quality','gain_HHN','gain_HHE','gain_HHZ','moment tensor','nodal_planes','event_class']
                                    
                                    
                                metrics_p_wave_df = pd.DataFrame(column_info, index=columns_header).T
                                metrics_p_wave_df['evtime'] = pd.to_datetime(metrics_p_wave_df['evtime'].apply(lambda x: x.isoformat() if isinstance(x, UTCDateTime) else x))

                                # ----------------------------------------------------------------------------------------------------
                                # Convert from pandas to Arrow and saving in feather formart file
                                os.makedirs(output_FEATHER_FILES_METRICS,exist_ok=True)
                                feather.write_feather(metrics_p_wave_df, file_feather_name)