# Importing modules

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,ScalarFormatter,FuncFormatter
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects


import cartopy.crs as ccrs
import cartopy.feature as cfeature

from obspy.signal.rotate import rotate_ne_rt
from obspy.taup import TauPyModel
from obspy.imaging.beachball import beach

import glob
import os
import numpy as np
import pandas as pd

from scipy.stats import circmean, circstd,gaussian_kde
from kneed import KneeLocator
import pyarrow.feather as feather

from sklearn.cluster import DBSCAN
from sklearn.preprocessing import RobustScaler
from sklearn.metrics import silhouette_score


# Importing SSPARQ setup:
from parameters_py.config import (
					WAVEFORM_DIR,CATALOG_FILE,XML_DIR,SSPARQ_OUTPUT,num_processes,TIME_FINAL_P,TIME_WINDOW
				   )

# Importing SSPARQ functions:
from src.utils import (
					quakeml_to_dataframe,moment_tensor_to_nodal_planes,calculate_plunge,mecclass,adjust_baz_for_ZEN,rms,energy,
                    format_y_ticks,calculate_quartis_mask
				   )

# Sets the global font size
plt.rcParams.update({'font.size': 14}) 

# --------- #
# Functions #
# --------- #


def plotting_event_orientation(df_row,SSPARQ_OUTPUT=SSPARQ_OUTPUT,TIME_FINAL_P=TIME_FINAL_P,TIME_WINDOW=TIME_WINDOW):

    """
    Plot a comprehensive overview of a seismic event and station orientation analysis.

    This function generates a multi-panel figure visualizing different aspects of a seismic event 
    recorded at a station, including waveform components (vertical, radial, transverse), signal 
    and noise windows, AIC changepoint detection, orientation metrics (e.g., signal strength, 
    energy ratios), a global map showing station and event locations, a beachball diagram for the 
    focal mechanism (if available), and ray path projections. The final figure is saved to disk.

    Parameters
    ----------
    df_row : pd.Series
        A row of a DataFrame containing pre-processed waveform data, metadata, and orientation analysis results.
    SSPARQ_OUTPUT : str, optional
        Path to the base output directory where figures will be saved.
    TIME_FINAL_P : float, optional
        Time window (in seconds) after the P-arrival used for signal/noise segmentation.
    TIME_WINDOW : float, optional
        Full time window (in seconds) around the origin for plotting waveform data.

    Notes
    -----
    - Assumes that waveform data (`trZ_data`, `tr1_data`, `tr2_data`) and timing (`trZ_time`) 
      are included in `df_row`.
    - Rotation from NE to RT is performed using the `phi` angle in `df_row`.
    - The moment tensor is used to draw the beachball plot if present.
    - Ray paths are computed using the TauPyModel with the IASP91 velocity model.

    Saves
    -----
    A PNG file is created in the path:
        SSPARQ_OUTPUT/FIGURES/EARTHQUAKES/<network>.<station>/METRICS_<station>_<event>_<quality>.png
    """
    
    df_row = df_row[1]

    # -------------------------------------------
    # Function estimates the Akaike Information directly from data 
    # For a given time series
    # The Summed Log Likelihood section implies that a natural 
    # changepoint estimate is the sample index that minimizes 
    # the AIC in equation
    
    aic_curve = df_row['aic_curve']

    k_min_index = np.argmin(aic_curve)
                            
    amp_P_arr = aic_curve[k_min_index]
    time_P_arr = df_row['clock_error']
    
    # Time instability:
    time_ins = df_row['clock_error']

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
            
    signal_window = (df_row['trZ_time'] >= signal_window_start) & (df_row['trZ_time'] <= signal_window_final)
    noise_window = (df_row['trZ_time'] >= noise_window_final) & (df_row['trZ_time'] <= noise_window_start)
                                    
    noise = df_row['trZ_data'][noise_window]
    trZ_noise_time = df_row['trZ_time'][noise_window]
            
    tr2 = df_row['tr2_data'][signal_window]
    tr1 = df_row['tr1_data'][signal_window]
    trZ = df_row['trZ_data'][signal_window]
    trZ_signal_time = df_row['trZ_time'][signal_window]
                
    # --------------------
    # figure 
    fig = plt.figure(figsize=(15, 7),constrained_layout=True)

    qc_symbol = '[✔]' if df_row['quality'] == 'good' else '[✘]'

    if not df_row.get('event_class'):

        fig.suptitle(
        f"{qc_symbol} Evento: {df_row['evname']} "
        f"(Δ: {round(df_row['gcarc'])}° | M: {round(df_row['evmag'],1)} {df_row['evtype']} | "
        f"D: {round(df_row['evdp'])} km) \n "
        f"SNR: {df_row['SNR']} dB | BAZ: {round(df_row['baz'])}° | "
        f"PHI: {df_row['phi']}° | theta: {df_row['theta']}° | "
        f"TI: {df_row['clock_error']}s | "
        f"GN: {df_row['gain_HHN']:.2e} | GE: {df_row['gain_HHE']:.2e}  | GZ: {df_row['gain_HHZ']:.2e}", 
        fontsize=15)    
    else:
        fig.suptitle(
        f"{qc_symbol} Evento: {df_row['evname']} "
        f"(Δ: {round(df_row['gcarc'])}° | M: {round(df_row['evmag'],1)} {df_row['evtype']} | "
        f"D: {round(df_row['evdp'])} km | C: {df_row['event_class']}) \n "
        f"SNR: {df_row['SNR']} dB | BAZ: {round(df_row['baz'])}° | "
        f"PHI: {df_row['phi']}° | theta: {df_row['theta']}° | "
        f"TI: {df_row['clock_error']}s | "
        f"GN: {df_row['gain_HHN']:.2e} | GE: {df_row['gain_HHE']:.2e}  | GZ: {df_row['gain_HHZ']:.2e}", 
        fontsize=15)  
        
    # creating grid
    gs = fig.add_gridspec(1, 2,width_ratios=[5,1])

    gs0 = gs[0].subgridspec(4, 1,height_ratios=[5,5,5,1])
    gs1 = gs[1].subgridspec(4, 1)
                                
    # Rotating components
    new_R, new_T = rotate_ne_rt(df_row['tr1_data'], df_row['tr2_data'], df_row['phi'])

    # Transversal data
    ax1 = fig.add_subplot(gs0[0, 0])
    ax1.plot(df_row['trZ_time'],new_T,'-k',lw=2,label='HHT')
    ax1.plot(trZ_signal_time,tr2,c='gray',ls='--',lw=1,label='HHE')
    ax1.axvspan(xmin=signal_window_start, xmax=signal_window_final, ymin=0, ymax=1,facecolor='none', edgecolor='blue', linestyle='--', lw=2,alpha=0.25)
    ax1.annotate(df_row['network']+'.'+df_row['station']+'..HHT', (0.95, 0.85),xycoords='axes fraction',fontsize=15, va='center',ha='right',bbox=dict(boxstyle="round", fc="white"))
    ax1.set_xlim(-TIME_WINDOW,TIME_WINDOW)
    ax1.tick_params(axis="x", labelbottom=False)
    ax1.grid(which='major',linestyle=':')
    ax1.legend(loc='lower left')

    # Radial data
    ax2 = fig.add_subplot(gs0[1, 0], sharex=ax1, sharey=ax1)
    ax2.plot(df_row['trZ_time'],new_R,'-k',label='HHR')
    ax2.plot(trZ_signal_time,tr1,c='gray',ls='--',lw=1,label='HHN')
    ax2.axvspan(xmin=signal_window_start, xmax=signal_window_final, ymin=0, ymax=1,facecolor='none', edgecolor='blue', linestyle='--', lw=2,alpha=0.25)
    ax2.annotate(df_row['network']+'.'+df_row['station']+'..HHR', (0.95, 0.85),xycoords='axes fraction',fontsize=15, va='center',ha='right',bbox=dict(boxstyle="round", fc="white"))
    ax2.tick_params(axis="x", labelbottom=False)
    ax2.grid(which='major',linestyle=':')
    ax2.legend(loc='lower left')

    # Vertical data and noise and signal window
    ax3 = fig.add_subplot(gs0[2, 0], sharex=ax1, sharey=ax1)
    ax3.plot(df_row['trZ_time'],df_row['trZ_data'],'-k')
    ax3.tick_params(axis="x", labelbottom=False)
    ax3.annotate(df_row['network']+'.'+df_row['station']+'..HHZ', (0.95, 0.85),xycoords='axes fraction',fontsize=15, va='center',ha='right',bbox=dict(boxstyle="round", fc="white"))
    ax3.axvspan(xmin=signal_window_start, xmax=signal_window_final, ymin=0, ymax=1,facecolor='none', edgecolor='blue', linestyle='--', lw=2,alpha=0.25,label='signal')
    ax3.axvspan(xmin=noise_window_start, xmax=noise_window_final, ymin=0, ymax=1,facecolor='none', edgecolor='red', linestyle='--', lw=2,alpha=0.25,label='noise')
    ax3.grid(which='major',linestyle=':')
    ax3.legend(loc='lower left')

    # Vertical data AIC curve
    ax4 = fig.add_subplot(gs0[3, 0], sharex=ax1)
    ax4.plot(df_row['trZ_time'],aic_curve, 'k')
    ax4.scatter(time_P_arr,amp_P_arr, marker='v',s=100,c='k',lw=1,ec='r',label='changepoint')
    ax4.annotate('AIC', (0.95, 0.85),xycoords='axes fraction',fontsize=15, va='center',ha='right',bbox=dict(boxstyle="round", fc="white"))
    ax4.set_xlabel('Timelag (s)',fontsize=15)
    ax4.grid(which='major',linestyle=':')                          
    ax4.legend(loc='lower left')
    ax4.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax4.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    
    # --------------------
    # Search Space of BAZ

    # Step size for the azimuth search (in degrees).
    dphi = 0.1

    # Array of azimuth angles to search through (in degrees).
    ang = np.arange(0., 360., dphi)
    
    # --------------------  
    # Transversal signal strength
    
    ax5 = fig.add_subplot(gs1[0, 0])
    ax5.plot(ang,df_row['signal_strength'],'.k')
    ax5.plot(df_row['phi'],df_row['SS_best'],'*r',ms=10)
    ax5.set_ylim(0,1)
    ax5.set_xlim(0,360)
    ax5.tick_params(axis="both", which='both',labelbottom=False, labelright=True, labelleft=True, labeltop=False,bottom=True, top=True, left=True, right=True)
    ax5.set_title(r'$SS_{T}$',fontsize=15)
    ax5.grid(which='major',linestyle=':')                          
                                
    # Similarity between vertical and radial
    ax6 = fig.add_subplot(gs1[1, 0],sharex=ax5)
    ax6.plot(ang,df_row['similarity_vertical_radial'],'.k')
    ax6.plot(df_row['phi'],df_row['SZR_best'],'*r',ms=10)
    ax6.set_ylim(-1,1)
    ax6.tick_params(axis="both", which='both',labelbottom=False, labelright=True, labelleft=True, labeltop=False,bottom=True, top=True, left=True, right=True)
    ax6.set_title(r'$CC_{RZ}$',fontsize=15)
    ax6.grid(which='major',linestyle=':')                          

    # Transverse-to-Radial Energy Ratio
    ax7 = fig.add_subplot(gs1[2, 0],sharex=ax5)
    ax7.plot(ang,df_row['energy_transverse_radial'],'.',c='k')
    ax7.plot(df_row['phi'],df_row['ERTR_best'],'*r',ms=10)
    ax7.tick_params(axis="both", which='both',labelbottom=False, labelright=True, labelleft=True, labeltop=False,bottom=True, top=True, left=True, right=True)
    ax7.set_title(r'$E_{T}/E_{R}$',fontsize=15)
    ax7.grid(which='major',linestyle=':')                          
    ax7.yaxis.label.set_color('k')
                                                                
    # Radial-to-Vertical Energy Ratio
    ax8 = fig.add_subplot(gs1[3, 0],sharex=ax5)
    ax8.plot(ang,df_row['energy_radial_vertical'],'.',c='k')
    ax8.plot(df_row['phi'],df_row['ERRZ_best'],'*r',ms=10)
    ax8.set_title(r'$E_{R}/E_{Z}$',fontsize=15)
    ax8.set_xlabel('Orientation Angle (deg)',fontsize=15)
    ax8.grid(which='major',linestyle=':')                          
    ax8.tick_params(axis="both", which='both',labelbottom=True, labelright=True, labelleft=True, labeltop=False,bottom=True, top=True, left=True, right=True)
    ax8.yaxis.label.set_color('k')

    # --------------------------
    # Adding global location map

    ax_map = plt.axes([0.025, 0.76, 0.12, 0.12], projection=ccrs.Orthographic(central_latitude=df_row['stla'],central_longitude=df_row['stlo']))
    ax_map.set_global()

    # ---------------------
    # Adding background map 

    ax_map.add_feature(cfeature.LAND)
    ax_map.add_feature(cfeature.OCEAN)
    ax_map.add_feature(cfeature.COASTLINE)
                            
    ax_map.scatter(df_row['evlo'],df_row['evla'],color="y",marker='*',s=200,ec='k',transform=ccrs.PlateCarree())
    ax_map.scatter(df_row['stlo'],df_row['stla'],color="r",marker='^',s=50,transform=ccrs.PlateCarree())
    ax_map.plot([df_row['stlo'], df_row['evlo']], [df_row['stla'], df_row['evla']], c='gray',ls='-',lw=2, transform=ccrs.Geodetic())

    mt = df_row.get('moment tensor')

    if mt is not None and len(mt):
    # continue com o plot
    #if df_row.get('moment tensor'):

        # ===================================================================================
        # focal mechanism (https://docs.obspy.org/tutorial/code_snippets/beachball_plot.html)
        # ===================================================================================
                                        
        # ---------------------------
        # Plotting: adding a new axes
                                        
        newax = fig.add_axes([-0.06, 0.34, 0.3,  0.3])
        
        # -----------------------------------------
        # Computing: Retrieving moment tensor info
                                                
        mt = df_row['moment tensor'] 
        
        # -----------------------------------------------------------------------------------------------------------------------------------------
        # Plotting: graphical representation of a focal mechanism (https://docs.obspy.org/packages/autogen/obspy.imaging.beachball.beachball.html)
    
        # Normalize event depth values between 0 and 600 km:
        min_val = 0
        max_val = 600
        normalized_values = [(x - min_val) / (max_val - min_val) for x in np.arange(min_val, max_val,10)]
    
        # Colormap "Plasma" for each value
        colors = [plt.cm.Spectral(value) for value in normalized_values]
                                                                            
        # Convert colors RGB to hexadecimal:
        hex_colors = [mcolors.rgb2hex(color) for color in colors]
    
        # Find the color for a given depth
        diff_ev_depth = [np.abs(numero - df_row['evdp']) for numero in np.arange(min_val, max_val,10)]
                                            
        # Find the min index for a given depth
        index_min_ev_depth = diff_ev_depth.index(min(diff_ev_depth))
    
        # Plotting the hexcolor
        bball = beach(fm=mt, xy=(0, 0.5),size=200, width=0.75, facecolor=hex_colors[index_min_ev_depth])
            
        # -------------------------
        # Plotting: axes parameters 
                                    
        newax.add_collection(bball)
        newax.set_xlim(-1, 1)
        newax.set_ylim(-1, 1)
        newax.set_aspect('equal')
        newax.axis('off')
            
    # ===========================================================
    # ray paths (https://docs.obspy.org/packages/obspy.taup.html)
    # ===========================================================

    # ---------------------------------------------------------------------------------------------------
    # Computing: The paths travelled by the rays to the receiver for a given phase and 1D velocity model 
    
    model = TauPyModel(model="iasp91")
    arrivals_ray_path = model.get_ray_paths(source_depth_in_km=df_row['evdp'], distance_in_degree=df_row['gcarc'], phase_list=['P','PKP','PKIKP'])

    # -------------------------
    # Plotting: axes parameters 
                                    
    ax_raypath = fig.add_axes([0.022, 0.24, 0.13,  0.13], projection='polar')
    arrivals_ray_path.plot_rays(ax=ax_raypath)
        
    # ==========================================
    
    output_figure_SSPARQ = SSPARQ_OUTPUT+'FIGURES/EARTHQUAKES/'+df_row['network']+'.'+df_row['station']+'/'
    os.makedirs(output_figure_SSPARQ,exist_ok=True)
    fig.savefig(output_figure_SSPARQ+'METRICS_'+df_row['station']+'_'+df_row['evname']+'_'+df_row['quality']+'.png',dpi=100)
    plt.close()

# ---------------------------------------------------------------------------------------------------

def station_overview_metrics(net_sta,SSPARQ_OUTPUT=SSPARQ_OUTPUT):

    """
    Generate comprehensive station metrics visualization including orientation, clock error, 
    and gain analysis over time using DBSCAN clustering or quartile-based outlier detection.

    This function processes seismic station data to create a multi-panel visualization showing:
    - Event counts and orientation angles over time
    - Clock error measurements
    - Gain ratios between components
    - Kernel density estimation of orientations
    
    The analysis uses DBSCAN clustering only when the silhouette score (measuring cluster separation) exceeds 0.22; 
    otherwise, it defaults to quartile-based outlier detection.

    Parameters:
    -----------
    net_sta : tuple
        A tuple containing (network_code, station_code) to identify the station.
    SSPARQ_OUTPUT : str, optional
        Path to the SSPARQ output directory containing the metrics feather files.
        Default is the global SSPARQ_OUTPUT variable.

    Returns:
    --------
    None
    
    Outputs:
    --------
    Saves a PNG visualization to:
    {SSPARQ_OUTPUT}/FIGURES/FINAL_RESULT/{network_code}/METRICS_TOTAL_{network_code}_{station_code}.png

    The visualization includes:
    - Top panel: Event count histogram by month (colored by quality/cluster)
    - Middle panel: Orientation angle scatter plot with SNR scaling
    - Lower middle panel: Clock error measurements
    - Bottom panel: Gain ratio comparisons (E/N, E/Z, N/Z)
    - Right panel: Kernel density estimate of orientation angles

    Notes:
    ------
    - The function skips stations that already have output figures
    - DBSCAN parameters are automatically tuned based on data characteristics
    - Color coding distinguishes between:
        * Good quality clusters (spectral colors)
        * Outliers (purple)
        * Bad quality data (gray)
    - The figure includes statistical annotations (mean±std) for orientation clusters
    """
   
    net = net_sta[0]
    sta = net_sta[1]

    colnames = ['network', 'station','evtime','SNR', 'phi', 'theta','clock_error', 'quality', 'gain_HHN', 'gain_HHE', 'gain_HHZ','event_class']
    
    feather_files_lst = [pd.read_feather(i,columns=colnames) for i in glob.glob(SSPARQ_OUTPUT+'FEATHER_FILES/METRICS/'+net+'.'+sta+'/*')]

    # Check if figure exists:
    
    output_figure_SSPARQ = SSPARQ_OUTPUT + 'FIGURES/FINAL_RESULT/'+net+'/'
    
    if len(feather_files_lst) > 1 and os.path.isfile(output_figure_SSPARQ + f'METRICS_TOTAL_{net}_{sta}.png') == False:
    
        station_df = pd.concat(feather_files_lst)
        station_df['year_month'] = station_df['evtime'].dt.to_period('M').astype(str)
        station_df['year_month'] = pd.to_datetime(station_df['year_month'], format='%Y-%m').dt.to_period('M')
        
        df_sta = station_df[station_df['station'] == sta].copy()
    
        # ------------------------------------
        # Data
        
        orientations_all_good = df_sta[df_sta['quality'] == 'good']['theta'].values # Orientations
        time_all_good = df_sta[df_sta['quality'] == 'good']['evtime'].values  # Time in datetime
        time_all_stamp = df_sta[df_sta['quality'] == 'good']['evtime'].apply(lambda x: int(x.timestamp()))  # Time in Timestamp
    
        if len(orientations_all_good) > 50:
    
            # ================================= #
            # START: DBSCAN clusters estimation #
            # ================================= #
            
            # ------------------------------------
            # This Scaler removes the median and scales the data according to the quantile range (defaults to IQR: Interquartile Range). 
            # The IQR is the range between the 1st quartile (25th quantile) and the 3rd quartile (75th quantile).
        
            data_all = np.array([time_all_stamp,orientations_all_good]).T
            scaler = RobustScaler()
            scaler.fit(data_all)
            data_scale = scaler.transform(data_all)
    
            # The size of the radius is specified by the distance threshold parameter (epsilon)
            
            eps_range = np.arange(0.23,0.28,0.005)
            per_samples = 20 # percentage of total de samples per group
    
            # -----------------------------
            # Silhouette DBSCAN estimation
            
            eps_lst = []
            silhouette_score_lst = []
    
            for n_eps in eps_range:
                try:
                    clustering = DBSCAN(eps=round(n_eps,3), min_samples=len(orientations_all_good)//per_samples).fit(data_scale)
                    ss = silhouette_score(data_scale, clustering.fit_predict(data_scale))
                    
                    eps_lst.append(n_eps)
                    silhouette_score_lst.append(ss)
                except:
                    pass 
    
            if len(silhouette_score_lst) > 0:
                # ------------------------------------    
                # Elbow DBSCAN estimation
                
                kn = KneeLocator(eps_lst, silhouette_score_lst, curve='concave',polynomial_degree=3,S=5)
                elbow_point = kn.knee
                
                if not elbow_point:
                    kn = KneeLocator(eps_lst, silhouette_score_lst, curve='concave',polynomial_degree=3)
                    elbow_point = kn.knee
            
                silhouette_score_elbow_point = silhouette_score_lst[eps_lst.index(elbow_point)]
    
    
                # ------------------------
                # DBSCAN clustering result
    
                clustering = DBSCAN(eps=round(elbow_point,2), min_samples=len(orientations_all_good)//per_samples).fit(data_scale)
                
                df_sta['class'] = -10  # standart value
                df_sta.loc[df_sta['quality'] == 'good', 'class'] = clustering.labels_
                unique_labels = list(set(clustering.labels_))
    
                if silhouette_score_elbow_point > 0.22 and len(unique_labels) < 6:
    
                    colors = [mcolors.to_rgba(plt.cm.Spectral(each),alpha=0.5) for each in np.linspace(0, 1, len(unique_labels))]
                    label_dbscan = ['g'+str(1+w)+': ' if w != -1 else 'ol: ' for w in unique_labels]
                    
                    label_dbscan_lst_end = [] 
                    for idx,labe in enumerate(unique_labels):
                        label_dbscan_lst_end.append(label_dbscan[idx]+str(len(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == labe)]['theta'].values)))
                
                    # =============================== #
                    # END: DBSCAN clusters estimation #
                    # =============================== #
                
                    # Cria o range de meses como Periods e converte para datetime (timestamp)
                    years = pd.period_range(start=df_sta['year_month'].min(), end=df_sta['year_month'].max(), freq='M')
                    
                    fig = plt.figure(figsize=(10, 10))
                    gs = gridspec.GridSpec(4, 2, width_ratios=[10,1], height_ratios=[1,10,1,1],hspace=0.05, wspace=0.05)
                
                    # Definindo os eixos
                    ax0 = fig.add_subplot(gs[0, 0])  #  Number of events axis
                    ax1 = fig.add_subplot(gs[1, 0],sharex=ax0)  # Orientation axis
                    ax2 = fig.add_subplot(gs[2, 0],sharex=ax0)  # Histogram axis
                    ax3 = fig.add_subplot(gs[3, 0],sharex=ax0)  # Time axis
                    ax4 = fig.add_subplot(gs[1, 1],sharey=ax1)  # KDE axis
    
                    if silhouette_score_elbow_point:
                        ax0.annotate('ss:'+str(round(silhouette_score_elbow_point,2)), (pd.to_datetime(df_sta['evtime'].values).max(), +15),fontsize=10, va='center', ha='center',bbox=dict(boxstyle="round", fc="white", ec='k', alpha=0.5))
    
                    for ye in years:
                        # Filtering according to Period (year_month)
                        df_sta_year = df_sta[df_sta['year_month'] == ye]
                
                        ye_num = mdates.date2num(ye)  # converte timestamp para número
    
                        # ---------------- #
                        # Orientation plot #
                        # ---------------- #
                        
                        if df_sta_year[df_sta_year['quality'] == 'good']['theta'].empty:
                            
                            # BAD VALUES #
    
                            orientations_bad = df_sta_year[df_sta_year['quality'] == 'bad']['theta'].values
                            snr_bad = df_sta_year[df_sta_year['quality'] == 'bad']['SNR'].abs().values
    
                            ax0.bar(ye_num, len(orientations_bad), color='gray', width=20, alpha=0.25,zorder=-10)
                            ax1.scatter([ye_num]*len(orientations_bad), orientations_bad, marker='.', c='gray',s=snr_bad*10, alpha=0.05, ec='k', label='bad')
                            
                            # time plot #
    
                            clock_error_bad = df_sta_year[df_sta_year['quality'] == 'bad']['clock_error'].values # Clock instability
    
                            ax2.scatter([ye_num]*len(clock_error_bad), clock_error_bad,color='gray',marker='.',edgecolor='none',alpha=0.1)
    
                            # gain plot #
    
                            gain_HHE_bad = df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHE'].values # Gain HHE bad
                            gain_HHN_bad= df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHN'].values # Gain HHN bad
                            gain_HHZ_bad = df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHZ'].values # Gain HHZ bad
        
                            g1 = ax3.scatter([ye_num]*len(gain_HHZ_bad),np.log(gain_HHE_bad/gain_HHN_bad),c='gray',marker='p',edgecolor='none',s=10,alpha=0.25,label='E/N')
                            g2 = ax3.scatter([ye_num]*len(gain_HHZ_bad),np.log(gain_HHE_bad/gain_HHZ_bad),c='gray',marker='>',edgecolor='none',s=10,alpha=0.25,label='E/Z')
                            g3 = ax3.scatter([ye_num]*len(gain_HHZ_bad),np.log(gain_HHN_bad/gain_HHZ_bad),c='gray',marker='^',edgecolor='none',s=10,alpha=0.25,label='N/Z')
                        
                        else:
                        
                            # GOOD VALUES #
                            
                            orientations_bad = df_sta_year[df_sta_year['quality'] == 'bad']['theta'].values  # Orientations bad
    
                            snr_bad = df_sta_year[df_sta_year['quality'] == 'bad']['SNR'].abs().values # SNR bad
    
                            ax0.bar(ye_num, len(orientations_bad), color='gray', width=20, alpha=0.25,zorder=-9) # Plot histogram bad
                            ax1.scatter([ye_num]*len(orientations_bad), orientations_bad, marker='.', c='gray',s=snr_bad*10, alpha=0.05, ec='k', label='bad') # Plot orientations bad
    
                            clock_error_bad = df_sta_year[df_sta_year['quality'] == 'bad']['clock_error'].values# Clock instability bad
                            ax2.scatter([ye_num]*len(clock_error_bad), clock_error_bad,c='gray',marker='.',edgecolor='none',alpha=0.1) # Plot clock bad
    
                            for uni,col in zip(unique_labels,colors):
                                orientations_good_cluster = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == uni)]['theta'].values # Orientations cluster good
                                
                                snr_good_cluster = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == uni)]['SNR'].abs().values # SNR cluster good
    
                                clock_error_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == uni)]['clock_error'].values # Clock instability cluster good
    
                                gain_HHE_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == uni)]['gain_HHE'].values # Gain HHE cluster good
                                gain_HHN_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == uni)]['gain_HHN'].values # Gain HHN cluster good
                                gain_HHZ_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == uni)]['gain_HHZ'].values # Gain HHZ cluster good
    
    
                                ax0.bar(ye_num, len(orientations_good_cluster), color=col, width=10, alpha=0.5,zorder=1) # Plot histogram cluster good
                                ax1.scatter([ye_num]*len(orientations_good_cluster), orientations_good_cluster, marker='.', color=col,s=snr_good_cluster*10, alpha=0.25, ec='k') # Plot orientations cluster good                          
                                ax2.scatter([ye_num]*len(clock_error_good), clock_error_good,color=col,marker='.',edgecolor='none',alpha=0.25) # Plot clock instability cluster good
    
                                g1 = ax3.scatter([ye_num]*len(gain_HHZ_good),np.log(gain_HHE_good/gain_HHN_good),color=col,marker='p',edgecolor='none',s=10,alpha=0.25,label='E/N')  # Plot E\N cluster good
                                g2 = ax3.scatter([ye_num]*len(gain_HHZ_good),np.log(gain_HHE_good/gain_HHZ_good),color=col,marker='>',edgecolor='none',s=10,alpha=0.25,label='E/Z')  # Plot E\Z cluster good
                                g3 = ax3.scatter([ye_num]*len(gain_HHZ_good),np.log(gain_HHN_good/gain_HHZ_good),color=col,marker='^',edgecolor='none',s=10,alpha=0.25,label='N/Z')  # Plot N\Z cluster good
            
                                if uni != -1:
    
                                    ax0.bar(ye_num, len(orientations_good_cluster), color=col, width=20,zorder=10) # Plot histogram cluster good
                                    ax1.scatter([ye_num]*len(orientations_good_cluster), orientations_good_cluster, marker='.', color=col,s=snr_good_cluster*10, alpha=0.5, ec='k') # Plot orientations cluster good
                                    ax1.annotate(f"{round(np.mean(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == uni)]['theta'].values), 1)}±{round(np.std(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == uni)]['theta'].values), 2)}°",(pd.to_datetime(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == uni)]['evtime'].values).mean(), round(np.mean(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == uni)]['theta'].values), 1)),fontsize=12, va='center', ha='center',path_effects=[path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
                                    ax2.scatter([ye_num]*len(clock_error_good), clock_error_good,color=col,marker='.',edgecolor='none',alpha=0.5) # Plot clock instability cluster good
    
                                    g1 = ax3.scatter([ye_num]*len(gain_HHZ_good),np.log(gain_HHE_good/gain_HHN_good),color=col,marker='p',edgecolor='none',s=10,alpha=0.5,label='E/N')  # Plot E\N cluster good
                                    g2 = ax3.scatter([ye_num]*len(gain_HHZ_good),np.log(gain_HHE_good/gain_HHZ_good),color=col,marker='>',edgecolor='none',s=10,alpha=0.5,label='E/Z')  # Plot E\N cluster good
                                    g3 = ax3.scatter([ye_num]*len(gain_HHZ_good),np.log(gain_HHN_good/gain_HHZ_good),color=col,marker='^',edgecolor='none',s=10,alpha=0.5,label='N/Z')  # Plot E\N cluster good
    
                    bad_values = len(df_sta[df_sta['quality'] == 'bad']['theta'].values)
                    label_dbscan_lst_end.append('bd: '+str(bad_values))
    
                    colors.append(mcolors.to_rgba('gray',alpha=0.25))  # add 'gray' as RGBA
                    
                    # Personalised handles
                    handles = [Line2D([0], [0], marker='o', color='w', label=grupo,markerfacecolor=cor,markeredgecolor='k',markersize=8) for grupo, cor in zip(label_dbscan_lst_end, colors)]
    
                    # ------------------------------ #
                    # Kernel density estimation plot #
                    # ------------------------------ #
                    
                    kde = gaussian_kde(orientations_all_good)
                    x_vals = np.linspace(min(orientations_all_good), max(orientations_all_good), 1000)
                    density = kde(x_vals)
                            
                    ax4.plot(density,x_vals, '-k')
                
                    # Histogram parameters
    
                    ax0.yaxis.set_major_locator(MultipleLocator(5))
                    ax0.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, rotation=30)
                    ax0.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax0.set_title(f'{net}.{sta}', fontsize=20)
                    ax0.set_ylim(0, 20)
                    ax0.set_ylabel("n")
                    ax0.grid(True)
                    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax0.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax0.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    
                    # Orientation parameters
                    
                    ax1.set_ylabel(r'Orientation($\theta$)')
                    ax1.set_ylim(-195, 195)
                    ax1.yaxis.set_major_locator(MultipleLocator(40))
                    ax1.yaxis.set_minor_locator(MultipleLocator(10))
                    ax1.yaxis.set_major_formatter(FuncFormatter(format_y_ticks))
                    ax1.grid(True)
                    ax1.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, top=True, rotation=30)
                    ax1.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax1.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
                    ax1.legend(handles=handles,loc='lower left',ncols=3)
                    
                    # Clock parameters
                    
                    ax2.set_ylim(-100, 100)
                    ax2.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, rotation=30)
                    ax2.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax2.grid(True)
                    ax2.yaxis.set_major_locator(MultipleLocator(100))
                    ax2.yaxis.set_minor_locator(MultipleLocator(20))
                    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax2.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
                    ax2.set_ylabel("Time")
    
                    # Gain parameters
                    
                    ax3.figure.legend(handles=[g1, g2, g3],loc='center',bbox_to_anchor=(0.85, 0.14),frameon=False,ncol=1,fontsize=10,borderaxespad=0.)
                    ax3.tick_params(axis="x", which='both', labelbottom=True, labeltop=False, rotation=30)
                    ax3.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax3.grid(True)
                    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax3.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax3.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
                    ax3.set_ylabel("Gain")
    
                    # KDE parameters
                    
                    ax4.set_xlabel('KDE')
                    ax4.xaxis.set_label_position('top')
                    ax4.yaxis.set_major_formatter(FuncFormatter(format_y_ticks))
                    ax4.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, rotation=30)
                    ax4.tick_params(axis="y", which='both', labelright=True, labelleft=False,labelbottom=False, labeltop=False, left=True, right=True,top=True, bottom=False)
                    ax4.grid(True)
                    
                    
                    # Saving the figure
                    
                    output_figure_SSPARQ = SSPARQ_OUTPUT + 'FIGURES/FINAL_RESULT/'+net+'/'
                    os.makedirs(output_figure_SSPARQ, exist_ok=True)
                    fig.savefig(output_figure_SSPARQ + f'ORIENTATION_TOTAL_{net}_{sta}.png',facecolor='w',dpi=300)
                    plt.close()
                
                else:
                
                    # ---------
                    # NO DBSCAN
          
                    # -----------------------
                    # Outlier analysis result
                            
                    df_sta['class'] = -10  # standart value
                                 
                    # Filter mask
                    mask_good,mask_outliers = calculate_quartis_mask(orientations_all_good)
        
                    df_sta.loc[df_sta['quality'] == 'good', 'class'] = [1 if i == True else -1 for i in mask_good]
                
                    # Creating mounth range datetime (timestamp) Periods
                    years = pd.period_range(start=df_sta['year_month'].min(), end=df_sta['year_month'].max(), freq='M')
                    
                    fig = plt.figure(figsize=(10, 10))
                    gs = gridspec.GridSpec(4, 2, width_ratios=[10,1], height_ratios=[1,10,1,1],hspace=0.05, wspace=0.05)
                
                
                    # Definindo os eixos
                    ax0 = fig.add_subplot(gs[0, 0])  #  Number of events axis
                    ax1 = fig.add_subplot(gs[1, 0],sharex=ax0)  # Orientation axis
                    ax2 = fig.add_subplot(gs[2, 0],sharex=ax0)  # Histogram axis
                    ax3 = fig.add_subplot(gs[3, 0],sharex=ax0)  # Time axis
                    ax4 = fig.add_subplot(gs[1, 1],sharey=ax1)  # KDE axis
                    
                    if silhouette_score_elbow_point:
                        ax0.annotate('ss:'+str(round(silhouette_score_elbow_point,2)), (pd.to_datetime(df_sta['evtime'].values).max(), +15),fontsize=10, va='center', ha='center',bbox=dict(boxstyle="round", fc="white", ec='k', alpha=0.5))
                    
                    label_handles_bad = []
                    label_handles_out = []
                    label_handles_dat = []
                    
                    for ye in years:
                        # Filtering according to Period (year_month)
                        df_sta_year = df_sta[df_sta['year_month'] == ye]
                
                        ye_num = mdates.date2num(ye)  # converte timestamp para número
                        
                        # ---------------- #
                        # Orientation plot #
                        # ---------------- #
                        
                        if df_sta_year[df_sta_year['quality'] == 'good']['theta'].empty:
                            orientations_bad = df_sta_year[df_sta_year['quality'] == 'bad']['theta'].values
                            snr_bad = df_sta_year[df_sta_year['quality'] == 'bad']['SNR'].abs().values
    
                            ax0.bar(ye_num, len(orientations_bad), color='gray', width=20, alpha=0.25,zorder=-10) # Plot histogram bad
                            ax1.scatter([ye_num]*len(orientations_bad), orientations_bad, marker='.', c='gray',s=snr_bad*10, alpha=0.05, ec='k') # Plot orientations bad
        
                            label_handles_bad.append(len(orientations_bad))
                            
                            # time plot #
    
                            clock_error_bad = df_sta_year[df_sta_year['quality'] == 'bad']['clock_error'].values # Clock instability
    
                            ax2.scatter([ye_num]*len(clock_error_bad), clock_error_bad,c='gray',marker='.',edgecolor='none',alpha=0.1) # Plot clock instability bad
    
                            # gain plot #
    
                            gain_HHE_bad = df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHE'].values.mean() # Gain HHE bad
                            gain_HHN_bad= df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHN'].values.mean() # Gain HHN bad
                            gain_HHZ_bad = df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHZ'].values.mean() # Gain HHZ bad
        
                            g1 = ax3.scatter(ye_num,np.log(gain_HHE_bad/gain_HHN_bad),c='gray',marker='p',edgecolor='none',s=10,alpha=0.25,label='E/N') # Plot E/N bad
                            g2 = ax3.scatter(ye_num,np.log(gain_HHE_bad/gain_HHZ_bad),c='gray',marker='>',edgecolor='none',s=10,alpha=0.25,label='E/Z') # Plot E/Z bad
                            g3 = ax3.scatter(ye_num,np.log(gain_HHN_bad/gain_HHZ_bad),c='gray',marker='^',edgecolor='none',s=10,alpha=0.25,label='N/Z') # Plot N/Z bad
       
                        else:
                            orientations_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == 1)]['theta'].values
                            orientations_out = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == -1)]['theta'].values
                            orientations_bad = df_sta_year[df_sta_year['quality'] == 'bad']['theta'].values
                                        
                            snr_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == 1)]['SNR'].abs().values
                            snr_out = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == -1)]['SNR'].abs().values
                            snr_bad = df_sta_year[df_sta_year['quality'] == 'bad']['SNR'].abs().values
        
                            # Bad orientations
                            
                            ax0.bar(ye_num, len(orientations_bad), color='gray', width=20, alpha=0.25,zorder=-9) # Plot histogram bad
                            ax1.scatter([ye_num]*len(orientations_bad), orientations_bad, marker='.', c='gray',s=snr_bad*10, alpha=0.05, ec='k') # Plot orientations bad
                            
                            label_handles_bad.append(len(orientations_bad))
    
                            # Outliers orientations
    
                            ax0.bar(ye_num, len(orientations_out), color='mediumpurple', width=20, alpha=0.5,zorder=-5) # Plot histogram out
                            ax1.scatter([ye_num]*len(orientations_out), orientations_out, marker='.', c='mediumpurple',s=snr_out*10, alpha=0.25, ec='k') # Plot orientations out
                            
                            # Good orientations
                            
                            ax0.bar(ye_num, len(orientations_good), color='#9e0039', width=20, alpha=0.75,zorder=1) # Plot histogram good                                       
                            ax1.scatter([ye_num]*len(orientations_good), orientations_good, marker='.', c='#9e0039',s=snr_good*10, alpha=0.75, ec='k') # Plot orientations good
                            
                            label_handles_out.append(len(orientations_out))
                            label_handles_dat.append(len(orientations_good))
    
                            # time plot #
    
                            clock_error_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == 1)]['clock_error'].values # Clock instability good
                            clock_error_out = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == -1)]['clock_error'].values # Clock instability out
                            clock_error_bad = df_sta_year[(df_sta_year['quality'] == 'bad')]['clock_error'].values # Clock instability bad
                            
                            ax2.scatter([ye_num]*len(clock_error_bad), clock_error_bad,c='gray',marker='.',edgecolor='none',alpha=0.1) # Plot clock instability bad
                            ax2.scatter([ye_num]*len(clock_error_out), clock_error_out,c='mediumpurple',marker='.',edgecolor='none',alpha=0.1) # Plot clock instability out
                            ax2.scatter([ye_num]*len(clock_error_good), clock_error_good,c='#9e0039',marker='.',edgecolor='none',alpha=0.5) # Plot clock instability out
    
                            gain_HHE_bad = df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHE'].values.mean() # Gain HHE bad
                            gain_HHN_bad= df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHN'].values.mean() # Gain HHN bad
                            gain_HHZ_bad = df_sta_year[df_sta_year['quality'] == 'bad']['gain_HHZ'].values.mean() # Gain HHZ bad
        
                            g1 = ax3.scatter(ye_num,np.log(gain_HHE_bad/gain_HHN_bad),c='gray',marker='p',edgecolor='none',s=10,alpha=0.25,label='E/N') # Plot E/N bad
                            g2 = ax3.scatter(ye_num,np.log(gain_HHE_bad/gain_HHZ_bad),c='gray',marker='>',edgecolor='none',s=10,alpha=0.25,label='E/Z') # Plot E/Z bad
                            g3 = ax3.scatter(ye_num,np.log(gain_HHN_bad/gain_HHZ_bad),c='gray',marker='^',edgecolor='none',s=10,alpha=0.25,label='N/Z') # Plot N/Z bad
                            
                            gain_HHE_out = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == -1)]['gain_HHE'].values.mean() # Gain HHE out
                            gain_HHN_out = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == -1)]['gain_HHN'].values.mean() # Gain HHN out
                            gain_HHZ_out = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == -1)]['gain_HHZ'].values.mean() # Gain HHZ out
    
                            g1 = ax3.scatter(ye_num,np.log(gain_HHE_out/gain_HHN_out),c='mediumpurple',marker='p',edgecolor='none',s=10,alpha=0.25,label='E/N') # Plot E/N out
                            g2 = ax3.scatter(ye_num,np.log(gain_HHE_out/gain_HHZ_out),c='mediumpurple',marker='>',edgecolor='none',s=10,alpha=0.25,label='E/Z') # Plot E/Z out
                            g3 = ax3.scatter(ye_num,np.log(gain_HHN_out/gain_HHZ_out),c='mediumpurple',marker='^',edgecolor='none',s=10,alpha=0.25,label='N/Z') # Plot N/Z out
    
                            gain_HHE_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == 1)]['gain_HHE'].values.mean() # Gain HHE good
                            gain_HHN_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == 1)]['gain_HHN'].values.mean() # Gain HHN good
                            gain_HHZ_good = df_sta_year[(df_sta_year['quality'] == 'good') & (df_sta_year['class'] == 1)]['gain_HHZ'].values.mean() # Gain HHZ good
    
                            g1 = ax3.scatter(ye_num,np.log(gain_HHE_good/gain_HHN_good),c='#9e0039',marker='p',edgecolor='none',s=10,alpha=0.5,label='E/N') # Plot E/N good
                            g2 = ax3.scatter(ye_num,np.log(gain_HHE_good/gain_HHZ_good),c='#9e0039',marker='>',edgecolor='none',s=10,alpha=0.5,label='E/Z') # Plot E/Z good
                            g3 = ax3.scatter(ye_num,np.log(gain_HHN_good/gain_HHZ_good),c='#9e0039',marker='^',edgecolor='none',s=10,alpha=0.5,label='N/Z') # Plot N/Z good
                            
                    ax1.annotate(f"{round(np.mean(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == 1)]['theta'].values), 1)}±{round(np.std(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == 1)]['theta'].values), 2)}°",(pd.to_datetime(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == 1)]['evtime'].values).mean(), round(np.mean(df_sta[(df_sta['quality'] == 'good') & (df_sta['class'] == 1)]['theta'].values), 1)),fontsize=12, va='center', ha='center',path_effects=[path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
                           
                    label_handles = ['bd: '+str(sum(label_handles_bad)),'ol: '+str(sum(label_handles_out)),'g1: '+str(sum(label_handles_dat))]
                    alphas = [0.05,0.25,0.75] 
                    colors = ['gray','mediumpurple','darkred']
                    colors_with_alpha = [mcolors.to_rgba(cor, alpha=a) for cor, a in zip(colors, alphas)]
                    
                    # Personalised handles
                    handles = [Line2D([0], [0], marker='o',color='w', label=grupo,markerfacecolor=cor,markeredgecolor='k',markersize=8) for grupo, cor in zip(label_handles, colors_with_alpha)]
    
                    # ------------------------------ #
                    # Kernel density estimation plot #
                    # ------------------------------ #
                    
                    kde = gaussian_kde(orientations_all_good)
                    x_vals = np.linspace(min(orientations_all_good), max(orientations_all_good), 1000)
                    density = kde(x_vals)
                            
                    ax4.plot(density,x_vals, '-k')
                
                    # Histogram parameters
    
                    ax0.yaxis.set_major_locator(MultipleLocator(5))
                    ax0.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, rotation=30)
                    ax0.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax0.set_title(f'{net}.{sta}', fontsize=20)
                    ax0.set_ylim(0, 20)
                    ax0.set_ylabel("n")
                    ax0.grid(True)
                    ax0.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax0.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax0.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
    
                    # Orientation parameters
                    
                    ax1.set_ylabel(r'Orientation($\theta$)')
                    ax1.set_ylim(-195, 195)
                    ax1.yaxis.set_major_locator(MultipleLocator(40))
                    ax1.yaxis.set_minor_locator(MultipleLocator(10))
                    ax1.yaxis.set_major_formatter(FuncFormatter(format_y_ticks))
                    ax1.grid(True)
                    ax1.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, top=True, rotation=30)
                    ax1.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax1.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
                    ax1.legend(handles=handles,loc='lower left',ncols=3)
                    
                    # Clock parameters
                    
                    ax2.set_ylim(-100, 100)
                    ax2.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, rotation=30)
                    ax2.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax2.grid(True)
                    ax2.yaxis.set_major_locator(MultipleLocator(100))
                    ax2.yaxis.set_minor_locator(MultipleLocator(20))
                    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax2.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
                    ax2.set_ylabel("Time")
    
                    # Gain parameters
                    
                    ax3.figure.legend(handles=[g1, g2, g3],loc='center',bbox_to_anchor=(0.85, 0.14),frameon=False,ncol=1,fontsize=10,borderaxespad=0.)
                    ax3.tick_params(axis="x", which='both', labelbottom=True, labeltop=False, rotation=30)
                    ax3.tick_params(axis="y", which='both', labelright=False, labelleft=True, left=True, right=True)
                    ax3.grid(True)
                    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
                    ax3.xaxis.set_major_locator(mdates.MonthLocator(interval=18))
                    ax3.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
                    ax3.set_ylabel("Gain")
    
                    # KDE parameters
                    
                    ax4.set_xlabel('KDE')
                    ax4.xaxis.set_label_position('top')
                    ax4.yaxis.set_major_formatter(FuncFormatter(format_y_ticks))
                    ax4.tick_params(axis="x", which='both', labelbottom=False, labeltop=False, rotation=30)
                    ax4.tick_params(axis="y", which='both', labelright=True, labelleft=False,labelbottom=False, labeltop=False, left=True, right=True,top=True, bottom=False)
                    ax4.grid(True)
                    
                    # Salvando a figura
                    output_figure_SSPARQ = SSPARQ_OUTPUT + 'FIGURES/FINAL_RESULT/'+net+'/'
                    os.makedirs(output_figure_SSPARQ, exist_ok=True)
                    fig.savefig(output_figure_SSPARQ + f'METRICS_TOTAL_{net}_{sta}.png',facecolor='w',dpi=300)
                    plt.close()
        else:
            print('station:',sta,'size:',len(df_sta))