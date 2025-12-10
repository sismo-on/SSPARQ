"""
--------------------------------------------------------------------------------
         Module that parses global parameters from a configuration file
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 08/2025

Description:
Module that parses global parameters from a configuration file at first import,
to make them available to the other parts of the program.

More information in:
https://wiki.python.org/moin/ConfigParserExamples

Input:
Configuration file, wherein global paths and parameters are defined.

Outputs:
The module provides a parser for simple configuration files consisting of groups
of named values.

"""

import configparser
import os
import glob


def select_and_parse_config_file(basedir='.', ext='cnf', verbose=True):
    """
    Reads a configuration file and returns an instance of ConfigParser:
    First, looks for files in *basedir* with extension *ext*.
    Asks user to select a file if several files are found,
    and parses it using ConfigParser module.
    @rtype: L{ConfigParser.ConfigParser}
    """
    config_files = glob.glob(os.path.join(basedir, u'*.{}'.format(ext)))


    if not config_files:
        raise Exception("No configuration file found!")

    if len(config_files) == 1:
        # only one configuration file
        config_file = config_files[0]
    else:
        print("Select a configuration file:")
        for i, f in enumerate(config_files, start=1):
            print("{} - {}".format(i, f))
        res = int(input(''))
        config_file = config_files[res - 1]

    if verbose:
        print("Reading configuration file: {}".format(config_file))

    conf = configparser.ConfigParser(allow_no_value=True)
    conf.read(config_file)

    return conf

# ==========================
# parsing configuration file
# ==========================

config = select_and_parse_config_file(basedir='.', ext='cnf', verbose=True)

# -----
# paths
# -----

## ------------------------
## Directory with waveforms (SeisComP Data Structure)

WAVEFORM_DIR = config.get('paths', 'WAVEFORM_DIR')

## ---------------------------
## Directory with the catalog (.CSV file of the National Earthquake Information Center (NEIC))

CATALOG_FILE = config.get('paths', 'CATALOG_FILE')

## ----------------------------
## Directory of the StationXML:

XML_DIR = config.get('paths', 'XML_DIR')

## -----------------------
## Directory of the output (Figures and Feathers file)

SSPARQ_OUTPUT = config.get('paths', 'SSPARQ_OUTPUT')

# ------
# event
# ------

## -------------------------------------------------------------------
## Taup_time model to calculate travel times
TAUPY_MODEL = config.get('event', 'TAUPY_MODEL')

## -------------------------------------------------------------------
## Apply band-pass filtering to the seismograms using the range above:

PERIOD_BANDS_MAX = config.getfloat('event', 'PERIOD_BANDS_MAX')

PERIOD_BANDS_MIN = config.getfloat('event', 'PERIOD_BANDS_MIN')

## ===================================================================================
## Default parameters to define the signal and noise windows used to estimate the SNR:

## ------------------------------------------------------------------------------
## Duration of the signal window before and after the P-wave arrival (in seconds)

TIME_WINDOW = config.getfloat('event', 'TIME_WINDOW')

## ---------------------------------------------------------------
## P-wave time window for events (in seconds after P-wave arrival)

TIME_FINAL_P = config.getfloat('event', 'TIME_FINAL_P')

## ---------------------------------------------------------------
## If True, the instrumental response is removed from the trace.
## Default is False.

REMOVE_RESPONSE = config.getboolean('event', 'REMOVE_RESPONSE')

## ---------------------------------------------------------------
## Physical unit of the output signal after response removal.
## Accepted values are DISP, VEL, or ACC (Default is VEL).

OUTPUT_UNIT = config.get('event', 'OUTPUT_UNIT')

## ---------------
## MULTIPROCESSING

num_processes = config.getint('event', 'num_processes')

# --------
# algoritm
# --------
## -------------------------------------------------------------
## Minimum required similarity of vertical and radial components (default is 0.5).

CVR_MIN = config.getfloat('algoritm', 'CVR_MIN')

## --------------------------------------
## Minimum required signal-to-noise ratio (default is 10).

SNR_MIN = config.getfloat('algoritm', 'SNR_MIN')

## --------------------------------------------------
## Minimum required transverse-to-radial energy ratio (default is 0.2).

TRR_MIN = config.getfloat('algoritm', 'TRR_MIN')

## -----------------------------------------------
## Maximum allowed radial-to-vertical energy ratio (default is 2).

RVR_MAX = config.getfloat('algoritm', 'RVR_MAX')

## -------------------------------------------------
## Optimization method (default: grid_search)
## Valid options:
## - grid_search   : exhaustive fine search with a 0.1-degree step, without refinement
## - hybrid_search : coarse search with 10-degree steps, followed by Newton's method refinement (faster)

ORIENTATION_METHOD = config.get('algoritm', 'ORIENTATION_METHOD')

## -----------------------------------------------------------------
## Minimum number of measurements required to compute the final plot (default = 1).

MIN_ORIENTATION_STATION = config.getfloat('algoritm', 'MIN_ORIENTATION_STATION')

## -----------------------------------------------------------------------
## DBSCAN CLUSTERING PARAMETERS
## -----------------------------------------------------------------------
## This section defines the main configuration parameters used in the
## DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
## algorithm. DBSCAN identifies clusters based on the density of points
## in the feature space, requiring two key parameters: `epsilon` and 
## `min_samples`. The configuration below sets the search limits and 
## scaling conditions for cluster estimation.
## 
## 
## Percentage of total samples per group (default = 20)
## Defines the fraction of total data points required to
## form a dense region. Used to compute the `min_samples`
## parameter internally, as:
##     min_samples = total_samples * (PER_SAMPLES / 100)
## Larger values increase the density requirement, making
## clusters harder to form and increasing the number of
## points labeled as noise.

PER_SAMPLES = config.getfloat('algoritm', 'PER_SAMPLES') 

## -----------------------------------------------------------------------
## Epsilon parameter range (distance threshold)
## The `epsilon` value defines the neighborhood radius —
## the maximum distance between two samples for them to be
## considered part of the same cluster.
## Smaller epsilon values generate more fragmented clusters,
## while larger ones merge nearby clusters into broader groups.
## These limits are used to scan the parameter space and 
## identify the most suitable value for the dataset.
##
## EPSILON_LOW  - lower bound of the epsilon range  (default = 0.2)
## EPSILON_UP   - upper bound of the epsilon range  (default = 0.3)

EPSILON_LOW = config.getfloat('algoritm', 'EPSILON_LOW') 
EPSILON_UP  = config.getfloat('algoritm', 'EPSILON_UP') 

## -----------------------------------------------------------------------
## The silhouette score evaluates the quality of a DBSCAN clustering. 
## The Silhouette index ranges from -1 to 1, where negative values 
## indicate higher within-cluster dissimilarity than between-cluster 
## dissimilarity, while positive values suggest a higher likelihood 
## of correct clustering. The MIN_SILHOUETTE_SCORE is the optimal 
## minimal value for distinguishing well-defined clusters from noise.
## (default = 0.2)

MIN_SILHOUETTE_SCORE = config.getfloat('algoritm', 'MIN_SILHOUETTE_SCORE') 