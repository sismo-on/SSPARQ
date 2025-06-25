"""
--------------------------------------------------------------------------------
         Module that parses global parameters from a configuration file
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (diogoloc@on.br)


Last Date: 06/2025

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

## ---------------
## MULTIPROCESSING

num_processes = config.getint('event', 'num_processes')

