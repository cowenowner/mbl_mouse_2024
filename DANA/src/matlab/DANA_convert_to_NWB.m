% NWB (Neurodata Without Borders) is the format of data for uploading that
% UPDATE: UTTER FAILURE - dies on first line of code (see below). Damn
% academics and thier pre-alpha software. Going to try the python
% interface.
%
% is required for the Brain Initiative Projects. NWB is foundational for
% the Distributed Archives for Neurophysiology Data Integration (DANDI)
% data repository to enable collaborative data sharing and analysis.
% (https://dandiarchive.org)
% Key instructions on data conversion.
% https://neurodatawithoutborders.github.io/matnwb/tutorials/html/ecephys.html
% https://nwb-overview.readthedocs.io/en/latest/conversion_tutorial/
%
% https://www.braininitiative.org/toolmakers/resources/neurodata-without-borders-nwb/
% https://www.biorxiv.org/content/10.1101/2021.03.13.435173v2.full
%
% This script will demonstrate how to enter that data. There are Matlab and
% Python interfaces. If you understand the Matlab version, Python is nearly
% identical. PyNWB is the Python API.
%
% Interfacing: NWBWidgets is a library for interactive web-based visualizations of NWB data, that provides tools for navigating and combining data across modalities.
% neurodata types are divided into modules such as ecephys (extracellular
% electrophysiology), icephys (intracellular electrophysiology), behavior
%
% an ElectricalSeries is a neurodata type that defines the data and metadata for an intracranially recorded voltage time series in an extracellular electrophysiology experiment.
% TimeSeries neurodata type, which is a generic structure designed for any measurement that is sampled over time, and defines fields, such as, data, units of measurement, and sample times (specified either with timestamps or sampling rate and start time)
%
% DANDI can handle large volumes (TBs) of data and host them for free, and
% provides a command-line interface (CLI) that is built to handle this
% volume of data. DANDI also automatically parses NWB files to extract
% metadata that makes these datasets easier for others to find.
%
% Cowen 2022
%

% STEP 1: Make sure matnwb and its subdirs is in your path.
% 
% This line of code will add it to your path.
addpath(genpath(fullfile(Git_dir,'matnwb')))

% Step 2: decide what data directory we are going to analyze.
data_dir = 'Z:\Data\ACUTE_DANA\220408_MFB_NAc';

% Step 2: set up some meta data and nwb file. 
% DAMMIT - IT CRASHES HERE ALREADY!!!! 
% Error using NwbFile
% The specified superclass 'types.core.NWBFile' contains a parse error, cannot be found on MATLAB's search path, or is shadowed by
% another file with the same name.
nwb = NwbFile(  'session_description', 'FSCV and single unit anesthetized',...
    'identifier', '220408_MFB_NAc', ...
    'session_start_time', datetime(2022, 4, 8, 10, 00, 00), ...
    'general_experimenter', 'Abhi Vishwanath', ... % optional
    'general_session_id', 'FSCV_and_Ephys', ... % optional
    'general_institution', 'University of Arizona'); 

nwb

