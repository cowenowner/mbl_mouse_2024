%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NPXL_Post_Process_Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The high-level code for running all of the post processing.
% Feel free to copy this and make a version for your own data
% post-processing. This code runs CatGT, extracts events, and
% filters/proceses .ap and .lf files for spike sorting.
%
% Expects:
% An event_codes.csv file in the top directory. It contains the event
% channel, variable name, and notes.
% C:\CatGT-win
% C:\TPrime-win
% Assumes CowenLib
% Assumes in your GitHub folder: https://github.com/djoshea/neuropixel-utils
%

% You must run kilosort manually throught the GUI.
% You must run Tprime after spike sorting to ensure aligned spikes and
% event times.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2022
%% SS EDITS 2023/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; fclose all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change variables here for your analysis.
  %  C:\SGL_DATA\myRun8_crus2_g0\myRun8_crus2_g0_imec0
%%DataDirs
%D:\Data\M521_2024_07_04_R_fiber_lightstim_wheel_pixel1-0_g0
PRM_RAT_ROOT='D:\Data\'; %'C:\SGL_DATA\'  - the root directory above all the individual data directories.
PRM_RAT_SUBDIRS={'M521_2024_07_04_R_fiber_lightstim_wheel_pixel1-0_g0' };%'mPFC_L5_bank0_g0','mPFC_L5_bank0_g0', 'myRun8_crus2_g0'
PRM_BAD_CHANNEL0_LIST = []; % This is ZERO based as you would see in SpikeGLX so be sure the first channel is zero.
PRM_TEMP_FOLDER_BASE = 'C:\Temp\kilosort'; % This needs to be a SSD.

%PROBE Type Neuropixels 2.0 or 1.0

%Functions to run
PRM_CREATE_TCAT_FILE =true; % make false if you already created this file on a previous run to save some time.
PRM_CREATE_LF_FILE=true; %make false if you don't want to run the LF files. Saves time
PRM_COPY_FILES=false; %If you want a copy of the tcat files copied from the Base dir to denoising/kilosort dir
PRM_CREATE_CHANNELMAP=true; %Will generate a new channel map if you are using a custom config
PRM_SPLIT_SHANKS=false;% NP2 This will split the split the channel map into separate shanks

PRM_RUN_DENOISE=false;
PRM_RUN_KILOSORT=false; %DONNOT RUN

% NP2
% FOLL TWO FUNCTIONS NP2 and STIM ONLY SS EDITS 
PRM_NP2=false;%Default set to false if this is NP1.0

%Note: if this is set to true and CREATE LF is set to True then it will not
%run AP CATGT Separately. The LF_NP2 CatGT  runs both since it has to use
%the AP file to generate the LF file for NP2
%OVERSTRIKE TO REMOVE ARTIFACTS
PRM_RUN_OVERSTRIKE=false;
%IMEC To REMOVE NOISE %Removes stim noise for NP2 files
PRM_RUN_REM_NOISE=false;

%Need to add the following to the run functions
PRM_STIM_FILE_EXT='*.xa_5_0.txt'; %This uses stim files from the xa stream of ch5 (after catGT runs) 
zero_time_win_in_s=0.003;%  Removes a 2 ms window post stim as default
PRM_CAGE_TIMES=[]; %[2040 2160] in seconds to remove any values within this time window

%Note if you run Channel Map make sure the meta file is in the same
%location where data was originally collected!!
%Otherwise uncomment the following two lines and add the parameters to the run command below

%PRM_chanMapPath='';
%PRM_chanMapFile'';

%%%%%%%%%%CatGT Parameters
%AP FILES

%LF FILES

%ANALOG FILES
PRM_EVENT_CHANNELS=[0 1 2 3 4 5]; %Which analog channels to extract
PRM_XA_Duration = [0]; %Time of pulses in ms - creates separate files. Good to keep 0 as a default
PRM_INAROW=2;
PRM_XA_THR1=0.8;%Threshold for detecting xa (in V)
PRM_XA_THR2=0.5; %Secondary threshold for detetcting xa (in V) -keep lower than THR1 for square pulses
%Default Parameters for extracting 
% PRM_EVENT_CHANNELS = [1:7];% by default, get everything to avoid mistakes.
% PRM_XA_THR1=3.0;%Threshold for detecting xa (in V)
% PRM_XA_THR2=2.0;%Secondary threshold for detetcting xa (in V) -keep lower than THR1 for square pulses
% PRM_XA_Duration=[0.0];%Default timestamp detects all edges (in ms)
% PRM_XIA_THR1=-3.0;%inverse pulse Threshold for detecting xa (in V) inverse pulse
% PRM_XIA_THR2=-2.0;% inverse pulseSecondary threshold for detetcting xa (in V) -keep lower than THR1 for square pulses
% PRM_INAROW=3; %Min No. of continuous samples to detect for analog pulse

%Note if you uncomment any of the above parameters to change from default,
%you will need to declare the variable in the command below (varargin) and
%in the NPXL_Post_Process_SS file. 

PRM_CatGT_dir = 'C:\CatGTWinApp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(PRM_RAT_SUBDIRS) %Specify which subfolders here
    PRM_ROOT_DATA_DIR=fullfile(PRM_RAT_ROOT,PRM_RAT_SUBDIRS{ii});
    PRM_TEMP_FOLDER_LOCATION=fullfile(PRM_TEMP_FOLDER_BASE,PRM_RAT_SUBDIRS{ii});
    
    %Run the PostProcess 
    NPXL_Post_Process_SS('PRM_ROOT_DATA_DIR',PRM_ROOT_DATA_DIR,...
        'PRM_BAD_CHANNEL0_LIST',PRM_BAD_CHANNEL0_LIST,...
        'PRM_TEMP_FOLDER_LOCATION',PRM_TEMP_FOLDER_LOCATION,...
        'PRM_CREATE_TCAT_FILE',PRM_CREATE_TCAT_FILE,...
        'PRM_CREATE_LF_FILE',PRM_CREATE_LF_FILE,...
        'PRM_COPY_FILES',PRM_COPY_FILES,...
        'PRM_CREATE_CHANNELMAP',PRM_CREATE_CHANNELMAP,...
        'PRM_RUN_DENOISE',PRM_RUN_DENOISE,...
        'PRM_RUN_KILOSORT',PRM_RUN_KILOSORT,...
        'PRM_EVENT_CHANNELS',PRM_EVENT_CHANNELS,...
        'PRM_XA_Duration',PRM_XA_Duration,...
        'PRM_INAROW',PRM_INAROW,...
        'PRM_XA_THR1',PRM_XA_THR1,...
        'PRM_XA_THR2',PRM_XA_THR2,...
        'PRM_NP2',PRM_NP2,...
        'PRM_RUN_OVERSTRIKE',PRM_RUN_OVERSTRIKE, ...
        'PRM_CAGE_TIMES',PRM_CAGE_TIMES,...
        'PRM_RUN_REM_NOISE',PRM_RUN_REM_NOISE);
end

        