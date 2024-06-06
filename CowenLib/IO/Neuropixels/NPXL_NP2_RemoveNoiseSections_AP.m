function NPXL_NP2_RemoveNoiseSections_AP(varargin)
% NPXL_NP2_RemoveNoiseSections_AP(varargin)
% Removes sections of the AP file which are noisy. Will have to rescale timestamps later after kilosort
% KILOSORT 2.5 doesn't like zerod out data. 
% Need to save excised TimeStamps in a text file 
% Used to remove the stim artifact
% BE SURE TO RUN THIS ON A tcat.ap.bin file - not the original data file as
% CatGT does some subtle but important inter-channel aligments.
% 
% NEED: 
%        Directory where the AP.Bin File is located
%       List of stim times
%       Time interval post stim to zero out 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SS 2023
%Load the neuropixel 2 utils. This has some difference from the original
%that make it easier to read and reads NP2 files correctly
PRM_NEUROPIXEL_UTILS_DIR_NP2 = fullfile(Git_dir,'neuropixel-utils_SS'); % if this works for you, go for it
% Otherwise, do this if you do not have the folder.
% PRM_NEUROPIXEL_UTILS_DIR = 'C:\Users\cowen\Documents\GitHub\neuropixel-utils_SS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial Parameters
addpath(PRM_NEUROPIXEL_UTILS_DIR_NP2)
PRM_ROOT_DATA_DIR =''; %Default assumes current folder
PRM_STIM_FILE_EXT='*.xa_5_0.txt';
zero_time_win_in_s=0.002;
PRM_CAGE_TIMES=[]; %in seconds to remove any values within this time window
PRM_BAD_CHANNEL0_LIST=[];
CHANNEL_MAP_FILE=[];
PRM_TEMP_FOLDER_LOCATION=pwd;
Extract_varargin; % overrides the defaults specified above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Get the  stim file names and times 
stim_file=dir(fullfile(PRM_ROOT_DATA_DIR,PRM_STIM_FILE_EXT));
stim_times = readmatrix(fullfile(stim_file.folder,stim_file.name));
zero_times=[stim_times-zero_time_win_in_s stim_times+zero_time_win_in_s];

%Remove Cage Stim Times if needed
if ~isempty(PRM_CAGE_TIMES)
   zero_times(end+1,:)=PRM_CAGE_TIMES;
   zero_times = sortrows(zero_times, 1); %Reorder 
end
%Get data file
ap_dir = dir(fullfile(PRM_ROOT_DATA_DIR,'*imec0*'));
PRM_AP_DIR=fullfile(ap_dir.folder,ap_dir.name);
d = dir(fullfile(PRM_AP_DIR,'*tcat.*.ap.bin'));
ap_path=fullfile(d.folder,d.name);
meta_path=strrep(ap_path,'.ap.bin','.ap.meta');

%Output file
[~,n,ext] = fileparts(ap_path);
dest_ap_fname = fullfile(PRM_TEMP_FOLDER_LOCATION,[ 'denoise_' n ext] );
dest_meta_fname = strrep(dest_ap_fname,'.ap.bin','.ap.meta');
if isempty(d)
    error('could not find tcat file.')
end

%Load Imec Data set
imec = Neuropixel.ImecDataset(ap_path, 'channelMap',CHANNEL_MAP_FILE);
% Read meta file
meta=SGLX_readMeta.ReadMeta(d.name ,d.folder);
%imec.inspectAP_timeWindow([zero_times(1,1) zero_times(1,2)+0.2])
%Times
start_time_all=0; %Starts at 0secs 
end_time_all=str2double(meta.fileTimeSecs); %Get File times secs from meta file
%Get time windows
start_times = zeros(length(zero_times) + 1, 1);
end_times = zeros(length(zero_times) + 1, 1);
start_times(1) = start_time_all;
end_times(1) = zero_times(1, 1);
for ii=1:length(zero_times)-1
    start_times(ii + 1) = zero_times(ii, 2);
    end_times(ii + 1) = zero_times(ii + 1, 1);
end
start_times(end)=zero_times(end, 2);
end_times(end) = end_time_all;
% Convert Times to idx
start_idx=imec.closestSampleAPForTime(start_times);
end_idx=imec.closestSampleAPForTime(end_times);
%Get the time shifted data
timeShifts = Neuropixel.TimeShiftSpec.buildToExciseGaps(start_idx, end_idx);
imecOut = imec.saveTransformedDataset(dest_ap_fname, 'timeShiftsAP', timeShifts ,'writeAP', true);
% Save the TimeShift Info
time_shift_file=strrep(dest_meta_fname,'.imec0.ap.meta','_timeShift.mat');
timeShifts_save=struct(timeShifts);
save(time_shift_file,'timeShifts_save')
end