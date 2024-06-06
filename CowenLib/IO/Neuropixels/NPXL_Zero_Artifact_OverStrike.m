function NPXL_Zero_Artifact_OverStrike(varargin)
% function NPXL_Zero_Artifact_OverStrike(varargin)
% Zeros out signal in ap file between specific time points. 
%
% Used to remove the stim artifact
% BE SURE TO RUN THIS ON A tcat.ap.bin file - not the original data file as
% CatGT does some subtle but important inter-channel aligments.
% 
%
% Download Overstrike from https://billkarsh.github.io/SpikeGLX/#post-processing-tools
% NEED: Directory where overstrike is located
%        Directory where the AP.Bin File is located
%       List of stim times
%       Time interval post stim to zero out 
% CAVEATS: THis can zero out only one time block per call of overstrike .
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SS 2023
PRM_OVERSTRIKE_DIR='C:\OverStrike-win';

PRM_ROOT_DATA_DIR =''; %Default assumes current folder
PRM_STIM_FILE_EXT='*.xa_5_0.txt';
zero_time_win_in_s=0.002;
PRM_CAGE_TIMES=[]; %in seconds to remove any values within this time window
PRM_BAD_CHANNEL0_LIST=[];
CHANNEL_MAP_FILE=[];
Extract_varargin; % overrides the defaults specified above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the  stim file names and times 
stim_file=dir(fullfile(PRM_ROOT_DATA_DIR,PRM_STIM_FILE_EXT));
stim_times = readmatrix(fullfile(stim_file.folder,stim_file.name));
zero_times=[stim_times stim_times+zero_time_win_in_s];

%Get data file
ap_dir = dir(fullfile(PRM_ROOT_DATA_DIR,'*imec0*'));
PRM_AP_DIR=fullfile(ap_dir.folder,ap_dir.name);
d = dir(fullfile(PRM_AP_DIR,'*tcat.*.ap.bin'));

if isempty(d)
    error('could not find tcat file.')
end

%Command line loop for over strike - remember it can only zero out one time
% window at a time 
ap_path=fullfile(d.folder,d.name);
%Needs to be in overstrike Dir
system(cd(PRM_OVERSTRIKE_DIR))

for ii=1:length(stim_times)
    cmd=sprintf('OverStrike -file=%s  -secs=%d,%d',ap_path,zero_times(ii,1),zero_times(ii,2));
    [status,cmdout] = system(cmd,'-echo');
end

%Remove Cage Stim Times if needed
if ~isempty(PRM_CAGE_TIMES)
    cmd=sprintf('OverStrike -file=%s  -secs=%d,%d',...
        ap_path,PRM_CAGE_TIMES(1),PRM_CAGE_TIMES(2))
    [status,cmdout] = system(cmd,'-echo');
end

%Zero all bad channels
if ~isempty(PRM_BAD_CHANNEL0_LIST)
    lst = sprintf('%d,', PRM_BAD_CHANNEL0_LIST);
    lst(end) = [];
    chn_cmd = sprintf(' -chnexcl=%s', lst');
    cmd=sprintf('OverStrike -file=%s  -secs=0,1E9',ap_path);
    cmd_full=[cmd chn_cmd];
    [status,cmdout] = system(cmd_full,'-echo');

    %Also remove it from the LFP file
    lf_cmd_full=strrep(cmd_full,'.ap.bin','.lf.bin');
    [status,cmdout] = system(cmd_full,'-echo');
end



fprintf('Overstrike Complete for %s',d.name)
end