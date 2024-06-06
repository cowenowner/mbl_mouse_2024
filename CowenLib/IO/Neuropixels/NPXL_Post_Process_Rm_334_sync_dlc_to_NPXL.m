function NPXL_Post_Process_Rm_334_sync_dlc_to_NPXL(NPXL_data_dir, varargin)
% Assumes that deeplab cut has been run and that the INTAN data has been
% post-processed. The deep lab cut directory (DLC) should be in the directory 
% ABOVE the NPXL data directory.
%  e.g.: 'Z:\Data\String_Pull_312a\RECORDING SESSIONS\Rat 333\3\Vids_DLC'
%
% INPUT:
%  NPXL_data_dir - the root directory where the INTAN data was stored
%     for example: 'Z:\Data\String_Pull_312a\RECORDING_SESSIONS\Rat_333\3\Rec_210818_101300'
%  event_name - the name of the sub-structure in the event file that
%     corresponds for example 'top_camera_frame_ID' from teh EVT structure.
%  dlc_csv_file_path - full path to the dlc file corresponding to the
%     event_name
%
% optional: PRM - a strcture that specifices the specific cameras that will
% be processed (and the EVT variables).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRM.camera_list = {'Top_Cam' 'Front_Cam'}; % names should correspond to the folder names in /DLC that correspond to each camera.
PRM.camera_unique_code = {'t' 'f'}; % names should correspond to the folder names in /DLC that correspond to each camera.
PRM.EVT_file_sync_variable = {'topcam' 'frontcam'};
PRM.DLC_folder = 'Vids_DLC';
% PRM.camera_list = {'Front_Cam'};
% PRM.EVT_file_sync_variable = {'front_camera_frame_uS'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Extract_varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root_dir = fileparts(NPXL_data_dir);
dlc_csv_file_path = fullfile(NPXL_data_dir,PRM.DLC_folder);

% load the events...
if ~exist(fullfile(NPXL_data_dir,'Events.mat'),'file')
    fullfile(NPXL_data_dir,'Events.mat')
    error('no EVT.mat file')
end
load(fullfile(NPXL_data_dir,'Events.mat'),'EVT')
% Use the '_fixed' camera variable if it exists
for ii = 1:length(PRM.EVT_file_sync_variable)
    vbl = [PRM.EVT_file_sync_variable{ii} '_fixed'];
    if any(contains(fieldnames(EVT),vbl))
        PRM.EVT_file_sync_variable{ii} = vbl;
    end
end
%% Sync each camera in PRM.camera_list
POS = [];
POS_info = [];

for iC = 1:length(PRM.camera_list)
    cam = PRM.camera_list{iC};
    cam_evt_vbl = PRM.EVT_file_sync_variable{iC};
    d = dir(fullfile(dlc_csv_file_path,['*' PRM.camera_unique_code{iC} 'DLC*.csv']));
    if isempty(d)
        error(['No .csv file from DeepLabCut found for ' cam ])
    elseif length(d) > 1
        d
        error('Multiple .csv files found')
    end
    new_dlc_csv_file_path = fullfile(dlc_csv_file_path,d(1).name);
    
    if exist(dlc_csv_file_path,'file')
        disp(['Found ' cam])
        % The following also creates the relevant .mat files in the .csv
        % directory.
        [POS{iC}, POS_info{iC}] = SPL_Sync_DeepLabCut_csv_to_NPXL(new_dlc_csv_file_path, EVT.(cam_evt_vbl).time_sec*1e6, cam);
    end
end
