function OUT = Q1_Does_LV_affect_DA_release(WCCV_dir, WCCV_file)
%% Operate on one file only. Put this function in a loop if you want to analyze multiple sessions.
% Cowen 2022.
if nargin == 0 
    WCCV_dir = 'C:\Users\Stephen Cowen\Documents\GitHub\DANA\WCCV\';
    WCCV_file = '20220401 data.xlsx';
end
WCCV_path_and_fname = fullfile(WCCV_dir,WCCV_file);
% 20220218 data.xlsx
% Load the stimulations...
[DATA, STIM] = DANA_load_wccv_excel_file(WCCV_path_and_fname);
% Sliding LV on the data...
[LV,LVR] = LocalVariance(diff(STIM(4).Stim_time_s*1000),[],[4000 1000]);
figure; plot(LVR); hold on; plot(LV)
%% T = readtable('patients.xls',...
%     'Range','C2:E6',...
%     'ReadVariableNames',false)