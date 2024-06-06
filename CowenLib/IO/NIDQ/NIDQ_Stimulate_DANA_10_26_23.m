%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DANA stimulations trains for neuropixels recordings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2022
% Edited by AV 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%% SINGLE PULSE STIM %%%%%%%%%%%%%%%%%%
%Save Stim file path and name input var
Data_path='C:\SGL_DATA\Acute_DANA\Rat482';
Params_fileName='params_from_stim_exp_01_12_24_Rat482.mat';

Stim_path = 'C:\Users\CowenLab\Documents\GitHub\DANA\Stimulation_trains\Fixed_burst_stimulation_trains\Single_point_stimulation_trains';

%Stimulation Parameters
PRM.biphasic = true; % set to false if using old WPI
PRM.biphasic_both_positive = false; % Both postive ONLY if using old WPI.
PRM.time_between_phases_ms = 0; % need some time - say .05ms if using the old orange WPI. 

% PULSE width CHANGE depending on stim location
PRM.pulse_duration_ms = 2; % mPFC
% PRM.pulse_duration_ms = 4; % MFB

PRM.time_between_bursts_sec = 280; %s

PRM.PLOT_IT=false;

%Amplitude CHANGE depending of stim location
PRM.output_current= 6.83; % This is equal to 600uA based on the translation - mPFC
% PRM.output_current= 3.48; % This is equal to 300uA based on the translation - MFB


PRM.n_trials = 40;

% Get all the LV stim times
[~,stim_txt_files] = find_files(fullfile(Stim_path,'*.txt'));
for ist = 1:length(stim_txt_files)
    fileID = fopen(stim_txt_files{ist},'r'); % opening text files
    pulse_sequence_sec{ist} = fscanf(fileID,'%f'); % storing the timestamps for each LV 
end
PRM.stim_txt_filenames = stim_txt_files;

% generate a randomised stimulation order of LVs 
Stim_by_trials = [ones(1,5),ones(1,5)*2,ones(1,5)*3,ones(1,5)*4,ones(1,5)*5 ...
    ones(1,5)*6,ones(1,5)*7,ones(1,5)*8]; % 8 stim patterns of 5 trials each
PRM.LV_stim_order_rand = Stim_by_trials(randperm(length(Stim_by_trials))); % randomizing the stim order

for ii = 1:length(PRM.LV_stim_order_rand)
    PRM.LV_stim_order_rand_txt_filenames{ii} = PRM.stim_txt_filenames{PRM.LV_stim_order_rand(ii)}; % LV filenames in thej randomized order
end

%SAVE PRM FILE
save(fullfile(Data_path,Params_fileName))
%%
%%%%%%%%%%%%%%%%%% ONLY RUN FROM HERE %%%%%%%%%%%%%%%%%%
tic;
for iT = 1:PRM.n_trials
    NIDQ_Stimulate_IO(pulse_sequence_sec{PRM.LV_stim_order_rand(iT)}, 'output_voltage', PRM.output_current, ...
        'individual_pulse_duration_sec', PRM.pulse_duration_ms/1000, ...
        'biphasic_pulses', PRM.biphasic, ...
        'biphasic_interval_sec', PRM.time_between_phases_ms/1000, ...
        'biphasic_both_positive', PRM.biphasic_both_positive,...
        'PLOT_IT',PRM.PLOT_IT);
    fprintf('%d - %s ',iT,PRM.LV_stim_order_rand_txt_filenames{iT});
    pause(PRM.time_between_bursts_sec)
end
toc;

