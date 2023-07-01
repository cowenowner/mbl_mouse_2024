%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for generating all pulse trains for experiment.
%
% Cowen 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp 1: All parameters are saved in the EXP structure.
% Where to save all of the stim files...
EXP.out_dir = 'C:\Users\Stephen Cowen\Documents\GitHub\DANA\Procedures_and_Experimental_Protocols\For_Experiment_On_2_11_2022\Stimulation_timing_files';
EXP.stim_train_duration_sec = 10; % NOTE: Hill et al showed an effect after 10s of stim at 20Hz
EXP.stim_mean_freq = 20; %
EXP.stim_pulse_dur_sec = 0.005; % indicates the duration of each stim pulse. This is important for ensruing that stim pulses do not overlap with the onset of a CV pulse. To address this, I just subtract this value from the CV onset times so that it is addressed automatically.
EXP.stim_min_isi_sec = 0.005; % minimum interval between adjacent stimulations.
EXP.ISI_distribution = 'optimize'; % will use optimization to find
% EXP.ISI_distribution = 'lognormal';
% EXP.ISI_distribution = 'exponential';
EXP.ISI_distribution = 'gamma';

EXP.stim_target_LV = [.5 1.0 1.5]; % If 0 then fixed interval stim.
%  EXP.stim_target_LV = [2.4 ]; % If 0 then fixed interval stim.
%   EXP.stim_target_LV = [0]; % If 0 then fixed interval stim.

EXP.CV_pulse_duration_sec =  0.0088; % the typical duration.
EXP.CV_freq = 5; % once every 200 ms.
%%%%%%%%%%%%%%%% CORE FUNCTION
switch EXP.ISI_distribution
    case 'optimize'
        FULL_EXP = DANA_create_stim_trains_through_optimization(EXP);
    otherwise
        FULL_EXP = DANA_create_stim_trains_from_ISI_distribution(EXP);
end
%%%%%%%%%%%%%%%%
exp_fname = sprintf('EXP_stim_%0.0fHz_for_%0.0fsec_%s',EXP.stim_mean_freq,EXP.stim_train_duration_sec,FULL_EXP.EXP.ISI_distribution);
save(fullfile(EXP.out_dir,exp_fname),"FULL_EXP")

% If you like these, then save them
for iT = 1:length(FULL_EXP.TRIAL_INFO)
    fname = sprintf('EXP_stim_%0.0fHz_for_%0.0fsec_LV_%0.2f_%s',...
        EXP.stim_mean_freq,EXP.stim_train_duration_sec,FULL_EXP.TRIAL_INFO{iT}.clean_stats.lv,FULL_EXP.EXP.ISI_distribution);
    % PLots and SAVES the data.
    DANA_plot_and_save_stim_train( FULL_EXP.TRIAL_INFO{iT}, fullfile(EXP.out_dir,fname));
end
