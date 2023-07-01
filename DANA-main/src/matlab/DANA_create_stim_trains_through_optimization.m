function [OUT] = DANA_create_stim_trains_through_optimization(EXP)
% function [OUTCV_train, stim_train, lv, blank_intervals] = DANA_create_stim_trains_through_optimization(EXP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO: Finish this so that it nicely creates stim trains with no overlap.
% Should not be hard. Just some more time.
if nargin == 0
    EXP.out_dir = 'C:\Users\Stephen Cowen\Documents\GitHub\DANA\Procedures_and_Experimental_Protocols\For_Experiment_On_1_28_2022\Stimulation_timing_files';
    %     EXP.stim_train_duration_sec = 20; %
    EXP.stim_train_duration_sec = 10; %
    EXP.stim_mean_freq = 20; %
    EXP.stim_pulse_dur_sec = 0.005; % indicates the duration of each stim pulse. This is important for ensruing that stim pulses do not overlap with the onset of a CV pulse. To address this, I just subtract this value from the CV onset times so that it is addressed automatically.
    EXP.stim_min_isi_sec = .005; % minimum interval between adjacent stimulations.
    EXP.stim_target_LV =  [0 .5, 1.0, 1.4]; % If 0 then fixed interval stim.
    EXP.CV_pulse_duration_sec =  0.0088; % the typical duration.
    EXP.CV_freq = 5; % once every 200 ms.
    EXP.ISI_distribution = 'optimization';
end
%%
n_spikes =  ceil(EXP.stim_mean_freq * EXP.stim_train_duration_sec);
% To convert the above to a specific frequency or interval, you can just
% multiply the ISI vector by a factor so that the sum of the ISIs comes out
% to EXP.stim_train_duration_sec
mean_stim_ISI_sec = 1/EXP.stim_mean_freq;
% start with a random set of n spikes.
% start_ISIs = abs(randn(n_spikes,1));
start_ISIs = repmat(mean_stim_ISI_sec,n_spikes-1,1);
% NOTE: we can always change a firing sequence post-hoc to any rate we want
% through mult/div.
INFO = [];
for iT = 1:length(EXP.stim_target_LV )
    [~, INFO{iT}] = Spike_train_of_given_local_variance_using_optimization(start_ISIs,  EXP.stim_target_LV(iT));
    % Convert to the actual timestamps used in the experiment.
    INFO{iT}.orig_timestamps_sec = (INFO{iT}.timestamps/INFO{iT}.timestamps(end))*EXP.stim_train_duration_sec;
    INFO{iT}.orig_stats = DANA_stim_sequence_stats_and_IFR(INFO{iT}.orig_timestamps_sec);
    % The following might screw stuff up... Need to compare timestamps
    % before and after. Tis really lowers the LV. By a lot.
    [INFO{iT}.clean_timestamps_sec,INFO{iT}.CV_train,INFO{iT}.blank_intervals] = ...
        DANA_eliminate_overlap_with_CV_pulses(INFO{iT}.orig_timestamps_sec,EXP.CV_freq,...
        EXP.CV_pulse_duration_sec, EXP.stim_pulse_dur_sec, EXP.stim_min_isi_sec);
    INFO{iT}.clean_stats = DANA_stim_sequence_stats_and_IFR(INFO{iT}.clean_timestamps_sec);

    DANA_plot_stim_train(INFO{iT});
end
OUT.NOTE = 'Structure containing key stim times and data for the DANA experiment.';
OUT.EXP = EXP;
OUT.TRIAL_INFO = INFO; % OUT should now have all you need for each trial in timestmap_sec.
