%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for generating all pulse trains for experiment.
% This is for EXCLUSIVELY generating bursts that repeat rather than varying
% all the time.
%
% Cowen 2022
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exp 1: All parameters are saved in the EXP structure.
% Where to save
% all of the stim files...
clearvars
close all
EXP.out_dir = fullfile(Git_dir,'DANA\Procedures_and_Experimental_Protocols\For_Experiment_On_11_22_2022b\Stimulation_timing_files');
EXP.stim_train_duration_sec = 10; % NOTE: Hill et al showed an effect after 10s of stim at 20Hz
EXP.stim_mean_freq = 20.2; %
EXP.stim_pulse_dur_sec = 0.005; % indicates the duration of each stim pulse. This is important for ensruing that stim pulses do not overlap with the onset of a CV pulse. To address this, I just subtract this value from the CV onset times so that it is addressed automatically.
EXP.stim_min_isi_sec = 0.005; % minimum interval between adjacent stimulations.
EXP.ISI_distribution = 'optimize'; % will use optimization to find

% EXP.stim_target_LV = [.5 1.0 1.5]; % If 0 then fixed interval stim.
% EXP.stim_target_LV = [.2:.01:1.5]; % If 0 then fixed interval stim.
EXP.stim_target_LV = [1.1:.005:1.4]; % If 0 then fixed interval stim.
EXP.n_pulses = 6; % in a burst

EXP.CV_pulse_duration_sec =  0.0088; % the typical duration.
EXP.CV_freq = 5; % once every 200 ms.
CV_train = 0:1/EXP.CV_freq:EXP.stim_train_duration_sec;
CV_train = CV_train(:);
blank_intervals = [CV_train(:)-EXP.stim_pulse_dur_sec CV_train(:) + EXP.CV_pulse_duration_sec];

%%%%%%%%%%%%%%%% CORE FUNCTION
INFO = [];
for iT = 1:length(EXP.stim_target_LV)
    INFO{iT}.blank_intervals = blank_intervals;
    INFO{iT}.orig_timestamps_sec = DA_generate_burst_pulse_trains_of_given_LV(EXP.stim_mean_freq, EXP.stim_target_LV(iT), EXP.n_pulses, EXP.stim_train_duration_sec + 24);
%     [~,LV(iT)] = LocalVariance( diff(INFO{iT}.orig_timestamps_sec) );

    %     INFO{iT}.orig_timestamps_sec = (INFO{iT}.timestamps/INFO{iT}.timestamps(end))*EXP.stim_train_duration_sec;
    INFO{iT}.orig_stats = DANA_stim_sequence_stats_and_IFR(INFO{iT}.orig_timestamps_sec);
    
    % The following might screw stuff up... Need to compare timestamps
    % before and after. This really lowers the LV. By a lot.

    % If a stim timestamp falls within an interval, move it back to right
    % before or right after. If a preceding pulse is within 10ms of each
    % other, move them both back.
    INFO{iT}.clean_timestamps_sec = INFO{iT}.orig_timestamps_sec;
%     INFO{iT}.clean_timestamps_sec = INFO{iT}.clean_timestamps_sec - (INFO{iT}.clean_timestamps_sec(1) - blank_intervals(1,2)+.002)
    BIX = false(size(INFO{iT}.clean_timestamps_sec));
    for iR = 1:Rows(blank_intervals)
        [tr,IX] = Restrict(INFO{iT}.clean_timestamps_sec,blank_intervals(iR,:));
        BIX = BIX | IX;
%         if ~isempty(tr)        % shift them all back to the start of the interval
%             shift = tr(end)-blank_intervals(iR,1);
%             INFO{iT}.cleaned_timestamps_sec(IX) = INFO{iT}.cleaned_timestamps_sec(IX)-shift-.001;
%             disp('shfty')
%         end
    end
    INFO{iT}.clean_timestamps_sec = INFO{iT}.clean_timestamps_sec(~BIX);
    INFO{iT}.clean_timestamps_sec(INFO{iT}.clean_timestamps_sec > EXP.stim_train_duration_sec) = [];
     figure
     plot( INFO{iT}.orig_timestamps_sec, ones(size( INFO{iT}.orig_timestamps_sec)),'o')
     hold on
     plot( INFO{iT}.clean_timestamps_sec, ones(size( INFO{iT}.clean_timestamps_sec)),'r+')
%      plot( INFO{iT}.clean_timestamps_sec(~BIX), ones(size( INFO{iT}.clean_timestamps_sec(~BIX))),'r+')
     plot( blank_intervals(:,1), ones(size( blank_intervals(:,1))),'g>')
     plot( blank_intervals(:,2), ones(size( blank_intervals(:,2))),'c<')

    INFO{iT}.clean_stats = DANA_stim_sequence_stats_and_IFR(INFO{iT}.clean_timestamps_sec);

end
%
%%%%%%%%%%%%%%%%
close all
exp_fname = sprintf('EXP_stim_fixedburst_%d_pulse_%0.0fHz_for_%0.0fsec',EXP.n_pulses,EXP.stim_mean_freq,EXP.stim_train_duration_sec);
if ~exist(EXP.out_dir, 'dir')
    mkdir(EXP.out_dir)
end

save(fullfile(EXP.out_dir,exp_fname),"INFO")
PLOT_IT = false;
% If you like these, then save them
for iT = 1:length(INFO)
    fname = sprintf('EXP_stim_fixedburst_%d_pulse_%0.0fHz_for_%0.0fsec_LV_%0.3f',...
        EXP.n_pulses, INFO{iT}.clean_stats.mean_fr,EXP.stim_train_duration_sec,INFO{iT}.clean_stats.lv);
    % PLots and SAVES the data.
    DANA_plot_and_save_stim_train( INFO{iT}, fullfile(EXP.out_dir,fname), PLOT_IT);
end
