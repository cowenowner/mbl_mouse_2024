function [OUT] = DANA_create_stim_trains_from_ISI_distribution(EXP)
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
    EXP.ISI_distribution = 'lognormal';

end
%%
n_spikes =  ceil(EXP.stim_mean_freq * EXP.stim_train_duration_sec);
% To convert the above to a specific frequency or interval, you can just
% multiply the ISI vector by a factor so that the sum of the ISIs comes out
% to EXP.stim_train_duration_sec
mean_stim_ISI_sec = 1/EXP.stim_mean_freq;
% start with a random set of n spikes.
% start_ISIs = repmat(mean_stim_ISI_sec,n_spikes-1,1);
% NOTE: we can always change a firing sequence post-hoc to any rate we want
% through mult/div.

% For a given parameter range - generate a bunch of ISIs - choose the
% one with the closest match to the target.
switch EXP.ISI_distribution
    case 'lognormal'
        coef =  0:.01:4;
    case 'exponential'
        coef = 0:.05:3;
    case 'gamma'
        coef = 0:.005:3;
    case 'fixed_interval_doublet'
        coef = 0:.005:3;
        coef2 = 0:.005:3;
end
lv = zeros(length(coef),1);
lvc = zeros(length(coef),1);
IS = zeros(length(coef),n_spikes-1); % one less iSI than spikes
ISc = zeros(length(coef),n_spikes-1); % one less iSI than spikes
TS = zeros(length(coef),n_spikes);
TSc = zeros(length(coef),n_spikes);

for iP = 1:length(coef)
    switch EXP.ISI_distribution
        case 'lognormal'
            tmpISI = random('LogNormal',EXP.stim_mean_freq,coef(iP),n_spikes-1,1);
        case 'exponential'
            tmpISI = exprnd(EXP.stim_mean_freq,n_spikes-1,1);
            tmpISI = tmpISI.^coef(iP);
        case 'gamma'
            tmpISI = gamrnd(coef(iP),EXP.stim_mean_freq,n_spikes-1,1);
        case 'fixed_interval_doublet'
            tmpISI = gamrnd(coef(iP),EXP.stim_mean_freq,n_spikes-1,1);
    end
    a = (EXP.stim_mean_freq*n_spikes)/sum(tmpISI);
    IS(iP,:)  = tmpISI*a;
    tmpISI = IS(iP,:);
    tmpISI(tmpISI<EXP.stim_min_isi_sec) = EXP.stim_min_isi_sec;
    a = (EXP.stim_mean_freq*n_spikes)/sum(tmpISI);
    ISc(iP,:)  = tmpISI*a;

    TS(iP,:) = [0 cumsum(IS(iP,:))];
    TSc(iP,:) = [0 cumsum(ISc(iP,:))];

    [~,lv(iP)] = LocalVariance(IS(iP,:)); % Changed to the more recent LV measure 3/7/2022
    [~,lvc(iP)] = LocalVariance(ISc(iP,:));
end
figure
for iP = 1:length(coef)
    plot(TS(iP,:), repmat(lv(iP),1,length(TS(iP,:))),'.')
    hold on
    plot(TSc(iP,:), repmat(lvc(iP),1,length(TSc(iP,:))),'+')

end
axis tight
ylabel('LV')

INFO = [];
for iT = 1:length(EXP.stim_target_LV )
    ix = Closest(lvc,EXP.stim_target_LV(iT));
    ISI = ISc(ix,:)';
    orig_timestamps = TSc(ix,:)';
    % Convert to the actual timestamps used in the experiment.
    INFO{iT}.orig_timestamps_sec = (orig_timestamps/orig_timestamps(end))*EXP.stim_train_duration_sec;
    INFO{iT}.orig_stats = DANA_stim_sequence_stats_and_IFR(INFO{iT}.orig_timestamps_sec);
    
    % The following might screw stuff up... Need to compare timestamps
    % before and after. Tis really lowers the LV. By a lot.
    [INFO{iT}.clean_timestamps_sec,INFO{iT}.CV_train,INFO{iT}.blank_intervals] = ...
        DANA_eliminate_overlap_with_CV_pulses(INFO{iT}.orig_timestamps_sec,EXP.CV_freq,...
        EXP.CV_pulse_duration_sec, EXP.stim_pulse_dur_sec, EXP.stim_min_isi_sec);
    INFO{iT}.clean_stats = DANA_stim_sequence_stats_and_IFR(INFO{iT}.clean_timestamps_sec);
    %
    DANA_plot_stim_train(INFO{iT});
end
OUT.NOTE = 'Structure containing key stim times and data for the DANA experiment.';
OUT.EXP = EXP;
OUT.TRIAL_INFO = INFO; % OUT should now have all you need for each trial in timestmap_sec.

