function [CV_train, stim_train, lv, blank_intervals] = DANA_create_CV_and_stim_trains(duration_sec, CV_freq, CV_duration_s, stim_params, out_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [CV_train, stim_train, lv] = DANA_create_CV_and_stim_trains(duration_sec, CV_freq, CV_duration_s, stim_params, out_dir)
%
% INPUT:
%  duration_sec = duration of stim pulse train.
%  CV_freq = the frequency of delivering the CV pulses (typically 10 or 5
%  Hz)
%  CV_duration_s = duration of a single CV pulse (typically 8ms)
%  
%  stim_params = a structure with the following form...
%     stim_params.mean_freq = 40; % frequency of the stimulation.
%     stim_params.wts = [0 1 0];
%     stim_params.exponent = 2;
%     stim_params.n_sec_stim = duration_sec;
%     stim_params.stim_pulse_dur_sec = 0.005; % 
%     stim_params.min_isi_sec = .005;
%
% out_dir = where to dump the text files. leave empty if you do not want to
% do this.
%
% OUTPUT:
%  CV_train = the START time of each CV pulse (typically 8ms long, but in CV_duration_s).
%  stim_train = the START time of each stim pulse.
%  lv = the local variance measurement for the train.
%
% NOTE: the CV train won't actually be used by WCCV as this is generated
% at a very low level, but it help confirm that the stim pulses are all
% good. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_IT = false;
if nargin == 0
    % for testing
    out_dir = 'C:\Temp\'; % or 'C:\Users\Stephen Cowen\Documents\GitHub\DANA\Procedures_and_Experimental_Protocols\For_Experiment_On_1_28_2022\Stimulation_timing_files';
    duration_sec = 20; % duration of the pulse train.
    CV_duration_s = 0.0088; % the typical duration.
    CV_freq = 5; % once every 200 ms.
    stim_params.mean_freq = 20; % frequency of the stimulation.
    stim_params.wts = [1 1 0]; % used to generate LV values
    stim_params.exponent = 3; % this indicates the decay or tail of the exponential distribution used for inter-stim interval generation. The bigger, the more bursty.
    stim_params.stim_pulse_dur_sec = 0.005; % indicates the duration of each stim pulse. This is important for ensruing that stim pulses do not overlap with the onset of a CV pulse. To address this, I just subtract this value from the CV onset times so that it is addressed automatically.
    stim_params.n_sec_stim = duration_sec;
    stim_params.min_isi_s = .005;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a variable pulse train.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stim_train, lv_before] = Generate_spikes_with_different_burst_stats(stim_params.wts, ...
    stim_params.exponent, stim_params.mean_freq, stim_params.n_sec_stim, stim_params.min_isi_s);

stim_train = stim_train(:);
% 
[stim_train, CV_train, blank_intervals] = DANA_eliminate_overlap_with_CV_pulses(stim_train,CV_freq, CV_duration_s, stim_params.stim_pulse_dur_sec, stim_params.min_isi_s);
stim_hz = 1/median(diff(stim_train));
lv = LocalVariance(diff(stim_train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write text files that LabView can read (hopefully).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(out_dir)
    % write the text files
    fileID = fopen(fullfile(out_dir,'CV_train.txt'),'w');
    fprintf(fileID,'%f\n',CV_train);
    fclose(fileID);
    
    fname = sprintf('stim_train_Hz%1.0f_LV%1.2f.txt',stim_hz,lv);
    fileID = fopen(fullfile(out_dir,fname),'w');
    fprintf(fileID,'%f\n',stim_train);
    fclose(fileID);
end

if PLOT_IT || nargout == 0
    figure
    plot(CV_train, zeros(size(CV_train)),'k.')
    hold on
    plot(blank_intervals(:,1), zeros(size(blank_intervals(:,1))),'g>')
    plot(blank_intervals(:,2), zeros(size(blank_intervals(:,2))),'r<')
    plot(stim_train, zeros(size(stim_train)),'c+')
    title(sprintf('LV = %f',lv))
    axis tight
end