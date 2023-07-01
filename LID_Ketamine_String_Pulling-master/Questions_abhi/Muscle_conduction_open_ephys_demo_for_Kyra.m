% Data analysis for Kyra....
%
% REQUIRES https://github.com/open-ephys/analysis-tools/blob/master/load_open_ephys_data.m
% REQUIRES https://github.com/CowenLab/CowenLib
%% Read in the data...
[data, timestamps, info] = load_open_ephys_data('E:\NerveStim_Kyra\25 NT D2_2022-02-21_11-41-16\Record Node 103\100_5.continuous');
% [dat_ch3, timestamps, info] = load_open_ephys_data('E:\NerveStim_Kyra\EPhys\2-24-2022\35 NT P2_2022-02-24_14-23-59\Record Node 103\100_25.continuous');
[dat_stim, timestamps, info] = load_open_ephys_data('E:\NerveStim_Kyra\25 NT D2_2022-02-21_11-41-16\Record Node 103\100_33.continuous');
info.header
% rereferencing is another way to reduce noise
% reref_data = data-dat_ch3;

% Sanity check: determine sampling rate.
sFreq = 1/median(diff(timestamps));
sFreq_Andy = 1/GrassStimulator_03_31_2022_Ch1.interval;
timestamps_Andy = (0:(GrassStimulator_03_31_2022_Ch1.length-1))/sFreq_Andy;

% If the data needs to be filtered (e.g., to get rid of 60Hz noise), do it
% here.
[d] = Notch_filter_cowen( sFreq, 59.6,  60.5 );
data2 = filtfilt(d,data);
% reref_data = filtfilt(d,reref_data);

% Stim TTLs getting the start and end times 
dat_stim_pos = dat_stim-(min(dat_stim));
stim_int = find_intervals([timestamps,dat_stim_pos],0.1);
stim_int = stim_int(1:length(stim_int)-1,:);
stim_int_Andy = find_intervals([timestamps_Andy',GrassStimulator_03_31_2022_Ch2.values],0.5);

% plot the data.
figure;
plot(timestamps,data)
hold on
plot(timestamps,data2)
plot_markers_simple(stim_int(:,1))
legend('orig','notch')

xlabel('seconds')
ylabel('bit values (not volts)')

figure
plot(timestamps_Andy,GrassStimulator_03_31_2022_Ch1.values)
plot_markers_simple(stim_int_Andy(:,1))

% Peri event histogram for the stim TTL end times
figure;
[M_reref,ix,x_sec] = PETH_EEG_simple([timestamps,data2],stim_int(:,1),sFreq/16,sFreq/8,sFreq,true);

figure
[M_reref_Andy,ix_Andy,x_sec_Andy] = PETH_EEG_simple([timestamps_Andy',GrassStimulator_03_31_2022_Ch1.values],stim_int_Andy(:,1),sFreq_Andy/16,sFreq_Andy/8,sFreq_Andy,true);


% Find the start time of muscle conduction post stim  
% t_0 = find(x_sec == 0)+10; %finding time 0 i.e. end of TTL pulse  
% mn_before_stim = mean(M_reref(:,521:921)')'; % finding the mean for 500 ms before stim pulse end i.e from -0.01 to -0.06 sec
% [mn,ci] = conf_int(mn_before_stim,'normal',.01); % finding the CI for the 500ms before stim
% find_below_CI = find(mean(M_reref(:,t_0:end)) < ci(2,1)); 
% start_mus_condu = x_sec(t_0+find_below_CI(1));

% Find the peak amplitude of the muscle conduction spike

% NEED VOLTS ????

% % Ginput asks for user input. Before doing this, zoom into a single clear
% % compound evoked AP...
% time = ginput(2);
% time_to_compond_A_sec = time(2,1)- time(1,1);