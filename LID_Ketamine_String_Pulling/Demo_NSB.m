%%%%%%% Demo script NS&B %%%%%%%%
%%% Abhi 2023 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local field potential analysis; data collected from a Dyskinetic animal 
% Motor cortex and dorsomedial striatum regions from both hemispheres

%****** Make sure to run the script with the session data in the current
%folder ******
%% Initiallize variables 
fqs = 1:.5:170; % a vector of frequencies to look at using power spectral density or spectrogram plots

%% Load stuff
% load basic files: like injection times: E, depths of the tetrodes:DEPTHS
% META data that has some important bit to uV coversion info
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();

% lfp .mat files: The unlesioned i.e. left hem is a good control so load
% left and right hem lfp and M1 and Str regions
LFPfile{1} = find_files(fullfile('Left_M1_best_nonreref_ratio*.mat'));
LFPfile{2} = find_files(fullfile('Right_M1_best_nonreref_ratio*.mat'));
LFPfile{3} = find_files(fullfile('Left_STR_best_nonreref_ratio*.mat'));
LFPfile{4} = find_files(fullfile('Right_STR_best_nonreref_ratio*.mat'));

% determine what events are available for alignment.
intervals_around_evt_min = [];
inj1_start_end_uS = [];
inj2_start_end_uS = [];

if any(E.EventID == 'KetInjectionStart') && any(E.EventID == 'SalineInjectionStart')
        inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
        inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80];
        big_peri_event_min = [-25 110];
        OUT.inj1 = 'Sal';
        OUT.inj2 = 'Ket';

elseif any(E.EventID == 'LDOPAInjectionStart')
        inj1_start_end_uS(1) = E.MinFromStart(E.EventID == 'LDOPAInjectionStart')*60*1e6;
        inj1_start_end_uS(2) = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd')*60*1e6;
        if any(E.EventID == 'SalineInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'SalineInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'SalineInjectionEnd')*60*1e6;
            OUT.inj2 = 'Sal';
        elseif any(E.EventID == 'KetInjectionStart')
            inj2_start_end_uS(1) = E.MinFromStart(E.EventID == 'KetInjectionStart')*60*1e6;
            inj2_start_end_uS(2) = E.MinFromStart(E.EventID == 'KetInjectionEnd')*60*1e6;
            OUT.inj2 = 'Ket';
        end
        intervals_around_evt_min =  [-25 -5; 2 22; 60 80];
        big_peri_event_min = [-25 180];
        OUT.inj1 = 'Ldo';
end

% Once you have determined the injections and the times around the evts 
% create start and end time points in uS to restrict the data later on
TIMES.Event1StartEndUsec(1,1) = inj1_start_end_uS(1);
TIMES.Event1StartEndUsec(1,2) = inj1_start_end_uS(2);
TIMES.Event2StartEndUsec(1,1) = inj2_start_end_uS(1);
TIMES.Event2StartEndUsec(1,2) = inj2_start_end_uS(2);
TIMES.IntervalUsec(1,1) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.IntervalUsec(1,2) = TIMES.Event1StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.IntervalUsec(2,1) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,1)*60*1e6;
TIMES.IntervalUsec(2,2) = TIMES.Event2StartEndUsec(1) + intervals_around_evt_min(1,2)*60*1e6;
TIMES.IntervalUsec(3,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,1)*60*1e6;
TIMES.IntervalUsec(3,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(2,2)*60*1e6;
TIMES.IntervalUsec(4,1) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,1)*60*1e6;
TIMES.IntervalUsec(4,2) = TIMES.Event2StartEndUsec(2) + intervals_around_evt_min(3,2)*60*1e6;


%% Now lets make some plots 
% we'll create a for loop within a for loop 
% First loop will go through the left and right hem and brain regions 
for ix = 1:length(LFPfile)
    
    % Load lfp using the filename
    filename = LFPfile{ix}; % get the filename for the current loop from the meta lfpfile cell array
    LFP = LK_Load_and_Clean_LFP_Abhi('',filename{1}); % this function converts LFP to uV and gives an aray of timestamps that will be useful
    
    % spectrogram to look at the power at diffferent frequencies
    spec_win_sec = 5; % the window or kernal in secs for the spectrogram function
    [S,w,t] = spectrogram(LFP.LFP,LFP.sFreq*spec_win_sec,LFP.sFreq*(spec_win_sec/2),fqs, LFP.sFreq);
    sFreq_of_spectrogram = 1/median(diff(t)); %getting the sampling frequency for the modified spectrogram output
    t_min = t/60; % converting the time to minutes better visualisation
    S = 10*log10(abs(S))'; % getting the real component of S (abs(S)) and log normalizing 
    
    % Do a conv filter using a hanning window to smooth the data
    han_size_sec = 40; % hanning window size in sec
    han_size_points = han_size_sec*sFreq_of_spectrogram; % converting the seconds to sample points
    Sn = conv_filter(S,hanning(han_size_points)/sum(hanning(han_size_points)));
    
    % The below does a baseline substarction and get the z score 
%     BASEIX = t_min < 10;
%     Sb = Sn-nanmean(Sn(BASEIX,:),1);
%     Sb = Sb./nanstd(Sn(BASEIX,:),[],1);
    
    figure
    subplot(2,4,1:4)
    imagesc(t_min,fqs,Sn');
    
    axis xy
    mkrs = cat2cell(E.EventID);
    plot_markers(E.MinFromStart,1)
    title(filename{1})
    ylabel('Hz')
    colorbar
    caxis(prctile(Sn(:),[1 98.6]))
    
    
    cnt = 5; %initialize this variable for use in teh next for loop
    
    % Second loop will loop through the different time intervals created
    % earlier
    for ii = 1:Rows(TIMES.IntervalUsec)
        IXI = LFP.t_uS > TIMES.IntervalUsec(ii,1) & LFP.t_uS < TIMES.IntervalUsec(ii,2); % get out only the period of interest
        
        L = []; % initialize this variable as it will change every loop
        L = LFP.LFP(IXI); % restrict the lf data to teh time period of interest using teh logical created above
        
        % create power spectral density plots
        
        subplot(2,4,cnt)
        pwelch(L,LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
        grid off
        pubify_figure_axis
        title(sprintf('Interval %d', ii))
        
        cnt = cnt + 1;
            
    end
end
   
%% Create an example of what the raw lfp looks like and the other cool stuff
% First few steps are loading and restricting data
filename = LFPfile{2}; 
LFP = LK_Load_and_Clean_LFP_Abhi('',filename{1});
ixi = LFP.t_uS > TIMES.IntervalUsec(2,1) & LFP.t_uS < TIMES.IntervalUsec(2,2); 
L = LFP.LFP(ixi);

F = SPEC_create_filters('gamma_80', LFP.sFreq); % create an 80 Hz filter
L_f = filtfilt(F{1},L); % filter signal
pow = abs(hilbert(L_f)); % get the power envelope
ang = angle(hilbert(L_f)); % convert to radians

IX = false(length(L),1);
IX(1:500) = true;

figure
plot(L(IX))
hold on
plot(L_f(IX),'LineWidth',2)
plot(pow(IX),'LineWidth',3)
plot(ang(IX),'LineWidth',3)
legend('Raw Signal', 'Filtered', 'Hilbert', 'Angle')
xlabel('Sample points (sfreq 500)')
title('Extracted features from the lfp signal')

    
    
    