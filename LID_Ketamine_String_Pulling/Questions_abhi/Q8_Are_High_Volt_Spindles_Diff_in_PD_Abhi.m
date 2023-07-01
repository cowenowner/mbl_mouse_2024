function OUT = Q8_Are_High_Volt_Spindles_Diff_in_PD()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are HVS altered in PD?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_IT = true;
SES = LK_Session_Info();

OUT = [];

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();

% Load the appropriate channel combo.
% Presumes that  LK_Find_best_LFP_for_a_band('gamma_80') was run in this
% directory.
load('./Processed_Data/best_channels_gamma_80.mat','OUT')
%% Load the best non-reref.
if exist(fullfile('LFP',OUT.best_non_reref),'file')
    % load if stored locally.
    LFP = LK_Load_and_Clean_LFP('./LFP',OUT.best_non_reref);
else
    global DIRS
    lfp_dir = fullfile(DIRS.LFP_Dir,SES.rat_str,SES.session_str,'LFP');
    LFP = LK_Load_and_Clean_LFP(lfp_dir, OUT.best_non_reref);
end
% copyfile(LFP.full_file_path,'./LFP');
t_LDOPA_inj_min = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd');
t_LDOPA_effect_win = [30 50] + t_LDOPA_inj_min;
t_LDOPA_base_win = [-25 -5] + t_LDOPA_inj_min;
IXeff = LFP.t_uS >t_LDOPA_effect_win(1)*60e6 & LFP.t_uS < t_LDOPA_effect_win(2)*60e6;
IXbase = LFP.t_uS >t_LDOPA_base_win(1)*60e6 & LFP.t_uS < t_LDOPA_base_win(2)*60e6;

%% Look for spindles in the baseline period.
figure;
plot(LFP.t_uS(IXbase),LFP.LFP(IXbase));

% [x] = ginput();
spindle_times = [4891649678.38487 4903675776.05835 ;5.2993e+09 5.3092e+09];
LS = Restrict([LFP.t_uS LFP.LFP], spindle_times);
figure
%plot(LS(:,1)/1e6,LS(:,2),'.-')
plot((1:Rows(LS))/LFP.sFreq, LS(:,2),'.-')
axis tight
xlabel('sec')

figure
cwt(LS(:,2),'bump',LFP.sFreq,'FrequencyLimits',[2 55],'VoicesPerOctave',24);

figure
spectrogram(LS(:,2),LFP.sFreq,round(LFP.sFreq/3),2:.25:55,LFP.sFreq,'yaxis');
a = caxis;
caxis([-40 a(2)])

filts = SPEC_create_filters({[13 18] [10 13]}, LFP.sFreq);
figure
plot(LS(:,1)/1e6,LS(:,2))
hold on
plot(LS(:,1)/1e6,filtfilt(filts{1},LS(:,2)))
plot(LS(:,1)/1e6,filtfilt(filts{2},LS(:,2)))

axis tight
xlabel('sec')

pow1 = abs(hilbert(filtfilt(filts{1},LS(:,2))));
pow2 = abs(hilbert(filtfilt(filts{2},LS(:,2))));

plot(LS(:,1)/1e6,LS(:,2))
hold on
plot(LS(:,1)/1e6,filtfilt(filts{1},LS(:,2)))
plot(LS(:,1)/1e6,pow1)
plot(LS(:,1)/1e6,pow2)

axis tight
xlabel('sec')


%% Find spindles for the entire datatset
pow_sp = abs(hilbert(filtfilt(filts{1},LFP.LFP)));
pow_cltr = abs(hilbert(filtfilt(filts{1},LFP.LFP)));

se = find_intervals([LFP.t_uS,pow_sp],49);

figure
plot(LFP.t_uS/60e6, LFP.LFP);
yyaxis right
plot(LFP.t_uS/60e6, conv_filter(pow_sp-pow_cltr,hanning(LFP.sFreq)));

plot(se(:,1)/60e6,ones(size(se(:,1))),'g>')
plot(se(:,2)/60e6,ones(size(se(:,1))),'r<')



