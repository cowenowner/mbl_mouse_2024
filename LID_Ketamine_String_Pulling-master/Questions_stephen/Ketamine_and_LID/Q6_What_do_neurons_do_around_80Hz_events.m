function OUT = Q6_What_do_neurons_do_around_80Hz_events()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do neurons correlate activity with 80 Hz events?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_IT = true;
SES = LK_Session_Info();
OUT = [];
bin_and_inc_size_msec = 2;
time_before_msec = 150;
time_after_msec = 150;

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();

load('./Processed_Data/gamma_80_events.mat','GE');

figure
histogram( GE.GAM.dur_s)

t_LDOPA_inj_min = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd');
t_LDOPA_effect_win_min = [30 50] + t_LDOPA_inj_min;
TSr = Restrict(TS,t_LDOPA_effect_win_min*60e6);
%%
% e = GE.GAM.t_maxpow_uS;
e = GE.GAM.max_peak_uS;
% e = GE.GAM.st_uS;
%  e = GE.GAM.start_peak_uS;
e_uS = Restrict(e,t_LDOPA_effect_win_min*60e6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iN = 1:length(TSr)
    figure
    %     [M, x_axis, A_msec, h ] = PETH_raster(TS{iN}/100,e/100,bin_and_inc_size_msec,time_before_msec,time_after_msec,{'mean'  'raster' });
    [M, x_axis, A_msec, h ] = PETH_raster(TSr{iN}/100,e_uS/100,bin_and_inc_size_msec,time_before_msec,time_after_msec,{'mean'  'raster' });
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iN = 1:length(TSr)
    figure
    [OUT.S(iN).Spike_PSD,OUT.S(iN).Spike_PSD_x_Hz] = Spike_psd(TSr{iN}/1000,2,1:.25:95,16);
end
