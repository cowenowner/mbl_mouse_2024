function OUT = Q3_What_is_LDOPA_80Hz()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does LDOPA Create 80-Hz oscillations?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Presume that this is an LDOPA day.

% global DIRS
PLOT_IT = true;
DO_CFC = false;
SES = LK_Session_Info();

OUT = [];
fqs = 1:.5:150;
OSC = Oscillation_frequencies;

[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();

POS = LK_Load_and_Clean_POS;
IMU = LK_Load_and_Process_IMU;

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

filts = SPEC_create_filters({'gamma_80' [65 75] [3 6.5]}, LFP.sFreq);
t_LDOPA_inj_min = E.MinFromStart(E.EventID == 'LDOPAInjectionEnd');
t_LDOPA_effect_win = [30 50] + t_LDOPA_inj_min;
t_LDOPA_base_win = [-25 -5] + t_LDOPA_inj_min;
IXeff = LFP.t_uS >t_LDOPA_effect_win(1)*60e6 & LFP.t_uS < t_LDOPA_effect_win(2)*60e6;
IXbase = LFP.t_uS >t_LDOPA_base_win(1)*60e6 & LFP.t_uS < t_LDOPA_base_win(2)*60e6;
IXeffPOS = POS.Time_uS >t_LDOPA_effect_win(1)*60e6 & POS.Time_uS  < t_LDOPA_effect_win(2)*60e6;
IXeffIMU = IMU.t_uS >t_LDOPA_effect_win(1)*60e6 & IMU.t_uS  < t_LDOPA_effect_win(2)*60e6;

TSr = Restrict(TS,t_LDOPA_effect_win*60e6);
binsize_ms = 5;
% [Q,Qt] = Bin_ts_array(TSr,binsize_ms*1e3);
[Q,Qt] = histcounts_cowen(TSr,'binsize',binsize_ms*1e3);
Qt = Qt(1:end-1) + binsize_ms*1e3/2;
Q = conv_filter(Q,hanning(10));
[pcs,sc,lat] = pca(Q);

% Create an index of 80 power.
L = zeros(length(LFP.LFP),length(filts));
for ii = 1:length(filts)
    L(:,ii) = filtfilt(filts{ii},LFP.LFP);
end
pow80 = abs(hilbert(L(:,1))) - abs(hilbert(L(:,2)))./ abs(hilbert(L(:,1))) + abs(hilbert(L(:,2)));
pow80 = (abs(hilbert(L(:,1))) - abs(hilbert(L(:,2))))./ (abs(hilbert(L(:,1))) + abs(hilbert(L(:,2))));

pow_low_th = abs(hilbert(L(:,3)));
%%
if PLOT_IT
    figure
    subplot(1,3,1)
    pwelch(LFP.LFP(IXeff),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
    grid off
    pubify_figure_axis
    title('80Hz during peak win')
    subplot(1,3,2)
    pwelch(LFP.LFP(IXbase),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
    grid off
    pubify_figure_axis
    title('prior to inj')
    subplot(1,3,3)
    cla
    pwelch(pow80(IXeff),LFP.sFreq,LFP.sFreq/2,1:.2:40,LFP.sFreq)
    hold on
    pwelch(pow80(IXbase),LFP.sFreq,LFP.sFreq/2,1:.2:40,LFP.sFreq)
    
    grid off
    pubify_figure_axis
    h = gca;
    h.Children(1).Color = [0 0 0 ];
    title('psd of power changes at 80hz')
    
end

figure
plot(LFP.t_uS(IXeff)/1e6, L(IXeff,1))
hold on
plot(LFP.t_uS(IXeff)/1e6, pow80(IXeff))
plot(LFP.t_uS(IXeff)/1e6, pow_low_th(IXeff))

plot(LFP.t_uS(IXeff)/1e6, LFP.LFP(IXeff))
legend('f80','p80','pthet','orig');

figure
plot(LFP.t_uS(IXeff)/1e6, L(IXeff,1))
hold on
plot(LFP.t_uS(IXeff)/1e6, pow80(IXeff))
yyaxis right
plot(Qt/1e6,mean(Q,2))
plot(Qt/1e6,sc(:,1),'r-')
plot(Qt/1e6,sc(:,2),'b-')

%% plot position data
figure
plot(POS.Time_uS(IXeffPOS)/60e6,POS.speed(IXeffPOS));
ylabel('pix/sec')
yyaxis right
plot(IMU.t_uS(IXeffIMU)/60e6,IMU.speed(IXeffIMU));
axis tight
xlabel('min')
legend('pos','imu')
ylabel('a.u.')
%% Instantaneous frequency
%  [fq_ga,t_sec] = instfreq(L(IXeff,1),LFP.sFreq,'Method','hilbert'); % works the SAME as below.
%  fq_th = instfreq(L(IXeff,3),LFP.sFreq,'Method','hilbert');
[fq_ga,t_sec] = instfreq(LFP.LFP(IXeff),LFP.sFreq,'FrequencyLimits',OSC.gamma_80);
[fq_th] = instfreq(LFP.LFP(IXeff),LFP.sFreq,'FrequencyLimits',[3 6.5]);

C = kmeans([fq_th fq_ga], 2);

figure
subplot(2,3,1)
scatter(fq_th,fq_ga,4,C)
xlabel('theta');ylabel('80Hz gamma');
lsline
mdl = fitlm(fq_th,fq_ga);
title(sprintf('p= %1.4f, R2adj= %1.2f', mdl.Coefficients.pValue(2),mdl.Rsquared.Adjusted));

subplot(2,3,2)
histogram(fq_th,40)
title('low theta fq')
subplot(2,3,3)
histogram(fq_ga,40)
title('gamma fq')
subplot(2,3,4:6)
plot(t_sec/60,fq_th)
yyaxis right
plot(t_sec/60,fq_ga)
xlabel('min')
axis tight

%% Do the same for power...
p1 = pow_low_th(IXeff);
p2 = pow80(IXeff);
p1(abs(p1)>150)=nan;
p2(abs(p2)>150)=nan;

C = kmeans([p1(1:20:end) p2(1:20:end)], 2);
figure
subplot(2,3,1)
scatter(p1(1:20:end),p2(1:20:end),4,C)
xlabel('pw theta');ylabel('pw 80Hz gamma');
lsline
mdl = fitlm(p1(1:20:end),p2(1:20:end));
title(sprintf('p= %1.4f, R2adj= %1.2f', mdl.Coefficients.pValue(2),mdl.Rsquared.Adjusted));

subplot(2,3,2)
histogram(p1(1:20:end),40)
title('low theta pow')
subplot(2,3,3)
histogram(p2(1:20:end),40)
title('gamma pow')
subplot(2,3,4:6)
plot(p1(1:20:end))
yyaxis right
plot(p2(1:20:end))
xlabel('sample')
axis tight


% hold on
% plot(fq_th)
% plot(fq_ga2)
% plot(fq_th2)
%

%% do spike stats correlate with lfp parameters?
Qmn_int = interp1(Qt, mean(Q,2), LFP.t_uS(IXeff));
Qsc1_int = interp1(Qt, sc(:,1), LFP.t_uS(IXeff));
Qsc2_int = interp1(Qt, sc(:,2), LFP.t_uS(IXeff));

[R,P] = corr([pow_low_th(IXeff) L(IXeff,1) pow80(IXeff) Qmn_int Qsc1_int Qsc2_int [0;diff(Qsc2_int)]],'Type','Spearman','Rows','complete');
figure
labs = {'lowth','80hz','pow80','mnQ','sc1','sc2','dsc2'};
imagesc(1:Cols(R),1:Rows(R),R)
set(gca,'XTickLabel',labs)
set(gca,'YTickLabel',labs)
caxis([0.01 .2])
axis xy
colorbar
title('correlation betwn vbls')


% Wavelet.

% cwt(LFP.LFP(IXeff),'bump',LFP.sFreq,'FrequencyLimits',[2 110])
% cwt(LFP.LFP(IXeff),'morse',LFP.sFreq,'FrequencyLimits',[2 110])
[cfs_tmp,frq_tmp] = cwt(LFP.LFP(IXeff),'bump',LFP.sFreq,'FrequencyLimits',[2 110],'VoicesPerOctave',24);
out_fqs = 2:1:110;
[cfs,~,A] = cwt_fix(cfs_tmp,frq_tmp,out_fqs);
% cfs = cfs(end:-1:1,:); frq = frq(end:-1:1);
THIX = out_fqs > 3.4 & out_fqs < 10;
GAIX = out_fqs > 74 & out_fqs < 94;
d_cfs = decimate_matrix(abs(cfs),10);
new_x_uS = linspace(min(LFP.t_uS(IXeff)),max(LFP.t_uS(IXeff)),Cols(d_cfs));
figure
subplot(2,1,1)

imagesc(new_x_uS(:)/60e6,out_fqs,d_cfs')
axis xy
caxis(prctile(d_cfs(1:10:end),[4 99]))
axis tight
% shading flat
xlabel('Time (min)')
ylabel('Frequency (Hz)')
colorbar

subplot(2,1,2)
spectrogram(LFP.LFP(IXeff),LFP.sFreq*1,LFP.sFreq*.5,2:.5:110,LFP.sFreq,'yaxis');
c = caxis;
caxis([-20 c(end)])

figure
subplot(2,1,1)
pwth = nanmean(cfs(:,THIX),2);
pwga = nanmean(cfs(:,GAIX),2);
plot(pwth(1:20:end),pwga(1:20:end),'.','MarkerSize',.5)
xlabel('theta pow');ylabel('gam pow');
lsline
subplot(2,1,2)
p = spectrogram(LFP.LFP(IXeff),LFP.sFreq*5,LFP.sFreq*2.5,[5 81],LFP.sFreq);
pw = 10*log10(abs(p));
plot(pw(1,:),pw(2,:),'.')
xlabel('theta pow');ylabel('gam pow');
lsline

% set(gca,'yscale','log')
% set(gca,'YLim',[2 120])
%% power power correlation
CC = corrcoef(cfs(1:20:end,:));
figure
subplot(1,2,1)
imagesc(out_fqs,out_fqs,CC)
axis xy
colorbar
caxis([-.1 .4])
ylabel('power power correlation')

% angle angle correlation
subplot(1,2,2)
CC = corrcoef(A(1:20:end,:));
imagesc(out_fqs,out_fqs,CC)
axis xy
colorbar
caxis([-.1 .2])
ylabel('angle angle correlation')
%%
% See if the single-unit activity correlates more with a particular power
% band...
SU = sc(:,1); tit='PC1';
%  SU = sc(:,2); tit='PC2';
%  SU = sc(:,3); tit='PC3';
% SU = sum(Q,2); tit='Summed';
% SU = sc(randperm(Rows(sc)),1); tit='shuf'; % sanity check

cfs2 = interp1_matrix(LFP.t_uS(IXeff),cfs,Qt);
mod = fitlm(cfs2,SU);
figure
plot(out_fqs,mod.Coefficients.Estimate(2:end))
plot_ref_line(0,'orientation','horiz')
GIX = mod.Coefficients.pValue(2:end)<0.001;
plot(out_fqs(GIX),zeros(1,sum(GIX)),'r*')
xlabel('Hz')
title(tit)
%% do a correlation.
r = [];p = [];
for iC = 1:Cols(cfs2)
    [r(iC),p(iC)] = corr(cfs2(:,iC),SU);
end
figure
plot(out_fqs, r)
hold on
GIX = p < 0.001;
plot(out_fqs(GIX),zeros(1,sum(GIX)),'r*')
xlabel('Hz')
title(tit)
plot_ref_line(0,'orientation','horiz')
%% do each neuron..
B = [];
P = [];
for iN = 1:Cols(Q)
    if 0
        mod = fitlm(cfs2,Q(:,iN));
        B(iN,:) = mod.Coefficients.Estimate(2:end);
        P(iN,:) = mod.Coefficients.pValue(2:end);
    else
        for iC = 1:Cols(cfs2)
            [B(iN,iC),P(iN,iC)] = corr(cfs2(:,iC),Q(:,iN));
        end
    end
end
figure
B(P > 0.001) = nan;
subplot(2,1,1)
imagesc(out_fqs,[],B);
colormap_nan_color(B)
colorbar('NorthOutside')
ylabel('Neuron')
title(tit)
subplot(2,1,2)
plot_confidence_intervals(out_fqs,B)
plot_ref_line(0,'orientation','horiz')


%%
if DO_CFC
    %% Do CFC - jeezus canolty duprelatour tort -
    met = 'duprelatour';
    low_rng = [2:.5:27];
    [CFC] = SPEC_cross_fq_coupling_comod_dupre2017(LFP.LFP(IXeff),LFP.sFreq,low_rng,met);
    
    figure
    imagesc(CFC.low_fq_range,CFC.high_fq_range,OUT.CM)
    title(met)
    xlabel('Low Frequency (Hz)')
    ylabel('High Frequency (Hz)')
    axis xy
    colorbar
    colormap(viridis)
    
    ix = find(IXeff);
    edges = round(linspace(ix(1),ix(end),8));
    mean(diff(edges))/60/LFP.sFreq
    CM = [];
    for ii = 1:length(edges)-1
        [CFC] = SPEC_cross_fq_coupling_comod_dupre2017(LFP.LFP(edges(ii):edges(ii+1)),LFP.sFreq,low_rng,met);
        CM(:,:,ii) = CFC.CM;
        fprintf('.')
    end
    
    figure
    for ii = 1:length(edges)-1
        subplot(length(edges)-1,1,ii)
        imagesc(CFC.low_fq_range,CFC.high_fq_range,CM(:,:,ii))
        axis xy
    end
    figure
    imagesc(CFC.low_fq_range,CFC.high_fq_range,nanmean(CM,3))
    title(met)
    xlabel('Low Frequency (Hz)')
    ylabel('High Frequency (Hz)')
    axis xy
    colorbar
    colormap(viridis)
    plot_horiz_line_at_zero(80);
end