% Q10_What_does_LID_and_ket_do_to_specpower_Abhi_poster()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variab-les.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_dir = 'E:\SfN_2022';
%data_dir = 'G:\SfN_LFP';
spec_fq = [2:.25:200];
% spec_fq = [1:.5:250]; % used to generate psd for spike field plot
sFreq = 500;
baseline_min = [-85 -65]; % LDO&Sal and LDO&Ket sessions
% baseline_min = [-50 -35]; % used only to genrate spike field psd plots; Sal&KEt sessions
% baseline_min = [-58 -32];
% spindle_min = [-58 -2];
% ana_range_min = [-58 80];
ana_range_min = [-88 100]; % LDO&Sal/Ket sessions
% ana_range_min = [-50 100]; % Sal and Ket sessions
spec_window_sec = 10;
% ranges_to_analyze_min = [-20 -5; 2 30; 60 80]; % used to generate for spike field plots only
% ranges_to_analyze_min = [-55 -35; -25 -5; 5 25; 130 150];
ranges_to_analyze_min = [-55 -35; -25 -5; 2 30; 60 80]; % LDO&Sal/Ket sessions
% ranges_to_analyze_min = [-20 -5; 2 30; 60 80]; Sal&Ket sessions
% freqs_to_ana = [15 25; 30 34; 40 55; 40 75; 63 66; 68 72; 78 90; 98 102; 110 114; 120 170; 178 182];
% freqs_to_ana = [5 11; 15 25; 35 39; 45 60; 40 75; 66 70; 65 69; 75 92; 98 102; 110 114; 120 170; 178 182]; % Use this
freqs_to_ana = [1 3; 5 11; 15 35; 35 75; 45 60; 75 95];
freq_labels = round(mean(freqs_to_ana,2));
ana_t_sec = ana_range_min(1)*60:1/sFreq:ana_range_min(2)*60;
ana_t_base = baseline_min(1)*60:1/sFreq:baseline_min(2)*60;
% ana_t_pos_sec = ana_range_min(1)*60:1/50:ana_range_min(2)*60;

% conditions = {'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1' 'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Saline_LID_Unlesion_M1' ...
%     'Saline&Ketamine_Control_M1' 'Saline&Ketamine_Control_Striatum' 'Saline&Ketamine_LID_Lesion_M1' 'Saline&Ketamine_LID_Lesion_Striatum' ''};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The names must correspond to the directories...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALL.SPEC = [];
ALL.SPEC_nobase = [];
ALL.BYFRQ = [];
ALL.PSD_base = [];
ALL.PSD_pre = [];
ALL.PSD_post_1 = [];
ALL.PSD_post_2 = [];
ALL.PSD_post_3 = [];
% conditions_to_compare = {'LDOPA&Ketamine_LID_UnLesion_M1' 'LDOPA&Ketamine_LID_Lesion_M1'};
% conditions_to_compare = {'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Lesion_M1' };
% conditions_to_compare = {'LDOPA&Saline_LID_Unlesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1' };
% conditions_to_compare = {'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1' };
% conditions_to_compare = {'Saline&Ketamine_LID_Unlesion_Striatum' 'Saline&Ketamine_LID_Lesion_Striatum' 'Saline&Ketamine_Control_Striatum'};
% conditions_to_compare = {'Saline&Ketamine_LID_Unlesion_M1' 'Saline&Ketamine_LID_Lesion_M1' 'Saline&Ketamine_Control_M1'};
% conditions_to_compare = {'LDOPA&Saline_Rat320_Lesion_M1' 'LDOPA&Saline_Rat337_Lesion_M1' 'LDOPA&Saline_Rat342_Lesion_M1' 'LDOPA&Saline_Rat378_Lesion_M1' 'LDOPA&Saline_Rat380_Lesion_M1'};
conditions_to_compare = {'LDOPA&Ketamine_Rat320_Lesion_M1' 'LDOPA&Ketamine_Rat342_Lesion_M1'};
% conditions_to_compare = {'LDOPA&Ketamine_Rat320_UnLesion_M1' 'LDOPA&Ketamine_Rat330_UnLesion_M1' 'LDOPA&Ketamine_Rat342_UnLesion_M1'};
% conditions_to_compare = {'LDOPA&Ketamine_Rat350_Lesion_M1' 'LDOPA&Ketamine_Rat352_Lesion_M1' 'LDOPA&Ketamine_Rat350_Unlesion_M1' 'LDOPA&Ketamine_Rat352_Unlesion_M1'};
% conditions_to_compare = {'Saline&Ketamine_Rat350_Lesion_M1' 'Saline&Ketamine_Rat352_Lesion_M1' 'Saline&Ketamine_Rat350_Unlesion_M1' 'Saline&Ketamine_Rat352_Unlesion_M1'};
% conditions_to_compare = {'LDOPA&Saline_Rat350_Lesion_M1' 'LDOPA&Saline_Rat352_Lesion_M1'};
for iD = 1:length(conditions_to_compare)
    files = dir(fullfile(data_dir,conditions_to_compare{iD},'*.mat'));
    for iF = 1:length(files)
        D = load(fullfile(data_dir,conditions_to_compare{iD},files(iF).name));
        t_sec = (0:(length(D.LFP.data)-1))/D.LFP.LFP_sFreqj;
        t_uS = 1e6*t_sec(:);
        if any(contains(fieldnames(D.LFP),'Sal_Start_min'))
            t_sec = t_sec - D.LFP.Sal_Start_min*60;
            t_uS = t_uS - D.LFP.Sal_Start_min*60e6;
            inj_time = D.LFP.Sal_Start_min;
        else
            t_sec = t_sec - D.LFP.Ket_Start_min*60;
            t_uS = t_uS - D.LFP.Ket_Start_min*60e6;
            inj_time = D.LFP.Ket_Start_min;
        end
        
        %Restrict the LFPdata to the analysis range so that the baseline is
        %the first point
        IX_ana = t_uS > ana_range_min(1)*60e6 & t_uS < ana_range_min(2)*60e6;
        IX_base = t_uS > baseline_min(1)*60e6 & t_uS < baseline_min(2)*60e6;
        IX_sec_base = t_sec > baseline_min(1)*60 & t_sec < baseline_min(2)*60;
%         t_sec_base = t_sec(IX_sec_base);
        base_time = t_uS(IX_base);
        LFP_ana = D.LFP.data(IX_ana);
        t_uS = t_uS(IX_ana);
        t_sec = t_sec(IX_ana);

        
        % remove spindles from baseline period
        if exist('spindle_min','var')
            IX_spindlet = t_uS > spindle_min(1)*60e6 & t_uS < spindle_min(2)*60e6;
            ForHVS_time = t_uS(IX_spindlet);
            res_data = LK_HVS_detector([ForHVS_time,double(D.LFP.data(IX_spindlet))],D.LFP.LFP_sFreqj);
        else
%             base_data = [];
            res_data = LK_HVS_detector([base_time,double(D.LFP.data(IX_base))],D.LFP.LFP_sFreqj);
        end
        sp_indices = find(res_data.sp_indicies == 1);
%         sp_indices_t = t_uS(res_data.sp_indicies);
%         sp_indices_sec = sp_indices_t./1000000;
        LFP = double(LFP_ana);
%         for ii = 1:length(sp_indices)
%             LFP(sp_indices(ii)) = 0;
%         end
        sp_times = [];
%         sp_times_all = t_uS(sp_indices(:,1));
        sp_times(:,:) = res_data.sp_times + t_uS(1);
        LFP = single(interp1(t_sec,LFP,ana_t_sec));
        if any(isnan(LFP))
            n = find(isnan(LFP));
            for ii = 1:length(n)
                LFP(n(ii)) = 0;
            end
        else
        end
%         sp_sec_ind = single(interp1(t_sec_base,res_data.sp_indicies,ana_t_base));
        LFP = [ana_t_sec(:) LFP(:)];
        
        O = get_specs(LFP,sFreq,spec_fq,ranges_to_analyze_min,baseline_min,spec_window_sec,sp_times);
        F = extract_freq(O.sz,O.spec_fq,freqs_to_ana);
        ALL.SPEC_nobase{iD}(:,:,iF) = O.s; 
        ALL.SPEC{iD}(:,:,iF) = O.sz;
        ALL.BYFRQ{iD}(:,:,iF) = F;
        ALL.PSD_base{iD}(iF,:) = O.psd_base;
        ALL.PSD_pre{iD}(iF,:) = O.psd(1,:);
        ALL.PSD_post_1{iD}(iF,:) = O.psd(2,:);
        ALL.PSD_post_2{iD}(iF,:) = O.psd(3,:);
        ALL.PSD_post_3{iD}(iF,:) = O.psd(4,:);
    end
end
%%
figure
a = [];
for ii = 1:length(conditions_to_compare)
    a(ii) = subplot(length(conditions_to_compare),1,ii);
    imagesc(O.t_sec/60,O.spec_fq,nanmean(ALL.SPEC{ii},3))
    axis xy
    xlabel('min')
    ylabel('Hz')
    title(conditions_to_compare{ii})
    colorbar
    plot_ref_line(0)
    pubify_figure_axis
end
equalize_color_axes(a)

pubify_figure_axis
ax=gca;
ax.FontSize=28;
set(gca,'ylim',[-1.25 3.25])
% set(gca,'xlim',[-60 120])
plot_ref_line(0,'line_width',2,'style',':')
plot_ref_line(-60,'line_width',2,'style',':')

clrs = lines(5);
figure
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(3,:,:) - (ALL.BYFRQ{ii}(2,:,:)+ALL.BYFRQ{ii}(5,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(4,:,:)- (ALL.BYFRQ{ii}(6,:,:)+ALL.BYFRQ{ii}(8,:,:))/2)';
    M = squeeze(ALL.BYFRQ{ii}(1,:,:))';
    Ms = movmedian(M,10,2);
    plot_confidence_intervals(O.t_sec/60,Ms,[],clrs(ii,:))
%     plot(O.t_sec/60,Ms)
%     Ms_m = Ms-mean(Ms);
%     figure
%     plot(O.t_sec/60,Ms_m)
%     figure
%     scatter(Ms_m(1,:),Ms_m(2,:))
%     SEM = std(Ms)./sqrt(length(Ms_m));
%     errorbar(O.t_sec/60,Ms_m,SEM)
end
pubify_figure_axis
% set(gca,'ylim',[-1 6.3])
plot_ref_line(0);
plot_ref_line(-60);
ylabel('1-3Hz zscore power')
xlabel('Time (min)')
title(conditions_to_compare)

clrs = lines(5);
figure
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(3,:,:) - (ALL.BYFRQ{ii}(2,:,:)+ALL.BYFRQ{ii}(5,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(4,:,:)- (ALL.BYFRQ{ii}(6,:,:)+ALL.BYFRQ{ii}(8,:,:))/2)';
    M = squeeze(ALL.BYFRQ{ii}(3,:,:))';
    Ms = movmedian(M,10,2);
    plot_confidence_intervals(O.t_sec/60,Ms,[],clrs(ii,:))
%     plot(O.t_sec/60,Ms)
%     Ms_m = Ms-mean(Ms);
%     figure
%     plot(O.t_sec/60,Ms_m)
%     figure
%     scatter(Ms_m(1,:),Ms_m(2,:))
%     SEM = std(Ms)./sqrt(length(Ms_m));
%     errorbar(O.t_sec/60,Ms_m,SEM)
end
pubify_figure_axis
% set(gca,'ylim',[-1 6.3])
plot_ref_line(0);
plot_ref_line(-60);
ylabel('15-30 Hz zscore power')
xlabel('Time (min)')
title(conditions_to_compare)


figure
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(3,:,:) - ALL.BYFRQ{ii}(4,:,:))';
    M = squeeze(ALL.BYFRQ{ii}(8,:,:) - (ALL.BYFRQ{ii}(7,:,:)+ALL.BYFRQ{ii}(9,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(7,:,:))';
    Ms = movmedian(M,10,2);
    plot_confidence_intervals(O.t_sec/60,Ms,[],clrs(ii,:))
end
pubify_figure_axis
set(gca,'ylim',[-1.25 3.25])
plot_ref_line(0);
plot_ref_line(-60);
title(conditions_to_compare)
ylabel('Focal 80Hz zscore power')
xlabel('Time (min)')

% Mean of 80Hz power across all rats
pow80=[];
pow80_sd =[];
peak_80 = find(O.t_sec/60 > 5 & O.t_sec/60 < 25);
% peak_inj2 = find(O.t_sec/60 > 5 & O.t_sec/60 < 25);
for ii = 1:length(conditions_to_compare)
    M = squeeze(ALL.BYFRQ{ii}(8,:,:) - (ALL.BYFRQ{ii}(7,:,:)+ALL.BYFRQ{ii}(9,:,:))/2)';
    mean_M = mean(M(:,peak_80));
    sd_M = std(M(:,peak_80));
    pow80 = cat(1,pow80,mean_M);
    pow80_sd = cat(1,pow80_sd,sd_M);
end


m_pow80 = mean(mean(pow80));
sd_pow80 = mean(mean(pow80_sd));


figure
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(9,:,:) - (ALL.BYFRQ{ii}(8,:,:)+ALL.BYFRQ{ii}(10,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(9,:,:) - (ALL.BYFRQ{ii}(7,:,:)))';
    M = squeeze(ALL.BYFRQ{ii}(10,:,:))';
    Ms = movmedian(M,10,2);
    plot_confidence_intervals(O.t_sec/60,Ms,[],clrs(ii,:))
end
pubify_figure_axis
% set(gca,'ylim',[-.7 1.5])
plot_ref_line(0);
plot_ref_line(-30);
title(conditions_to_compare)
ylabel('z 140Hz')
xlabel('min')

% plotting the psd
mn_psd_ket = mean(ALL.PSD_post_2{1});
amn_psd_sal = mean(ALL.PSD_post_2{2});
figure;plot(spec_fq,mn_psd_sal)
hold on 
plot(spec_fq,mn_psd_ket)
pubify_figure_axis
ax=gca;
ax.FontSize=35;
set(gca,'xlim',[0 150])
set(gca,'ylim',[0 50])

% For spike field cohenrance plots
mn_psd=[];
mn_psd(1,:) = mean(ALL.PSD_pre{1});
mn_psd(2,:) = mean(ALL.PSD_post_1{1});
mn_psd(3,:) = mean(ALL.PSD_post_2{1});
mn_psd(4,:) = mean(ALL.PSD_post_{1});

% mn_psd(4,:) = mean(ALL.PSD_post_3{1});
for ii = 1:length(mn_psd(:,1))
   figure;plot(spec_fq,mn_psd(ii,:)) 
   pubify_figure_axis
   ax=gca;
   ax.FontSize=35;
   set(gca,'ylim',[0 40])
   set(gca,'xlim',[0 250])
end
% ploting the psd mean of saline and ketamine sessions for control and LID
mn_psd_unlesion = mean(ALL.PSD_post_2{1});
mn_psd_lesion = mean(ALL.PSD_post_2{2});
mn_psd_control = mean(ALL.PSD_post_2{3});
mn_psd_rat = mean(ALL.PSD_post_2{4});
figure;plot(spec_fq,mn_psd_unlesion)
hold on
plot(spec_fq,mn_psd_lesion)
plot(spec_fq,mn_psd_control)
plot(spec_fq,mn_psd_rat)
set(gca,'ylim',[5 60])
title('Lesioned Sal&Ket 2 to 30 min post ketamine n=4')
set(gca,'ylim',[5 60])
title('Unlesioned Sal&Ket 2 to 30 min post ketamine n=4')
set(gca,'ylim',[5 60])
title('Unlesion,lesion,Control Sal&Ket 2 to 30 min post ketamine n=7,7,10')
%%
function F = extract_freq(s,frqs,bands)
F = [];
for iF = 1:Rows(bands)
    IX = frqs >= bands(iF,1) & frqs <= bands(iF,2);
    F(iF,:) = nanmean(s(IX,:));
end
end

function O = get_specs(LFP,sFreq,spec_fq,ranges_to_analyze_min,baseline_min,spec_window_sec,sp_times)
[s,~,t_sec] = spectrogram(LFP(:,2),sFreq*spec_window_sec,sFreq*spec_window_sec/2,spec_fq,sFreq);
% [sb,~,t_sec_b] = spectrogram(base_data(:,2),sFreq*(spec_window_sec),sFreq*(spec_window_sec)/2,spec_fq,sFreq);
t_sec = t_sec + LFP(1,1) + spec_window_sec/2;
t_uS = linspace(t_sec(1)*1e6,(t_sec(end)*1e6)-(spec_window_sec/2),Cols(t_sec));
s = 10*log10(abs(s));
IX_spec = t_sec > baseline_min(1)*60 & t_sec < baseline_min(2)*60;
% sp_indices_sec = find(s(1,:) == -inf); 
t_uS_sp = [];
for ii = 1:size(sp_times)
    idx_sp = find(t_uS(1,:) >= sp_times(ii,1) & t_uS(1,:) <= sp_times(ii,2));
    if min(idx_sp)==1
        idx_sp_adj = [idx_sp max(idx_sp)+1];
    else
        idx_sp_adj = [idx_sp min(idx_sp)-1 max(idx_sp)+1];
    end
    t_uS_sp = cat(2,t_uS_sp,idx_sp_adj);
end
% turn spindle times into nans for computing mean
for ii = 1:length(t_uS_sp)
    s(:,t_uS_sp(ii)) = nan;
end
sb = s(:,IX_spec);
mn = trimmean(sb,2,2);
% sd = trimstd(sb,2,2);
sz = (s-mn);%./sd;
% turn_to_0 = find(isnan(sz(1,:)));
for ii = 1:length(t_uS_sp)
    sz(:,t_uS_sp(ii)) = 0;
end
O.s = s;
O.sz = sz;
O.t_sec = t_sec;
O.spec_fq = spec_fq;


IX = LFP(:,1) >  baseline_min(1)*60 & LFP(:,1) < baseline_min(2)*60;
p = pwelch(LFP(IX,2),sFreq*spec_window_sec,sFreq*spec_window_sec/2,spec_fq,sFreq);
p = 10*log10(abs(p));
O.psd_base = p;

for iR = 1:Rows(ranges_to_analyze_min)
    IX = LFP(:,1) > ranges_to_analyze_min(iR,1)*60 & LFP(:,1) < ranges_to_analyze_min(iR,2)*60;
    p = pwelch(LFP(IX,2),sFreq*spec_window_sec,sFreq*spec_window_sec/2,spec_fq,sFreq);
    p = 10*log10(abs(p));
    O.psd(iR,:) = p;
end

if 0
    figure
    subplot(2,1,1)
    imagesc(O.t_sec/60,O.spec_fq,O.s)
    colorbar
    axis xy
    xlabel('min')
    ylabel('Hz')
    c = caxis;
    caxis([-20 c(end)])
    
    subplot(2,1,2)
    imagesc(O.t_sec/60,O.spec_fq,O.sz)
    colorbar
    axis xy
    
    figure
    plot(spec_fq,O.psd)
end


% [cfs_tmp,frq_tmp] = cwt(LFP(:,2),'bump',sFreq,'FrequencyLimits',spec_fq([1 end]),'VoicesPerOctave',24);
% [cfs,~,A] = cwt_fix(cfs_tmp,frq_tmp,spec_fq);


end