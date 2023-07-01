% Q10_What_does_LID_and_ket_do_to_specpower_Abhi_poster()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_dir = 'E:\SfN_LFP';
%data_dir = 'G:\SfN_LFP';
spec_fq = [2:.25:150];
sFreq = 500;
baseline_min = [-60 -40];
ana_range_min = [-60 120];
% ana_range_min = [-60 60];
spec_window_sec = 30;
ranges_to_analyze_min = [-35 -5; 5 35; 40 70];
% ranges_to_analyze_min = [-35 -5; 5 35; 40 60];
freqs_to_ana = [15 25; 30 34; 40 55; 63 66; 68 72; 78 90; 98 102; 110 114; 120 160; 168 172];
freq_labels = round(mean(freqs_to_ana,2));
ana_t_sec = ana_range_min(1)*60:1/sFreq:ana_range_min(2)*60;
ana_t_pos_sec = ana_range_min(1)*60:1/50:ana_range_min(2)*60;

% conditions = {'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1' 'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Saline_LID_Unlesion_M1' ...
%     'Saline&Ketamine_Control_M1' 'Saline&Ketamine_Control_Striatum' 'Saline&Ketamine_LID_Lesion_M1' 'Saline&Ketamine_LID_Lesion_Striatum' ''};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The names must correspond to the directories...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALL.SPEC = [];
ALL.BYFRQ = [];
ALL.PSD_base = [];
ALL.PSD_pre = [];
ALL.PSD_post_1 = [];
ALL.PSD_post_2 = [];
% conditions_to_compare = {'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1'};
% conditions_to_compare = {'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Lesion_M1' };
conditions_to_compare = {'LDOPA&Saline_LID_Unlesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1' };
% conditions_to_compare = {'Saline&Ketamine_LID_Unlesion_Striatum' 'Saline&Ketamine_LID_Lesion_Striatum' 'Saline&Ketamine_Control_Striatum'};
% conditions_to_compare = {'Saline&Ketamine_LID_Unlesion_M1' 'Saline&Ketamine_LID_Lesion_M1' 'Saline&Ketamine_Control_M1'};
for iD = 1:length(conditions_to_compare)
    files = dir(fullfile(data_dir,conditions_to_compare{iD},'*.mat'));
    for iF = 1:length(files)
        D = load(fullfile(data_dir,conditions_to_compare{iD},files(iF).name));
        t_sec = (0:(length(D.LFP.data)-1))/D.LFP.LFP_sFreqj;
        if any(contains(fieldnames(D.LFP),'Sal_Start_min'))
            t_sec = t_sec - D.LFP.Sal_Start_min*60;
        else
            t_sec = t_sec - D.LFP.Ket_Start_min*60;
        end
        LFP = single(interp1(t_sec,double(D.LFP.data),ana_t_sec));
        LFP = [ana_t_sec(:) LFP(:)];
        O = get_specs(LFP,sFreq,spec_fq,ranges_to_analyze_min,baseline_min,spec_window_sec);
        F = extract_freq(O.sz,O.spec_fq,freqs_to_ana);
        ALL.SPEC{iD}(:,:,iF) = O.sz;
        ALL.BYFRQ{iD}(:,:,iF) = F;
        ALL.PSD_base{iD}(iF,:) = O.psd_base;
        ALL.PSD_pre{iD}(iF,:) = O.psd(1,:);
        ALL.PSD_post_1{iD}(iF,:) = O.psd(2,:);
        ALL.PSD_post_2{iD}(iF,:) = O.psd(3,:);
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

clrs = lines(5);
figure
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(3,:,:) - (ALL.BYFRQ{ii}(2,:,:)+ALL.BYFRQ{ii}(4,:,:))/2)';
    M = squeeze(ALL.BYFRQ{ii}(3,:,:))';
    Ms = movmedian(M,10);
    plot_confidence_intervals(O.t_sec/60,Ms,[],clrs(ii,:))
end
pubify_figure_axis
plot_ref_line(0);
title(conditions_to_compare)
ylabel('z 50Hz')
xlabel('min')

figure
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(9,:,:) - (ALL.BYFRQ{ii}(8,:,:)+ALL.BYFRQ{ii}(10,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(9,:,:) - (ALL.BYFRQ{ii}(7,:,:)))';
    M = squeeze(ALL.BYFRQ{ii}(9,:,:))';
    Ms = movmedian(M,10);
    plot_confidence_intervals(O.t_sec/60,Ms,[],clrs(ii,:))
end
pubify_figure_axis
plot_ref_line(0);
title(conditions_to_compare)
ylabel('z 140Hz')
xlabel('min')

figure
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(3,:,:) - ALL.BYFRQ{ii}(4,:,:))';
    M = squeeze(ALL.BYFRQ{ii}(6,:,:) - (ALL.BYFRQ{ii}(5,:,:)+ALL.BYFRQ{ii}(7,:,:))/2)';
    Ms = movmedian(M,5);
    plot_confidence_intervals(O.t_sec/60,Ms,[],clrs(ii,:))
end
pubify_figure_axis
plot_ref_line(0);
title(conditions_to_compare)
ylabel('Focal 80Hz')
xlabel('min')

%%
function F = extract_freq(s,frqs,bands)
F = [];
for iF = 1:Rows(bands)
    IX = frqs >= bands(iF,1) & frqs <= bands(iF,2);
    F(iF,:) = nanmean(s(IX,:));
end
end

function O = get_specs(LFP,sFreq,spec_fq,ranges_to_analyze_min,baseline_min,spec_window_sec)
[s,~,t_sec] = spectrogram(LFP(:,2),sFreq*spec_window_sec,sFreq*spec_window_sec/2,spec_fq,sFreq);
t_sec = t_sec + LFP(1,1) + spec_window_sec/2;
s = 10*log10(abs(s));
IX = t_sec > baseline_min(1)*60 & t_sec < baseline_min(2)*60;
mn = trimmean(s(:,IX),2,2);
sd = trimstd(s(:,IX),2,2);
sz = (s-mn)./sd;
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