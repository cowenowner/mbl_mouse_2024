function Q11_What_does_LID_and_ket_do_to_CFC_Abhi_poster(conditions_to_compare)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars
if 0
    % Do it all. Takes a while.
    conditions_to_compare{1} = {'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1'};
    conditions_to_compare{2} = {'LDOPA&Saline_LID_Lesion_M1'  'Saline&Ketamine_LID_Unlesion_Striatum' };
    conditions_to_compare{3} = {'Saline&Ketamine_Control_M1' 'Saline&Ketamine_LID_Lesion_M1' };
    conditions_to_compare{4} = {'Saline&Ketamine_Control_Striatum' 'Saline&Ketamine_LID_Lesion_Striatum' };
    for ii = 1:length(conditions_to_compare)
        Q11_What_does_LID_and_ket_do_to_CFC_Abhi_poster(conditions_to_compare{ii});
        ca
        ii
    end
    
end
if nargin == 0
    % conditions_to_compare = {'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1'};
    % conditions_to_compare = {'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Lesion_M1' };
    conditions_to_compare = {'Saline&Ketamine_Control_M1' 'Saline&Ketamine_LID_Lesion_M1' };
    % conditions_to_compare = {'Saline&Ketamine_Control_Striatum' 'Saline&Ketamine_LID_Lesion_Striatum' };
    % conditions_to_compare = {'Saline&Ketamine_LID_Lesion_Striatum' 'Saline&Ketamine_LID_Unlesion_Striatum' };
end
clrs = lines(15);

data_dir = 'G:\SfN_LFP';
cfc_method = 'duprelatour'; % this seems to work well in practice. it is a version of DAR

% cfc_method = 'tort'; % seems super slow.
% cfc_method = 'canolty'; % seems super slow.
spec_fq = 1:.25:190;
% low_fqs_hz = 1:.5:25;
low_fq_range_hz = [1 20];
step = 1;
low_fqs_hz = low_fq_range_hz(1):step:low_fq_range_hz(end);

spec_window_sec = 5;
sFreq = 500;
ana_range_min = [-80 120];
block_size_sec = 4; % for cfc
ranges_to_analyze_min = [-70 -55; -25 -5; 5 25; 40 60; 80 100];
ana_t_sec = ana_range_min(1)*60:1/sFreq:ana_range_min(2)*60;

% conditions = {'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1' 'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Saline_LID_Unlesion_M1' ...
%     'Saline&Ketamine_Control_M1' 'Saline&Ketamine_Control_Striatum' 'Saline&Ketamine_LID_Lesion_M1' 'Saline&Ketamine_LID_Lesion_Striatum' ''};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The names must correspond to the directories...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALL.SPEC = [];


% conditions_to_compare = {'LDOPA&Ketamine_LID_Unlesion_M1' 'LDOPA&Ketamine_LID_Lesion_M1' };
% conditions_to_compare = {'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Lesion_M1' };
% conditions_to_compare = {'Saline&Ketamine_Control_Striatum' 'Saline&Ketamine_LID_Lesion_Striatum' };

for iD = 1:length(conditions_to_compare)
    files = dir(fullfile(data_dir,conditions_to_compare{iD},'*.mat'));
    for iF = 3:length(files)
        D = load(fullfile(data_dir,conditions_to_compare{iD},files(iF).name));
        t_sec = (0:(length(D.LFP.data)-1))/D.LFP.LFP_sFreqj;
        if any(contains(fieldnames(D.LFP),'Sal_Start_min'))
            t_sec = t_sec - D.LFP.Sal_Start_min*60;
        else
            t_sec = t_sec - D.LFP.Ket_Start_min*60;
        end
        LFP(:,1) = ana_t_sec(:);
        LFP(:,2) = interp1(t_sec,double(D.LFP.data)*D.LFP.to_uV_conversion,ana_t_sec);
        [bad_times_s,good_times_s,BIX] = Clean_LFP(LFP,sFreq,'remove_high_voltage_spindles',true,'artifact_thresh',2000);
        %         [hvs_times_s, PAR, no_hvs_times_s,HVS_IX,NOHVS_IX] = High_voltage_spindle_detector(LFP,D.LFP.LFP_sFreqj);
        % add a little more of a buffer
        
        [LFP_nohvs, GIX] = Restrict(LFP,good_times_s);
        LFP_zeroed = LFP;
        LFP_zeroed(~GIX,2) = 0;
        LFP_zeroed(isnan(LFP(:,2)),2:end) = 0;
        
        SPEC = get_specs(LFP_zeroed,D.LFP.LFP_sFreqj,spec_fq,ranges_to_analyze_min,spec_window_sec);
        plot_specs(SPEC,[conditions_to_compare{iD} ' ' files(iF).name])
        O = do_cfc(LFP_nohvs,D.LFP.LFP_sFreqj,low_fqs_hz,cfc_method,ranges_to_analyze_min,block_size_sec);
        plot_cfc(O,[conditions_to_compare{iD} ' ' files(iF).name])
        
        %         F = extract_freq_2D(O.sz,O.spec_fq,freqs_to_ana);
        ALL.SPEC{iD}(:,:,iF) = SPEC.sz;
        ALL.PSD{iD}(:,:,iF) = SPEC.psd;
        ALL.CFC{iD}(:,:,:,iF) = O.CM;
        fprintf('.')
    end
end
%%
figure
a = []; cnt = 1;
for iR = 1:length(conditions_to_compare)
    %     for iC = 1:Rows(ranges_to_analyze_min)
    for iC = 1:4
        %         if iC > 1
        
        a(cnt) = subplot_ij(length(conditions_to_compare), Rows(ranges_to_analyze_min), iR , iC);
        cnt = cnt + 1;
        
        imagesc(O.CFC{1}.low_fq_range,O.CFC{1}.high_fq_range,nanmean(ALL.CFC{iR}(:,:,iC,:),4))
        %          imagesc(O.CFC{1}.low_fq_range,O.CFC{1}.high_fq_range,nanmean(ALL.CFC{iR}(:,:,iC,:)-ALL.CFC{iR}(:,:,1,:),4))
        % end
        %          imagesc(O.CFC{1}.low_fq_range,O.CFC{1}.high_fq_range,nanmean(ALL.CFC{iR}(:,:,iC,:),4)./nanstd(ALL.CFC{iR}(:,:,iC,:),[],4))
        %         imagesc(O.CFC{1}.low_fq_range,O.CFC{1}.high_fq_range,nanstd(ALL.CFC{iR}(:,:,iC,:),[],4))
        title({conditions_to_compare{iR} sprintf('%d min',ranges_to_analyze_min(iC,:))})
        xlabel('Low Freq (Hz)')
        ylabel('High Freq (Hz)')
        axis xy
        colorbar
        %         colormap(viridis)
    end
end
% sgtitle(conditions_to_compare)
% equalize_color_axes(a)
%%
clear LFP D t_sec GIX LFP_nohvs LFP_zeroed ana_t_sec
fname = 'Q11_';
for ii = 1:length(conditions_to_compare)
    fname = [fname conditions_to_compare{ii}];
end
fname = [fname '.mat'];
save(fullfile('C:\Temp',fname));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function F = extract_freq_2D(CFC,frqs_x,frqs_y,bands_x,bands_y)
F = [];
for iX = 1:Rows(bands_x)
    for iY = 1:Rows(bands_y)
        IXx = frqs_x >= bands_x(iF,1) & frqs_x <= bands_x(iF,2);
        IXy = frqs_y >= bands_y(iF,1) & frqs_y <= bands_y(iF,2);
        F(iY,iX) = nanmean(CFC(IXx,IXy));
    end
end
end

function O = do_cfc(LFP,sFreq,low_fqs,method,ranges_to_analyze_min,block_size_sec)
%% Do CFC
%     method = 'duprelatour'; % jeezus canolty duprelatour tort -
CM = [];
CFC = [];
for iR = 1:Rows(ranges_to_analyze_min)
    IX = LFP(:,1) > ranges_to_analyze_min(iR,1)*60 & LFP(:,1) < ranges_to_analyze_min(iR,2)*60;
    CFC{iR} = SPEC_cross_fq_coupling_comod_dupre2017(LFP(IX,2),sFreq,low_fqs,method,sFreq*block_size_sec*60);
    %     tmp = SPEC_cross_fq_coupling_comod_dupre2017(LFP(IX,2),sFreq,low_fqs,method,[],100);
    CM(:,:,iR) = CFC{iR}.CM;
end
O.CFC = CFC;
O.CM = CM;
O.ranges_to_analyze_min = ranges_to_analyze_min;
O.method = method;
% plot_cfc(O,'')
end

function plot_cfc(O,tit)
figure
aa = [];
for ii = 1:length(O.CFC)
    aa(ii) = subplot(2,ceil(length(O.CFC)/2),ii);
    imagesc(O.CFC{ii}.low_fq_range,O.CFC{ii}.high_fq_range,O.CFC{ii}.CM)
    title(num2str(ii))
    xlabel('Low Frequency (Hz)')
    ylabel('High Frequency (Hz)')
    axis xy
    colorbar
    %     set(gca,'YLim',[O.CFC{ii}.high_fq_range(1) 180])
    %     coSlormap(viridis)
end
sgtitle(tit)
% equalize_color_axes(aa)

end


function O = get_specs(LFP,sFreq,spec_fq,ranges_to_analyze_min,spec_window_sec)
[s,~,t_sec] = spectrogram(LFP(:,2),sFreq*spec_window_sec,sFreq*spec_window_sec/2,spec_fq,sFreq);
t_sec = t_sec + LFP(1,1) + spec_window_sec/2;
s = 10*log10(abs(s));
IX = t_sec > ranges_to_analyze_min(1)*60 & t_sec < ranges_to_analyze_min(2)*60;
% mn = trimmean(s(:,IX),2,'round' ,2);
mn = trimmean_cowen(s(:,IX)',2)';
sd = trimstd(s(:,IX),2,2);
sz = (s-mn)./sd;
O.s = s;
O.sz = sz;
O.t_sec = t_sec;
O.spec_fq = spec_fq;
O.ranges_to_analyze_min = ranges_to_analyze_min;
O.spec_window_sec = spec_window_sec;

IX = LFP(:,1) >  ranges_to_analyze_min(1)*60 & LFP(:,1) < ranges_to_analyze_min(2)*60;
p = pwelch(LFP(IX,2),sFreq*spec_window_sec,sFreq*spec_window_sec/2,spec_fq,sFreq);
p = 10*log10(abs(p));
O.psd_base = p;

for iR = 1:Rows(ranges_to_analyze_min)
    IX = LFP(:,1) > ranges_to_analyze_min(iR,1)*60 & LFP(:,1) < ranges_to_analyze_min(iR,2)*60;
    p = pwelch(LFP(IX,2),sFreq*spec_window_sec,sFreq*spec_window_sec/2,spec_fq,sFreq);
    p = 10*log10(abs(p));
    O.psd(iR,:) = p;
end

% Power correlation
R = [];
for ii = 1:Rows(ranges_to_analyze_min)
    IX = O.t_sec/60 > ranges_to_analyze_min(ii,1) & O.t_sec/60 < ranges_to_analyze_min(ii,2);
    V = O.sz(:,IX)';
    V = V(~(isnan(sum(V,2)) | isinf(sum(V,2))) ,:);
    R = [R corr(V,'Type','Pearson')];
    %     R = [R corr(V,'Rows','pairwise','Type','Pearson')];
end
O.R = R;

end

function plot_specs(O,tit)
figure
subplot(3,2,1:2)
imagesc(O.t_sec/60,O.spec_fq,O.s)
colorbar
axis xy
xlabel('min')
ylabel('Hz')
c = caxis;
caxis([-20 c(end)])
plot_ref_line(0)
plot_markers_simple(O.ranges_to_analyze_min)
title(tit)

subplot(3,2,3:4)
imagesc(O.t_sec/60,O.spec_fq,O.sz)
colorbar
axis xy
plot_markers_simple(O.ranges_to_analyze_min)
plot_ref_line(0)
title('Z scored spec')

subplot(3,2,5)
plot(O.spec_fq,O.psd,'LineWidth',2)
legend('1','2','3','4','5'); legend boxoff
xlabel('Hz')
pubify_figure_axis

subplot(3,2,6)
imagesc([],O.spec_fq,O.R)
% caxis(prctile(R(:),[2 80]))
caxis([-.1 .75])
axis xy; colorbar
xlabel('epoch')
ylabel('Hz')
title('corr of power in spectrogram')
end



