% Directory that has the LFP subfolder...
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all;
fqs = 1:.5:170;
win_sec = 20;
skip = 1;
data_dir = 'C:\Data\NSandB_course2023\23240\1_7_10_23\vSTR_Ket_Neuro2_g0\vSTR_Ket_Neuro2_g0_imec0'; % high dose
 % data_dir =
 % 'C:\Data\NSandB_course2023\23239\01_7_10_23\vSTR_Ket_Neuro2_g0\vSTR_Ket_Neuro2_g0_imec0';
 % % low dos
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = dir(fullfile(data_dir,'LFP2\*.mat'));
% figure out the depth for each channel.
type = 'spec'; % 'cwt'
for iF = 1:skip:length(d)
    load(fullfile(data_dir,"LFP2",d(iF).name))
    if iF == 1
        DEPTHTBL = NPXL_Depth_From_Meta(LFP.original_meta);
        % figure
        % imagesc(zscore([DEPTHTBL.ChannelID DEPTHTBL.SiteID DEPTHTBL.x_uM DEPTHTBL.y_uM DEPTHTBL.Depth_uM] ))
    end
    % Get the channel and determine the depth.
    ix = find(DEPTHTBL.Ch == LFP.Channel);
    Depth_uM = DEPTHTBL.Depth_uM(ix);
    if isempty(Depth_uM)
        disp('woa')
        continue
    end
    switch type
        case 'cwt'
            SPEC = SPEC_cwt_cowen(LFP.Data,LFP.new_sFreq,fqs);
            SPEC = abs(SPEC);
            SPEC2 = log10(SPEC);
            t_min = (0:(Cols(SPEC)-1))/(LFP.new_sFreq*60);
        case 'spec'
            [SPEC,~,t_sec] = spectrogram(LFP.Data,round(LFP.new_sFreq*win_sec),round(LFP.new_sFreq*win_sec/2),fqs,LFP.new_sFreq);
            SPEC = abs(SPEC);
            SPEC2 = log10(SPEC);
            t_min = t_sec/60;
    end
    NIX = t_min < 10;
    mn = mean(SPEC(:,NIX),2);
    sd = std(SPEC(:,NIX),[],2);
    SPECz = (SPEC-mn)./sd;

    SPECz = movmedian(SPECz',10)';
    SPECz = movmean(SPECz',20)';
    SPECz = movmean(SPECz,2);

    ALL(iF).SPECz = SPECz;
    ALL(iF).SPEClog = SPEC2;
    ALL(iF).t_min = t_min;
    ALL(iF).fqs = fqs;
    ALL(iF).Depth_uM = Depth_uM;
 
end

[s,six] = sort([ALL.Depth_uM]);
%%
close all
for ii = 1:length(six)
    sorted_ix = six(ii); % for some reason the sorting is odd- some nunmber are out of order despite this very simple code.
    figure(1)
    clf
    % subplot(2,1,1)
    % imagesc(t_min,fqs,SPEC2)
    % axis xy
    % pubify_figure_axis
    % xlabel('min');ylabel('Hz')
    % % clim([-1.5 3])
    % title(sprintf('%s depth %d uM', d(iF).name, Depth_uM),'Interpreter','none','FontSize',8)
    % 
    % subplot(2,1,2)

    imagesc(t_min,fqs,ALL(sorted_ix).SPECz)
    axis xy
    pubify_figure_axis
    xlabel('min');ylabel('Hz')
    colorbar_label('z scores')
    clim([-1 5])
    set(gcf,'Position',[1          49        1536         837])
    title(sprintf('%s depth %d uM', d(sorted_ix).name, ALL(sorted_ix).Depth_uM+680),'Interpreter','none','FontSize',12)
    pause
    % saveas(gcf,strrep(d(ix).name,'.mat',[ type '.png']))

end

% for ii = 1:length(six)
%      ix = six(ii);
% 
%     ALL(ix).Depth_uM + 1
% end