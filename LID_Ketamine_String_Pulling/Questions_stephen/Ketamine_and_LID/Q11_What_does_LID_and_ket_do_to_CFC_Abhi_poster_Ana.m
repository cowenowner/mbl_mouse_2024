%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q11_What_does_LID_and_ket_do_to_CFC_Abhi_poster_Ana
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
clrs = lines(15);
conditions_to_load = {'LDOPA&Saline_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Lesion_M1'  'LDOPA&Ketamine_LID_Lesion_M1' 'LDOPA&Ketamine_LID_Unlesion_M1' 'Saline&Ketamine_Control_M1' 'Saline&Ketamine_LID_Lesion_M1' 'Saline&Ketamine_Control_Striatum' 'Saline&Ketamine_LID_Lesion_Striatum' 'Saline&Ketamine_LID_Lesion_Striatum' 'Saline&Ketamine_LID_Unlesion_Striatum'};
ddir = 'C:\Temp';
d = dir(fullfile(ddir,'Q11_*.mat'));
intervals_to_plot = 1:5;
cnt = 1;
ALL = [];
for iF = 1:length(d)
    fname = fullfile(ddir,d(iF).name);
    D = load(fname);
    % load and separate all of the data.
    for ii = 1:length(D.ALL.CFC)
        ALL.CFC{cnt} = D.ALL.CFC{ii};
        ALL.PSD{cnt} = D.ALL.PSD{ii};
        ALL.SPEC{cnt} = single(D.ALL.SPEC{ii});
        ALL.cond{cnt} = D.conditions_to_compare{ii};
        cnt = cnt + 1;
    end
end
%% %%%
figure
a = [];
for iF = 1:length(ALL.CFC)
    subplot(length(ALL.CFC),1,iF)
    ALL.SPEC{iF}(isinf(ALL.SPEC{iF})) = nan;
    ALL.SPEC{iF}(ALL.SPEC{iF}==0) = nan;
    M = nanmedian(ALL.SPEC{iF},3);
%     M = convn(M',hanning(40)/sum(hanning(40)),'same')';
%      M = movmean(M,20,'omitnan');
     M = movmedian(M,20,'omitnan');
    imagesc(D.SPEC.t_sec/60,D.SPEC.spec_fq,M)
    axis xy
    ylabel('Hz')
    caxis(prctile(M(:),[3,97]));
%     colormap_nan_color(M)
    title(ALL.cond{iF})
    plot_markers_simple(D.ranges_to_analyze_min)
    colorbar
end
%%
figure
a = [];
cnt = 1;
for iF = 1:length(ALL.CFC)
    for iC = 1:length(intervals_to_plot)
        intv = intervals_to_plot(iC);
        a(cnt) = subplot_ij(length(ALL.CFC), length(intervals_to_plot), iF , iC);
        
        imagesc(D.O.CFC{1}.low_fq_range,D.O.CFC{1}.high_fq_range,nanmean(ALL.CFC{iF}(:,:,intv,:),4))
        caxis_all(cnt,:) = caxis;
        cnt = cnt + 1;

        if iC == 1
            title([ALL.cond{iF} sprintf('%d min',D.ranges_to_analyze_min(intv,:))],'FontSize',6)
        else
            title({sprintf('%d min',D.ranges_to_analyze_min(intv,:))},'FontSize',8)
        end
        if iF == length(ALL.CFC)
            xlabel('Low Freq (Hz)')
        end
        if iC == 1
            ylabel('High Freq (Hz)')
        end
        
        axis xy
        set(gca,'YLim',[D.O.CFC{1}.high_fq_range(1) 200])
        colorbar
        %         colormap(viridis)
    end
end
% sgtitle('CFC')
equalize_color_axes(a,[prctile(caxis_all(:,1),25) prctile(caxis_all(:,2),75)] )
%%
figure
a = [];
cnt = 1;
for iF = 1:length(ALL.PSD)
    for iC = 1:length(intervals_to_plot)
        intv = intervals_to_plot(iC);
        a(cnt) = subplot_ij(length(ALL.PSD), length(intervals_to_plot), iF , iC);
        
        plot(D.spec_fq,nanmean(ALL.PSD{iF}(intv,:,:),3))
        caxis_all(cnt,:) = caxis;
        cnt = cnt + 1;

        if iC == 1
            title({ALL.cond{iF} sprintf('%d min',D.ranges_to_analyze_min(intv,:))},'FontSize',8)
        else
            title({sprintf('%d min',D.ranges_to_analyze_min(intv,:))})
        end
        if iF == length(ALL.PSD)
            xlabel('Freq (Hz)')
        end
    end
end
sgtitle('PSD')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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
    set(gca,'YLim',[O.CFC{ii}.high_fq_range(1) 160])
    %     coSlormap(viridis)
end
sgtitle('CFC')
% equalize_color_axes(aa)

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
imagesc(O.spec_fq,O.spec_fq,corr(O.sz))
axis xy; colorbar
xlabel('Hz')
ylabel('Hz')

end


