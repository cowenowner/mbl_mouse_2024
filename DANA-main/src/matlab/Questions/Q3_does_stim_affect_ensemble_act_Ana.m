% Q3_does_stim_affect_ensemble_act_Ana()
%%
clearvars;
close all
clrs = lines(200);
ana_dir = 'E:\Dropbox\Foldershare\Analysis_Results_Dropbox\Q3_does_stim_affect_ensemble_act';
d = dir(fullfile(ana_dir,'*.mat'));

for iD = 1:length(d)
    Dset = load(fullfile(d(iD).folder,d(iD).name));
    if iD == 1
        TBL = Dset.TBL;
    else
        TBL = [TBL; Dset.TBL];
    end
end
% Sanity check for depth...

% figure;plot(TBL.Depth_corrected_uM);ylabel('uM from surface')

% Cluster based on waveform shape.
% This is messed up - some waveforms look strange.
% Some depths are not consistent as sometimes the depth was entered
% correctly, sometimes not.
if 1
    % WV = Z_scores(TBL.WV')';
    % norm so that trough = -1;
    WV = TBL.WV;
    for iR = 1:Rows(WV)
        WV(iR,:) = WV(iR,:) / abs(min( WV(iR,:) ));
    end
    WV(isnan(WV)) = 0;
    F = Z_scores([WV]); % My favorite.
    [eva, C] = kmeans_optimal_k(F,9);
    TBL.CellTypeKM = C;
    [~,sc] = pca(F);
    figure
    uC = unique(C);
    gix = 30:60;
    good_cluster = ones(length(uC),1);
    for iU = 1:length(uC)
        mwv = mean(WV(C==uC(iU),gix));
        % mwv = mean(TBL.WV(C==uC(iU),gix))
        if max(mwv) > abs(min(mwv))
            good_cluster(iU) = 0;
        end
        plot(mwv,'LineWidth',4)
        hold on
    end
    axis off
    axis tight
    good_cluster_IDs = find(good_cluster==1)';
    figure
    % Just plot the good ones
    for iU = 1:length(good_cluster_IDs)
        mwv = mean(WV(C==good_cluster_IDs(iU),gix));
        plot(mwv,'LineWidth',4)
        hold on
    end
    axis off
    axis tight

    figure
    hist(C,1:3)
    pubify_figure_axis


    figure
    gscatter(sc(:,1),log(TBL.PeakRatesHz),C); xlabel('pc1');ylabel('Log Peak Rate')
    pubify_figure_axis
end
% save('C:\Temp\Q3_DANA.mat') % use this for the BRAIN poster. Need to clean up
% kmeans to be more consistent and based onw waveform width.
%%
tmp_clrs = [1 .2 .2;.2 .2 1];
uHz = unique(TBL.Hz);
uLVs = unique(TBL.LV);
smth_bin = 5;
TBL.mean_PETHnorm_smth = TBL.mean_PETHnorm;
for iR = 1:Rows(TBL)
    TBL.mean_PETHnorm_smth(iR,:) = conv(TBL.mean_PETHnorm(iR,:),hanning(smth_bin)/sum(hanning(smth_bin)),'same');
end


figure
hg = [];cnt = 1;

for iN = 1:length(good_cluster_IDs)
    hg(cnt) = subplot(1,length(good_cluster_IDs),iN);
    for ii = 1:length(uHz)
        IX = TBL.Hz == uHz(ii) & TBL.CellTypeKM == good_cluster_IDs(iN);

        plot_confidence_intervals(Dset.PETH_x_axis_ms/1000, TBL.mean_PETHnorm_smth(IX,:),[],tmp_clrs(ii,:))
        length(unique(TBL.CellID(IX)))
    end
    title(sprintf('Clust %d', good_cluster_IDs(iN)))
    pubify_figure_axis
    axis tight
    plot_horiz_line_at_zero;
    if iN == 1
        legend_color_text({'10' '20'},tmp_clrs(1:2,:))

    end
    ylabel('Norm rate'); xlabel('sec')
    cnt = cnt + 1;
end
% equalize_axes(hg)
for ii = 1:length(hg)
    subplot(hg(ii))
    plot_vert_line_at_zero; plot_vert_line_at_zero(-10)

end
set(gcf,'Position',[240         293        1172         465])


figure
hg = [];

for iN = 1:length(good_cluster_IDs)

    hg(iN) = subplot(1,length(good_cluster_IDs),iN);

    for ii = 1:length(uLVs)
        IX = TBL.LV == uLVs(ii);  %TBL.CellTypeKM
        IX = TBL.LV == uLVs(ii) & TBL.CellTypeKM == good_cluster_IDs(iN);
        plot_confidence_intervals(Dset.PETH_x_axis_ms/1000, TBL.mean_PETHnorm_smth(IX,:),[],clrs(ii,:))
    end
    title(sprintf('LV Clust %d', good_cluster_IDs(iN)))

    ylabel('Norm rate'); xlabel('sec')
    pubify_figure_axis
    axis tight
    plot_vert_line_at_zero; plot_vert_line_at_zero(-10)
    plot_horiz_line_at_zero;
    legend_color_text({'0' '.3' '1' '1.2'},clrs(1:4,:))
end
for ii = 1:length(hg)
    subplot(hg(ii))
    plot_vert_line_at_zero; plot_vert_line_at_zero(-10)

end
set(gcf,'Position',[240         293        1172         465])




figure
for iHz = 1:length(uHz)
    h(iHz) = subplot(1,length(uHz),iHz);
    for iLV = 1:length(uLVs)
        IX = TBL.LV == uLVs(iLV) & TBL.Hz == uHz(iHz);
        plot_confidence_intervals(Dset.INFO.PETH_x_axis_ms/1000, TBL.mean_PETHnorm(IX,:),[],clrs(iLV,:))
    end
    pubify_figure_axis
    axis tight
    title(sprintf('%d Hz',uHz(iHz))); ylabel('Norm rate'); xlabel('sec')
end
equalize_axes(h);
subplot(h(1))
plot_vert_line_at_zero; plot_vert_line_at_zero(-10)
plot_horiz_line_at_zero;
subplot(h(2))
plot_vert_line_at_zero; plot_vert_line_at_zero(-10)
plot_horiz_line_at_zero;
legend_color_text({'0' '.3' '1' '1.2'},clrs(1:4,:))
set(gcf,'Position',[240         293        1172         465])



figure
h = [];
cnt = 1;
for iHz = 1:length(uHz)
    for iLV = 1:length(uLVs)
        h(cnt) = subplot_ij(length(uLVs),length(uHz),iLV,iHz);
        IX = TBL.LV == uLVs(iLV) & TBL.Hz == uHz(iHz);
        imagesc(Dset.INFO.PETH_x_axis_ms/1000, [], TBL.mean_PETHnorm(IX,:))
        caxis([-2 4])
        plot_vert_line_at_zero
        plot_vert_line_at_zero(-10)
        cnt = cnt + 1;
        title(sprintf('%d Hz %1.2f LV',uHz(iHz),uLVs(iLV))); ylabel('Neuron'); xlabel('sec')
    end
    pubify_figure_axis
end

pubify_figure_axis
axis tight
plot_vert_line_at_zero; plot_vert_line_at_zero(-10)
plot_horiz_line_at_zero;
legend_color_text({'0' '.3' '1' '1.2'},clrs(1:4,:))



hg = []; pp = [];
figure
labs = {'1-15' '15-29' '29-34'};
for iN = 1:length(good_cluster_IDs)

    for ii = 1:length(uHz)
        IX = TBL.Hz == uHz(ii) & TBL.CellTypeKM == good_cluster_IDs(iN);
        M = TBL.LocalVarianceAroundStim(IX,:);
        M = M - M(:,1);
        subplot_ij(2,length(uHz),iN,ii)
        % boxplot(M(:,2:end),'notch','on')
        Boxplot_points(M(:,2:end))
        title(sprintf('WVClu %d  %d Hz',good_cluster_IDs(iN),uHz(ii)))
        ylabel('LV');set(gca,'XTickLabel',labs)
        pubify_figure_axis
        plot_horiz_line_at_zero
        set(gca,'YLim', [-1 1])
        for iP = 2:4
            [~,pp(iN,iP,ii)] = ttest(M(:,iP));
        end
    end

end
set(gcf,'Position',[363   195   918   614])

figure

for iN = 1:length(good_cluster_IDs)

    for ii = 1:length(uHz)
        IX = TBL.Hz == uHz(ii) & TBL.CellTypeKM == good_cluster_IDs(iN);;
        M = TBL.CVAroundStim(IX,:);
        M = M - M(:,1);
        subplot_ij(2,length(uHz),iN,ii)
        boxplot(M,'notch','on')
        ylabel('CV')
        % title(sprintf('%d Hz',uHz(ii)))
        title(sprintf('WVClu %d  %d Hz',good_cluster_IDs(iN),uHz(ii)))
        pubify_figure_axis
        plot_horiz_line_at_zero
    end
end
set(gcf,'Position',[363   195   918   614])

figure
for iN = 1:length(good_cluster_IDs)
    for ii = 1:length(uLVs)
        IX = TBL.LV == uLVs(ii);
        M = TBL.LocalVarianceAroundStim(IX,:);
        M = M - M(:,1);

        subplot_ij(2,length(uLVs),iN,ii)
        boxplot(M,'notch','on')
        ylabel('LV')
        % title(sprintf('%1.2f LV',uLVs(ii)))
        title(sprintf('WVClu %d  %1.2f LV',good_cluster_IDs(iN),uLVs(ii)))
        plot_horiz_line_at_zero
    end
end
set(gcf,'Position',[363   195   918   614])

figure
for ii = 1:length(uLVs)
    IX = TBL.LV == uLVs(ii);
    M = TBL.CVAroundStim(IX,:);
    M = M - M(:,1);
    subplot(1,length(uLVs),ii)
    boxplot(M,'notch','on')
    ylabel('CV')
    title(sprintf('%1.2f LV',uLVs(ii)))
end




