function [spikes_uS_ca, WV_uV] = plot_neurons(spikes_uS_ca, WV_uV, WV_x_ms, varargin)
% Simple plot of all of the neurons collected - wavefomrs and autocorrs.
% some simple clustering.
% Cowen 2023
SORT_ROWS = false;
AC_binsize_ms = 4;
AC_duration_ms = 200;
Groups = [];

Extract_varargin
ID = 1:length(spikes_uS_ca);

if SORT_ROWS
    % Sort it all so it looks nice.
    [WV_uV, ~, six]= sort_matrix(WV_uV);
    tmp = [];
    for iN = 1:length(spikes_uS_ca)
        tmp{iN} = spikes_uS_ca{six(iN)};
    end
    spikes_uS_ca = tmp;
    ID = six;
    if ~isempty(Groups)
        Groups = Groups(six);
    end
end

ac = [];
for iN = 1:length(spikes_uS_ca)
    [tmp,lag_x_ms] = AutoCorr(spikes_uS_ca{iN}/1e3, AC_binsize_ms, AC_duration_ms/AC_binsize_ms); % 100x faster than Cross_cor. Same result.
    ac(iN,:) = movmean(tmp,3);
end
ac_norm = ac./ sum(ac,2,'omitmissing');
ac_norm(isnan(ac_norm)) = 0;


subplot(4,2,1)
imagesc(WV_x_ms,[], WV_uV)
subplot(4,2,2)
imagesc(lag_x_ms,[], ac)
subplot(4,2,3)
plot_confidence_intervals(WV_x_ms,WV_uV)
subplot(4,2,4)
plot_confidence_intervals(lag_x_ms,ac)
subplot(4,2,7)
[Y,C]=tsne_plot(WV_uV,'Groups',Groups);
subplot(4,2,5)
u = unique(C);
clrs = lines(length(u));
for iG = 1:length(u)
    IXG = C==u(iG);
    plot_confidence_intervals(WV_x_ms,WV_uV(IXG,:),[],clrs(iG,:))
end
subplot(4,2,8)
[Ya,Ca]=tsne_plot(ac_norm,'Groups',Groups);
subplot(4,2,6)
u = unique(Ca);
clrs = lines(length(u));
for iG = 1:length(u)
    IXG = Ca==u(iG);
    plot_confidence_intervals(lag_x_ms,ac(IXG,:),[],clrs(iG,:))
end




