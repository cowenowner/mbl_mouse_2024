function [R] = Correlate_neurons_to_neurons_multiple_bin_widths(TS, base_bin_size, smooth_kernel_bins, varargin)
% TS - cell array of timestamsp.
% time_val = 2 col where 1st col is time, second is the param to correlate
% smooth_kernel_bins - movmean widths to compare.
% base_bin_size = the baseline bin size and then the average is computed
% after this base.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr_type = 'Pearson';
% corr_type = 'Spearman'; Can freeze
PLOT_IT = false;
Extract_varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = [];
st = Find_start(TS);
ed = Find_end(TS);
RIX = triu(ones(length(TS)),1)>0;
[R.rowix,R.colix] = ind2sub(size(RIX), find(RIX>0));
R.R = nan(sum(RIX(:)),length(smooth_kernel_bins));
R.P = nan(sum(RIX(:)),length(smooth_kernel_bins));

if length(TS) < 3 || isnan(st) || st==ed
    return
end


Q = histcounts_cowen(TS,'binsize',base_bin_size);
for iBinWidth = 1:length(smooth_kernel_bins)
    Qs = movmean(Q, smooth_kernel_bins(iBinWidth));
    % Qsz = Z_scores(Qs);
    [Rtmp,Ptmp] = corr(Qs,'Type',corr_type,'Rows','complete');
    R.R(:,iBinWidth) = Rtmp(RIX);
    R.P(:,iBinWidth) = Ptmp(RIX);
    % Rz = corr(Qsz,'Type',corr_type,'Rows','complete');
end

if PLOT_IT
    figure
    subplot(1,2,1)
    imagesc(R.R);
    caxis([0 .5])
    xlabel('window width'); colorbar;ylabel('cell pair')
    yyaxis right
    plot(mean(R.R,1,'omitnan'),'w')
    title('R')
    subplot(1,2,2)
    imagesc(R.P);colorbar
    xlabel('Bin width'); 
    title('P')
end