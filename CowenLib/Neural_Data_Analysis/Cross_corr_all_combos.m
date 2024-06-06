function [cc,x] = Cross_corr_all_combos(T, binsize, n_lags, varargin)
% function [ac,x] = Cross_corr(ts1, ts2, binsize, n_lags, varargin)
%
% Standard cross-correlogram for spike trains.
%
% INPUT: T = cell array of timestamp vectors (sorted)
%        binsize = binsize for binning timestamps. (same units as ts)
%        n_lags = n_lags in the xcorr. 
%
% OUTPUT: cc = cross correlation 
%         x = the x axis (whatever units ts and binsize are in).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_it = false;
scaleopt = 'none'; %'coeff'; 
include_autocorrs = true;
fast_crosscorr = true;
% scaleopt = 
Extract_varargin

if include_autocorrs
    shift = 0;
else
    shift = 1;
end

ONES = ones(length(T));
UPPER = triu(ONES,shift) > 0;
[i,j]=ind2sub(size(UPPER),find(UPPER));
cc = [];
for ii = 1:length(i)
    if fast_crosscorr
        [cc(ii,:),x] = CrossCorr(T{i(ii)},T{j(ii)}, binsize, n_lags*2); % SUPER FAST
        % [cc1,x1] = CrossCorrCount(T{i(ii)},T{j(ii)}, binsize, n_lags*2); % SUPER FAST
    else
        [cc(ii,:),x] = Cross_corr(T{i(ii)},T{j(ii)}, binsize, n_lags,'scaleopt',scaleopt); %SSLOOW
        % [cc2,x2] = Cross_corr(T{i(ii)},T{j(ii)}, binsize, n_lags,'scaleopt',scaleopt); %SSLOOW
    end
end

if plot_it
    % figure
    imagesc(x,[],cc)
    axis tight
    box off
    plot_vert_line_at_zero
    set(gca,'FontSize',8)
    % ylabel([scaleopt])
end