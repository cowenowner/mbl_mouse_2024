function R = Correlate_spikes_to_param_by_kernel_width(TS, time_val, base_bin_size, smooth_kernel_bins, varargin)
% Does not look at single neurons unless only a single neuron is passed in
% Will do PCA on the ensemble.
%
% parametrically calculate the correlation between spiking and some
% continuous variabls (time_val) for different time window (smooth_kernel_bins) sizes.
% A moving average is calculated with a window size according to smooth_kernel_bins.
%
% Time units are assumed to be consistent for all input vals - so units do
% not matter as long as they are the same for TS time_val smooth_kernel_bins.
%
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
time_val = double(time_val);
R.r_v_neuron = nan(length(TS), length(smooth_kernel_bins));
R.r_v2_neuron = nan(length(TS), length(smooth_kernel_bins));
R.r_v_pc = nan(2, length(smooth_kernel_bins));
R.r_v2_pc = nan(2, length(smooth_kernel_bins));
[Q, edges, Qt] = histcounts_cowen(TS,'binsize',base_bin_size);
% Interp method -  time resolution of the motion stays constant
% but the spike time resolution changes
v = interp1(time_val(:,1) , time_val(:,2), Qt);
% % Abhi method - both the time resolution of the spikes and the position
% % change together.
% vv = nan(length(edge_ix)-1,1);
% for iX = 1:(length(edge_ix)-1)
%     vv(iX) = mean(time_val(edge_ix(iX):edge_ix(iX+1),2));
% end


for iBinWidth = 1:length(smooth_kernel_bins)
    % Qs = conv_filter(Q, hanning(smooth_kernel_bins(iBinWidth)));
    Qs = movmean(Q, smooth_kernel_bins(iBinWidth));
    
    Qsz = Z_scores(Qs);
    [~,sc] = pca(Qsz);
    % edge_ix = binsearch_vector(time_val(:,1),edges);

    for iN = 1:Cols(Q)
        [R.r_v_neuron(iN,iBinWidth),R.r_v_neuron_p(iN,iBinWidth)] = corr(Qsz(:,iN),v,'Type',corr_type,'Rows','complete');
        % [R.r_v2_neuron(iN,iBinWidth),R.r_v2_neuron_p(iN,iBinWidth)] = corr(Qsz(:,iN),vv,'Type',corr_type,'Rows','complete');
    end

    if Cols(Q) > 2 && Rows(Q)>20 && sum(sum(abs(Q)) >0) > 2
        R.r_v_pc(1,iBinWidth) = corr(sc(:,1),v,'Type',corr_type,'Rows','complete');
        % R.r_v2_pc(1,iBinWidth) = corr(sc(:,1),vv,'Type',corr_type,'Rows','complete');
        R.r_v_pc(2,iBinWidth) = corr(sc(:,2),v,'Type',corr_type,'Rows','complete');
        % R.r_v2_pc(2,iBinWidth) = corr(sc(:,2),vv,'Type',corr_type,'Rows','complete');
    end
end

if PLOT_IT
    figure
    subplot(1,2,1)
    imagesc(R.r_v_neuron.*(R.r_v_neuron_p<0.05));colorbar
    xlabel('window width'); ylabel('neuron')
    yyaxis right
    hold on
    plot(mean(R.r_v_neuron,'omitmissing'),'w','LineWidth',4)
    subplot(1,2,2)
    ylabel('r')
    plot(smooth_kernel_bins,R.r_v_pc')
    yyaxis right
    plot(smooth_kernel_bins,R.r_v_pc.^2')
    ylabel('r2')

    xlabel('Bin width'); legend('pc1','pc2')
end