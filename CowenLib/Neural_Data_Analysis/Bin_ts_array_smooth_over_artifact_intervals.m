function [Q,bin_centers] = Bin_ts_array_smooth_over_artifact_intervals(t_cell_array, small_bin_size, large_bin_size, artifact_intervals)
% Geneates a Q matrix in 2 steps to remove (as well as can be) bins that
% have artifact.
% % Cowen 2023
start_time = Find_start(t_cell_array);
end_time = Find_end(t_cell_array);
edges = start_time:small_bin_size:end_time;
if mod(large_bin_size,small_bin_size)~=0
    error('large_bin_size must be a multiple of the small bin size')
end
win_size_bins = large_bin_size/small_bin_size;

for iC = 1:length(t_cell_array)
    [hc(:,iC),x] = histcounts(t_cell_array{iC}, edges);
end
x = x(:);
xctr = x(1:end-1) + diff(x)/2;
% Find bins in x within the artifact window.
BIX = false(size(xctr));
ix = binsearch_vector(xctr,artifact_intervals(:,1)); % SOO much faster than find or logicals
ix2 = binsearch_vector(xctr,artifact_intervals(:,2));
for ii = 1:length(ix)
    BIX(ix(ii):ix2(ii)) = true;
end
hc(BIX,:) = nan;
Qbig = movmean(hc,win_size_bins,'omitmissing');
bin_centers_big = movmean(xctr,win_size_bins);
% Since this has been 'smoothed' we can take the bin center estimate for
% evern n bins and be done in record time.
% Note: don't do decimate - filtering messes things up.
Q = Qbig(1:win_size_bins:end,:);
bin_centers = bin_centers_big(1:win_size_bins:end);



