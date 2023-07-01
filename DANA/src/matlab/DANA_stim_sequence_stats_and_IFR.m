function STATS = DANA_stim_sequence_stats_and_IFR(t_sec, varargin)
smooth_bin_sec = .05;
n_bins = 10;
Extract_varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = diff(t_sec);
STATS.mean_fr = length(t_sec)/(t_sec(end) - t_sec(1));
STATS.lv = LocalVariance(d);
STATS.cv = std(d)/mean(d);
[Q,Qs,bin_edges] = Instantaneous_Firing_Rate(t_sec, smooth_bin_sec, n_bins,false);
STATS.IFR = Q;
STATS.IFR_smoothed = Qs;
STATS.IFR_timestamps_sec = bin_edges(1:end-1) + median(diff(bin_edges))/2;
