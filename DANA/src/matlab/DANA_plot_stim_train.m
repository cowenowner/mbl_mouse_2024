function DANA_plot_stim_train(INFO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots all of the information hopefully that we'll need for the
% experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_time_sec = INFO.orig_timestamps_sec(end) - INFO.orig_timestamps_sec(1);
stim_rate = length(INFO.orig_timestamps_sec)/tot_time_sec;

figure
subplot(3,2,1:2)
plot(INFO.orig_timestamps_sec, zeros(size(INFO.orig_timestamps_sec)), '+',INFO.clean_timestamps_sec, zeros(size(INFO.clean_timestamps_sec)), 'o' )
hold on
axis tight
plot_markers_simple(INFO.blank_intervals,[],1)
yyaxis right
plot(INFO.clean_stats.IFR_timestamps_sec,INFO.clean_stats.IFR_smoothed )
ylabel('IFR')
% plot(INFO.blank_intervals(:,2), zeros(size(INFO.blank_intervals(:,2)))-1,'r<')
xlabel('sec')
pubify_figure_axis
title(sprintf('LV: %1.2f, preLV: %1.2f, nspk %d, %2.1f Hz, CV: %1.2f, pre CV: %1.2f ', INFO.clean_stats.lv, INFO.orig_stats.lv, ...
    length(INFO.orig_timestamps_sec), stim_rate, INFO.clean_stats.cv, INFO.orig_stats.cv))
subplot(3,2,3); histogram_cowen({diff(INFO.orig_timestamps_sec) diff(INFO.clean_timestamps_sec) },.01)
xlabel('sec'); title('ISI hist pre post')
subplot(3,2,4); Auto_corr(INFO.clean_timestamps_sec, .01, .5/.01,'plot_it',true);
xlabel('sec'); title('AutoCorr'); ylabel('count')
pubify_figure_axis
subplot(3,2,5); 
histogram_cowen({log2(1./diff(INFO.orig_timestamps_sec)) log2(1./diff(INFO.clean_timestamps_sec)) },.2)
xlabel('log2 rate Hz'); title('log Inst Rate hist pre post')
pubify_figure_axis

