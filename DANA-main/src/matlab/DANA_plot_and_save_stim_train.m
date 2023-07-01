function DANA_plot_and_save_stim_train(INFO, fname, PLOT_IT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots all of the information hopefully that we'll need for the
% experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    fname = [];
end
if nargin < 3
    PLOT_IT = true;
end
tot_time_sec = INFO.orig_timestamps_sec(end) - INFO.orig_timestamps_sec(1);
stim_rate = length(INFO.orig_timestamps_sec)/tot_time_sec;

if PLOT_IT
    DANA_plot_stim_train(INFO)
end
if ~isempty(fname)
    saveas(gcf,[fname '.fig']) % save the figure - will be useful.
    % Write a simple text file.
    fileID = fopen([fname '.txt'],'w');
    fprintf(fileID,'%.15f\n',INFO.clean_timestamps_sec); % need to do this %.16f or it gets rid of the precsion!!!
    fclose(fileID);
    save([fname '.mat'],'INFO')
end
