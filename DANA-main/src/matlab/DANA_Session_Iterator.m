%%
clearvars;
close all;
Analysis_function = @Q3_does_stim_affect_ensemble_act;

gd = Git_dir;
% fullfile(gd,'DANA\Data\Acute\Processed_Data\Rat000\20221101') % I sorted
% this last year and now I hate my sorting - ignore until re-sorted.
session_dirs = {   fullfile(gd,'DANA\Data\Acute\Processed_Data\Rat445\01')
    fullfile(gd,'DANA\Data\Acute\Processed_Data\Rat425\01')
    fullfile(gd,'DANA\Data\Acute\Processed_Data\Rat438\01')
    fullfile(gd,'DANA\Data\Acute\Processed_Data\Rat439\01')};

Analysis_function_string = func2str(Analysis_function); % The name of the function called.

output_directory = fullfile(Analysis_dir ,Analysis_function_string);
if ~exist(output_directory,'dir')
    mkdir(output_directory)
end

for iSes = 1:length(session_dirs)
    cd(session_dirs{iSes})
    out_fname = fullfile(output_directory,sprintf('Dset_%d.mat',iSes));

    fprintf('%s\n',out_fname);
    
    Dset = Analysis_function(); % Run the function (assumes it saves all of the relevant data and images).
    Dset.iSes = iSes;
    save(out_fname,'-struct','Dset','-v7.3')
end
