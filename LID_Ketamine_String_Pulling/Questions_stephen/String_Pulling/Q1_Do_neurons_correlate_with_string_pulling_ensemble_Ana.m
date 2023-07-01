function OUT = Q1_Do_neurons_correlate_with_string_pulling_ensemble_Ana()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run this function in the session directory.
% Determine if single-units oscillate more after ketamine injection.
% Compare autocorrs before and after injection.
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = [];
close all
% global DIRS % Assumes this variable was set up already. Do this in LK_Analyze_sessions (or somewhere).
GP = LK_Globals;
PLOT_IT = true;
mfile = 'Q1_Do_neurons_correlate_with_string_pulling_ensemble';
ana_dir = 'C:\Users\Stephen Cowen\Dropbox\Foldershare\Analysis_Results_Dropbox\Q1_Do_neurons_correlate_with_string_pulling_ensemble';
d = dir(fullfile(ana_dir,'Dset*.mat'));
ALL = [];
Rsq_norm = [];
sescnt = 1;
%%%%%%%%%%%%%%%%%%%%

for iF = 1:length(d)
    Dset = load(fullfile(ana_dir,d(iF).name));
    if Dset.aborted
        continue
    end
    XIX = Dset.smooth_win_ms < 500 & Dset.smooth_win_ms > 50;
    
    ALL(sescnt).rat = Dset.SES.rat;
    ALL(sescnt).session = Dset.SES.session;
    ALL(sescnt).title_str = Dset.SES.title_str;
    ALL(sescnt).n_neurons = Dset.n_neurons;

    tmp = Dset.Rsq - Dset.Rsq_sh;
    tmp(Dset.Rsq_p > 0.05) = nan;
    tmp(tmp < 0.001) = nan; % If it's less than shuffle- clearly not significant.
    
    ALL(sescnt).mean_Rsq500 = nanmean(tmp(:,XIX),2);

    Rsq_norm(:,:,sescnt) = tmp;  
    
    sescnt = sescnt + 1;
end
if PLOT_IT
    for iSes =  1:length(ALL)
        figure
        imagesc(Dset.smooth_win_ms,[],squeeze(Rsq_norm(:,:,iSes)))
        set(gca,'YTick',1:length(Dset.variables))
        set(gca,'YTickLabel',Dset.variables)
        colorbar_label('R^2-shuff')
        title(sprintf('%s n=%d',ALL(iSes).title_str,ALL(iSes).n_neurons))
    end
    figure
    bar([ALL.mean_Rsq500]')
    legend(Dset.variables); legend boxoff
    set(gca,'XTickLabel',[ALL.session])
    xlabel('Session')
    ylabel('Rsq < 500ms smooth win')
    pubify_figure_axis
end
% This converts the structre above to a useful table.
% TBL = struct2table(ALL);
%% Summarize the statistics.
% Let's assume that a binsize < 500 ms is reasonable. Bigger bin sizes make
% me think the neuron is responding to some far more long-term change in
% state (which could be interesting) - maybe need to separate the two.


