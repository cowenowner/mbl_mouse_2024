function OUT = Fit_to_spline_by_trial(D,x,trial_start_end,polyorder)
% INPUT: D = 3 cols 1) time, 2)data to fit, 3) xposition
% x = desired x for fitting.
% ST_ED = Start and end for each trial.

OUT = nan(Rows(trial_start_end),length(x));

for iTrial = 1:Rows(trial_start_end)
    Dtmp = Restrict(D,trial_start_end(iTrial,:));
    Dtmp = Dtmp(~isnan(sum(Dtmp,2)),:);
    Dtmp = sortrows(Dtmp,3);
    
    warning off
    p = polyfit(Dtmp(:,3),Dtmp(:,2),polyorder);
    warning on
    OUT(iTrial,:) = polyval(p,x);
end