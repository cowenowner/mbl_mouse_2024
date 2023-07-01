%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q4_Do_spikes_phase_lock_to_LFP_Ana';
ses_to_ana = 'Q4_Do_spikes_phase_lock_to_LFP';

GP = LK_Globals;
GP.Analysis_Dir = Analysis_dir;

PLOT_IT = false;
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
cm = lines(5);
TBL = table();

for iSes = 1:length(d)
    Dset = load(fullfile(adir,d(iSes).name));
    TBL = [TBL; Dset.TBL];
end
IX = TBL.Condition == 'LDO' & TBL.FreqBand == 'gamma_80' & TBL.Interval == 2 & TBL.Depth_uM < 2000;
sum(IX)
sum(TBL.Ang_p(IX) < 0.05)
sum(TBL.powAng_p(IX) < 0.05) % this is too high - very strange. I need to dig deeper.
