function tc = QtoTC_byCue(Q_in)
% function tc = QtoTC_byCue(Q_in)
%
% computes tuning curv
%
% INPUTS:
% Q_in: Q_matrix TSD (from MakeQFromS()) with usr fields cue, bin and trial
%
% OUTPUTS:
% tc{nCues}: cell array of length nCues containing [nCells x nBins] tuning curves
%
% MvdM Aug 2023, neuropixels odor decoding pilot

cues = unique(Q_in.usr.cue);
bins = unique(Q_in.usr.bin);

nBins = max(bins);
nCells = size(Q_in.data, 1);

for iC = 1:length(cues)

    tc{iC} = zeros(nCells, nBins);

    this_cue = cues(iC);
    this_Q = SelectTSD([], Q_in, Q_in.usr.cue == this_cue);

    trials = unique(this_Q.usr.trial);

    for iT = 1:length(trials)

        this_trial = trials(iT);
        temp = SelectTSD([], this_Q, this_Q.usr.trial == this_trial);

        tc{iC} = tc{iC} + temp.data; % could insert extra step to ensure bins are in ascending order

    end

    tc{iC} = tc{iC} ./ length(trials);

end