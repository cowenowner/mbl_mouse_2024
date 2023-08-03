function [out_times, out_labels] = TSToTable(ts_in)
% function [out_times, out_labels] = TSToTable(ts_in)
%
% converts a ts input to two vectors of the same length:
% out_times: ordered list of all times in ts_in
% out_labels: idxs into ts_in.t{idx} of event times
%
% EXAMPLE
%
% ts_in.t{1} = [1 3]; ts_in.t{2} = [2 4];
%
% out_times: [1 2 3 4]
% out_labels: [1 2 1 2]
%
% MvdM 2023

nT = length(ts_in.t);

out_times = [];
out_labels = [];

for iT = 1:nT

    this_t = ts_in.t{iT};

    out_times = cat(1, out_times, this_t);
    out_labels = cat(1, out_labels, iT.*ones(size(this_t)));

end

[out_times, sort_idx] = sort(out_times, 'ascend');
out_labels = out_labels(sort_idx);