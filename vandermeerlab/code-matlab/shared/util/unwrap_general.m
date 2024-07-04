function out = unwrap_general(in, low, high)
% function out = unwrap_general(in, low, high)
% 
% unwraps circular data [in] based on [low] and [high] wraparound points
%
% e.g. out = unwrap(in, -180, 180);

if ~isvector(in)
    error('Input must be 1-D')
end

out = nan(size(in)); out(1) = 0;

range = high - low;

diffs = diff(in);

diff_low_idx = find(diffs <= -range/2);
diffs(diff_low_idx) = diffs(diff_low_idx) + range;

diff_high_idx = find(diffs > range/2);
diffs(diff_high_idx) = diffs(diff_high_idx) - range;

out(2:end) = diffs;
out = cumsum(out);