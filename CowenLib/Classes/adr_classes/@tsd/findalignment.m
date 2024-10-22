function ix = findalignment(D, tstmp)
%function ix = findAlignment(D, tstmp)
% tsd/findAlignment  Returns index of specified timestamp in tsd
%
% ix = findAlignment(D, tstmp)
%
% INPUTS:
% 	    D - tsd object
%       tstmp - timestamp that you want to find the index of
% OUTPUTS:
%       ix - index of timestamp in D closest to tstmp
%
% ADR, version L4.0, last modified by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m


nIX = length(tstmp);
ix = zeros(size(tstmp));
for iIX = 1:nIX
   ix(iIX) = binsearch(D.t, tstmp(iIX));
end

