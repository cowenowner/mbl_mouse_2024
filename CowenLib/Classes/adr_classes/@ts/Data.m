function t = Data(TS, ts)

% ts/Data  Returns data (timestamps) in ts object
%
% t = Data(TS, ts)  
% t = Data(TS)      
%
% INPUTS:
%       TS - ts object
%       ts - (if specified) function returns largest timestamp in TS <= ts
% OUTPUTS:
%       t - (if ts specified) the largest timestamp in TS that is <= ts
%           OR
%           (if no ts input) the entire array of data (timestamps) in TS
%
% For ts objects, Data and Range methods are identical
%
% ADR, version L4.2, last modified by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m


if nargin == 1
   t = TS.t;
else
   t = TS.t(binsearch(TS.t, ts));
end