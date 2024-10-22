function tsa = ts(t)

% ts  Creates ts object (contains sequence of timestamps as might be stored in an NSMA t-file).
%
% tsa = ts(t)
%
% INPUTS: 
%       t - vector of timestamps
% OUTPUTS: 
%       tsa - a ts object
%
% Methods
%    ts/Data         - Returns the timestamps as a matlab array
%    ts/StartTime    - First timestamp
%    ts/EndTime      - Last timestamp
%    ts/Restrict      - Keep timestamps within a certain range
%
% ADR 1998, version L4.0, last modified '98 by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m


if nargin == 0
   ts.t = [];
   return
end

tsa.t = t;
tsa = class(tsa, 'ts');