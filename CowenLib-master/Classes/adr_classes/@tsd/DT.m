function dt = DT(tsa, tsflag)

% tsd/DT  Returns DT (mean timestep between sequential entries) from tsd
%
% dt = DT(tsa, tsflag)
%	
% INPUTS:
%       tsa - tsd object
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'ms' returns time in ms
% OUTPUTS:
%       dt - mean timestep between sequential entries in tsa
%
% ADR, version L1.0, last modified by ADR
% cowen (2009): Switched to median - I don't trust mean as outliers screw things
% up if you use segmented data.
% 
% status: promoted


dt = median(diff(tsa.t));

if nargin == 2
   switch tsflag
   case 'sec'
      dt = dt/10000;
   case 'ms'
      dt = dt/10;
   case 'ts'
   otherwise
      error('Unknown tsflag.');
   end
end