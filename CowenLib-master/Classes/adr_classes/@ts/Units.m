function R = Units(tsa,tUnits)
%
% R = Units(tsa)
%	returns units of tsa
%
% R = Units(tsa,tUnits)
%   assigns tUnits to tsa, if units exist converts to tUnits
%
% version L1.0
% 
% Status: IN PROGRESS


switch nargin
    case 1
        if isempty(strmatch('units',fieldnames(tsa)))
            warning('tsa units not specified' )
            R=[];
        else
            R=tsa.units;
        end
    case 2
        if isempty(strmatch('units',fieldnames(tsa)))
            R=ts(range(tsa),tUnits);
        else
            R=ts(range(tsa,tUnits),tUnits);
        end
        

    otherwise
        error('Unknown number of input arguments.');

end