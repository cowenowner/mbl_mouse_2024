function R = convert(R,unit,unitflag)
% ctsd/Range
% 
%
%  Converts t from unit to unitflag
%     unit       units of t
%     unitflag:  if 'ts' returns time in timestamps,
% JCJ

switch (unitflag)
    case 'sec'
        switch unit
            case 'sec'
                R = R;
            case 'ts'
                R = R/10000;
            case 'ms'
                R = R/1000;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
        
    case 'sec0'
        switch unit
            case 'sec'
                R = (R - min(R));
            case 'ts'
                R = (R - min(R))/10000;
            case 'ms'
                R = (R - min(R))/1000;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
        
    case 'ts'
        switch unit
            case 'sec'
                R = R*10000;
            case 'ts'
                R = R;
            case 'ms'
                R = R*10;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
        
    case 'ms'
        switch unit
            case 'sec'
                R = R*1000;
            case 'ts'
                R = R/10;
            case 'ms'
                R = R;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
    otherwise
        warning('Convert called with invalid unitflag: no conversion possible');
end