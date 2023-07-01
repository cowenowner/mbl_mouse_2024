function C = MaskByTSD(A,B)
%
% C = ts/MaskByTSD(A,B)
% 
% INPUTS:
%    A = tsd object
%    B = tsd object
%
% OUTPUTS:
%    C = tsd only including times in which nearest tsd sample from B
%          is not NaN
%
% ADR 2000
% v 4.0: JCJ 3/3/2003 includes support for time units

if isempty(strmatch('units',fieldnames(A)))
    warning('units not specified: assuming units = sec (converstions preformed)' )
    unit ='sec';
else
    unit = A.units;
    
end


TIMES = Range(A, unit);
DATAS = Data(A);
nT = length(TIMES);

tTime = Range(B, unit); % both A and B must have same units
tData = Data(B);

for iS = 1:nT
  if isnan(tData(binsearch(tTime, TIMES(iS))))
    TIMES(iS) = nan;
  end
end

f = find(~isnan(TIMES));
C = tsd(TIMES(f), DATAS(f), unit);
