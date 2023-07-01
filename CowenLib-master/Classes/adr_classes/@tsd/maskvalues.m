function mTSD = MaskValues(tsa, logic_op, val);
% mTSD = MaskValues(tsd, logic_op, val);
%
% Create a new tsd that only contains data specified by a logical
% operation and a value. For instance MaskValues(tsd,'>',100) will
% mask all values in the tsd that are greater than 100. Masking
% involves substituting nans for values.
%
% INPUT:  a tsd to be masked
%         a logical operation as a character string 
%             '>' 
%             '<' 
%             '==' 
%             '~='  
%
%         a criteria value to mask. Must be a scalar.
%
% OUTPUT: a new tsd that has the values speficied by the logical
%         operation masked out
%
%
% cowen Sat Jun 26 16:05:42 1999
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

% Construction output tsd
mTSD = tsa;

switch logic_op
  case '>'
    idx = find(tsa.data > val);
  case '<'
    idx = find(tsa.data < val);
  case '>='
    idx = find(tsa.data >= val);
  case '<='
    idx = find(tsa.data <= val);
  case '=='
    idx = find(tsa.data == val);
  case '~='
    idx = find(tsa.data ~= val);
  otherwise
    error('invalid operator')
end

% Set the masked values to nan.
mTSD.data(idx) = NaN;

