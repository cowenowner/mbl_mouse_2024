function S = Fill_Nans_with_interpolants(S)
% For each column, find the nans and replace with a linear interpolant.
S = double(S);
for iC = 1:Cols(S)
    IX = isnan(S(:,iC));
    ix = find(IX);
    ix_to_get = find(~IX);
    if ~isempty(ix_to_get)
        S(ix,iC) = interp1(ix_to_get,S(ix_to_get,iC),ix);
    end
end
% the above works fine if there are Nans in the middle of the data but not
% if there are nans at the start and end. The following fixes this
for iC = 1:Cols(S)
    nanIX = isnan(S(:,iC));
    if any(nanIX)
        st_ix = find(~isnan(S(:,iC)),1,'first');
        if ~isempty(st_ix)
            S(1:(st_ix-1),iC) = S(st_ix,iC);
        end
    end
    nanIX = isnan(S(:,iC));
    if any(nanIX)
        st_ix = find(isnan(S(:,iC)),1,'first');
        if ~isempty(st_ix)
            S(st_ix:end,iC) = S(st_ix-1,iC);
        end
    end
end