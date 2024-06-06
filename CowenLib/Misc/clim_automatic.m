function clim_automatic(reduction)
% Cowen 2023
if nargin < 1
    reduction = .7;
end
clim(clim*reduction);
