function S = ConvertSPtoS(SP)
% function S_ts = ConvertNPXtoMVDMLAB(SP_in)

nCells = length(SP);

% convert to s and align with LFPs
for iC = 1:nCells

    SP(iC).t = SP(iC).t_uS * 10^-6;

end

% sort by depth along probe
depths = [SP(:).neuropixels_depth_uM];
[sorted_depths, sort_idx] = sort(depths, 'descend');
SP = SP(sort_idx);


S = ts;
for iC = nCells:-1:1

    S.t{iC} = SP(iC).t;
    S.label{iC} = SP(iC).cluster_id;
    S.usr.depth(iC) = sorted_depths(iC);

end