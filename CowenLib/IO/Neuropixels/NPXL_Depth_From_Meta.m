function TBL = NPXL_Depth_From_Meta(meta)
% Pass in the meta structure (after running SGLXMetaToCoords()  to add
% depth ot the meta file)
%
% SGLXMetaToCoords is now located here: https://github.com/jenniferColonell/SGLXMetaToCoords
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = meta.snsChanMap;
tmp2 = strrep(tmp,':',',');
tmp2 = strrep(tmp2,';',',');
tmp2 = strrep(tmp2,'AP','');
tmp2 = strrep(tmp2,'LF','');
tmp2 = strrep(tmp2,'SY0','-1');
tmp2 = strrep(tmp2,')(',';');
tmp2 = strrep(tmp2,'(','[');
tmp2 = strrep(tmp2,')',']');
SNS= eval(tmp2);
% I think the last row of SNS is the sync so nothing.
% align with the SNS
SNS(1,:) = [];
SNS(end,:) = [];

tmp = meta.snsGeomMap;
ix = strfind(tmp,')');
tmp = tmp(ix+1:end);
tmp2 = strrep(tmp,':',',');
tmp2 = strrep(tmp2,';',',');
tmp2 = strrep(tmp2,')(',';');
tmp2 = strrep(tmp2,'(','[');
tmp2 = strrep(tmp2,')',']');
% ix = strfind (tmp2,',');
% tmp2(2:(ix(1)-2)) = [];
% tmp2 = strrep(tmp2,'C,',[]);
GEO= eval(tmp2);
% I think the first row of GEO is just the probe info so eliminate to
% align with the SNS
% GEO(1,:) = [];

Ch = SNS(:,1);
ChannelID = SNS(:,2);
SiteID = SNS(:,3);
shank = GEO(:,1);
x_uM = GEO(:,2); % within shank x
y_uM = GEO(:,3); 
kdist = GEO(:,4);
x_inter_shank_uM = x_uM + (shank*250); % presumes 250uM spacing between shanks)

TBL= table(Ch, ChannelID, SiteID,shank,x_uM,y_uM,kdist,x_inter_shank_uM);
TBL.Depth_uM = max(TBL.y_uM) - TBL.y_uM;

