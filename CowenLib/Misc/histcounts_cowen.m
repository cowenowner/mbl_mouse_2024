function [hc,x, xctr] = histcounts_cowen(IN, varargin)
% just do hi stcounts for a cell array - using edges.
edges = [];
nbins = [];
binsize = [];

Extract_varargin;

if ~isempty(binsize)
    if iscell(IN)
        edges = Find_start(IN):binsize:(Find_end(IN)+eps);
    else
        edges = IN(1):binsize:(IN(end)+eps);
    end
end
if ~isempty(edges)
    edges_or_nbins = edges;
end
if ~isempty(nbins)
    edges_or_nbins = nbins;
end


if ~isa(IN,'cell')
    [hc,x] = histcounts(IN, edges_or_nbins);
    return
end

if length(edges_or_nbins) == 1
    hc = zeros(edges_or_nbins,length(IN));
else
    hc = zeros(length(edges_or_nbins)-1,length(IN));
end

for iC = 1:length(IN)
    [hc(:,iC),x] = histcounts(IN{iC}, edges_or_nbins);
end
x = x(:);
if nargout > 2
    xctr = x(1:end-1) + diff(x)/2;
end
if nargout ==0
    figure
    imagesc([],x,hc)
end

