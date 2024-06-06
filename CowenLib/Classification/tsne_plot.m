function [Y,Groups] = tsne_plot(M, varargin)
n_clusters = 4;
Groups = [];
Extract_varargin;

if isempty(Groups)
    % [~,Groups] = kmeans_optimal_k(M,7);
     Groups = clusterdata(M,'Linkage','ward','SaveMemory','on','Maxclust',n_clusters);
end

Y = tsne(M);

gscatter(Y(:,1),Y(:,2),Groups);
