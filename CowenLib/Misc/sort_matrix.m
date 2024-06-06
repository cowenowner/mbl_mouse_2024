function [OUT, v, six]= sort_matrix(M,sort_type,param)
% Resorts matrix by peak or trough or pdist. 
% INPUT: M = matrix to sort (row = sample, col = param)
%        sort_type = how to sort('peak')
%        param = 1) used for kmeans - number of clusters OR, of sort_type = 'by_value' 
%        then this is the value upon which the matrix is sorted.
%
% OUTPUT: resorted matrix.
% Cowen 2021 - updated the kmeans to work with big-column matrices.
if nargin < 2
    sort_type = 'peak_loc'; % location of the peak
    sort_type = 'pc1'; % location of the peak
end
if nargin < 3
    param = [];
end

% sort_val2 = ones(Rows(M),1); % A second sort index if you need it.
[~,sort_val2] = max(M,[],2,'omitmissing');
switch(lower(sort_type))
    case {'mean'}
        [sort_val] = mean(M,2,'omitmissing');
    case {'peak_loc' 'max_loc'}
        [~,sort_val] = max(M,[],2);   
    case {'trough_loc' 'min_loc'}
        [~,sort_val] = min(M,[],2);   
    case {'peak' 'max'}
        [sort_val] = max(M,[],2);
    case {'by_value'}
        sort_val = param;
    case {'trough' 'min'}
        [sort_val] = min(M,[],2,'omitmissing');
    case 'sparsity'
        sort_val = Density(M')';
    case {'pc1','pc2'}
        IX = ~isnan(sum(M));
        [~,score] = pca(M(:,IX));
        if isempty(score)
            disp('PCA failed')
            sort_val = 1:Rows(M);
            sort_val = sort_val(:);
        else
            sort_val = score(:,1);
            if strcmpi(sort_type,'pc2')
                sort_val = score(:,2);
            end
        end
    case 'pdist'
        sort_val = clusterdata(M,1);
    case 'kmeans'
        [~,sort_val2] = max(M,[],2); % Have the sort within the cluster be the peak.
        if nargin < 3
            param = 3;
        end
        % pca and kmeans can't handle NANs
        MM = M;
        MM(isnan(MM)) = 0;
        if Cols(M) > Rows(M)
            [~,sc] = pca(MM);
            sort_val = kmeans(sc(:,1:min(round([Cols(M)/2 Rows(M)/2]))),param);
        else
            sort_val = kmeans(MM,param);
        end
    case 'hierarchical'
        MM = M;
        MM(isnan(MM)) = 0;
        klist=2:7;%the number of clusters you want to try
        eva = evalclusters(MM,'linkage','CalinskiHarabasz','klist',klist); % there are other measures like silhouette
        sort_val = clusterdata(MM,'Linkage','ward','SaveMemory','off','Maxclust',eva.OptimalK);
        
        
end
[v,six] = sortrows([sort_val'; 1:Rows(M);sort_val2']',[1 3]);
OUT = M(six,:);
