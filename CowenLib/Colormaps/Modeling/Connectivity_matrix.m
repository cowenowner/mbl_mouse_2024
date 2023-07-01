function C = Connectivity_matrix(nTotalNeurons,pNeurons,pExitatory,CP,dist_type)
% Creates a connectivity matrix for the neurons given the probability of
% connection for each region (e.g. CA3, CA1...)
% INPUT:
%   nTotalNeurons = total number of neurons (scalar)
%   pNeurons = proportion of neurons in each region (sums to 1). One
%     element in the vector for each region. (vector) nRegions =
%     length(pNeurons)
%   pExcitatory = proportion of neurons in each region that are excitatory.
%     (vector of nRegions length).
%   CP = Probability of connectivity for each region (4 x nRegions x
%     nRegions). The first index indicates the pair
%     (1=exex,2=inin,3=exin,4=inex), the second 2 indices indicate the
%     region pairs. For example CP(2,3,1) indicates the prob of connection
%     of pairs of excitatory cells from region 3 to region 2. This value is
%     then used in a random number generator to create the connectivity
%     matrix.
%   CW = The weights of the connection. Same form as CP, except this
%     describes the distribution of the stregth of connection weights between
%     the region-pairs. For instance, 2 regions may be connected quite
%     densely, but the strength of the connections between these regions may
%     be somewhat weak - e.g. they connect to the distal dendrites. Other
%     regions or pairs e.g. inhibitory neurons to pyr cells within a region,
%     may be connected somewhat sparsely, but have very strong connection
%     weights.
%   dist_type = probability distribution (e.g. small world, gaussian,log).
%   dist_params = parameters for each distribution (e.g. mean and sd).
% OUTPUT:
%   C = a structure that describes the connectivity between regions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cowen (2008).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_type_CP = 'exp';
dist_type_CW = 'unif';
wt_scale_factor = 0.5;
nRegions = length(pNeurons);
C.W = sparse(nTotalNeurons,nTotalNeurons); % Connectivity
last_ix = 0;
thresh = 0.15; % Threshold for having a synapse.
C.All_Nix_ex  = [];
C.All_Nix_inh = [];
for iRegion = 1:nRegions
    % Number of neurons in each region (all and inh and exc)
    nNeurons(iRegion) = pNeurons(iRegion)*nTotalNeurons;
    nExc(iRegion) = round(nNeurons(iRegion) * pExitatory(iRegion));
    nInh(iRegion) = nNeurons(iRegion) - nExc(iRegion); %#ok<NASGU>
    % Indices to the neurons in each region.
    Nix{iRegion}     = last_ix + 1:round( nNeurons(iRegion) + last_ix);
    Nix_ex{iRegion}  = Nix{iRegion}(1:nExc(iRegion));
    Nix_inh{iRegion} = Nix{iRegion}((nExc(iRegion)+1):end);
    % Keep and index of all of the inh and ex neurons.
    C.All_Nix_ex = [C.All_Nix_ex Nix_ex{iRegion}];
    C.All_Nix_inh = [C.All_Nix_inh Nix_inh{iRegion}];
    %
    last_ix = Nix{iRegion}(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connect the neurons. This is teh probability of 2 neurons being
% connected AT ALL. The next part assignes the connections a weight.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iRegion1 = 1:nRegions
    for iRegion2 = 1:nRegions
        nre = length(Nix_ex{iRegion1});
        nce = length(Nix_ex{iRegion2});
        nri = length(Nix_inh{iRegion1});
        nci = length(Nix_inh{iRegion2});
        % Index to cell pairs with synapses
        Six = find(random(dist_type_CP,CP(1,iRegion1,iRegion2),[nre nce]) > thresh);
        % Weights of the synapses.
        if ~isempty(Six)
            W = sparse(nre,nce);
            W(Six) = random(dist_type_CW,CP(1,iRegion1,iRegion2),1,size(Six));
            C.W(Nix_ex{iRegion1},Nix_ex{iRegion2})   = W*wt_scale_factor;
        end
        % Index to cell pairs with synapses
        Six = find(random(dist_type_CP,CP(2,iRegion1,iRegion2),[nri nci]) > thresh);
        % Weights of the synapses.
        if ~isempty(Six)
            W = sparse(nri,nci);
            W(Six) = random(dist_type_CW,CP(2,iRegion1,iRegion2),1,size(Six));
            C.W(Nix_inh{iRegion1},Nix_inh{iRegion2}) = W* wt_scale_factor;
        end
        % Index to cell pairs with synapses
        Six = find(random(dist_type_CP,CP(3,iRegion1,iRegion2),[nre nci]) > thresh);
        % Weights of the synapses.
        if ~isempty(Six)
            W = sparse(nre,nci);
            W(Six) = random(dist_type_CW,CP(3,iRegion1,iRegion2),1,size(Six));
            C.W(Nix_ex{iRegion1},Nix_inh{iRegion2})  = W* wt_scale_factor;
        end
        % Index to cell pairs with synapses
        Six = find(random(dist_type_CP,CP(4,iRegion1,iRegion2),[nri nce]) > thresh);
        % Weights of the synapses.
        if ~isempty(Six)
            W = sparse(nri,nce);
            W(Six) = random(dist_type_CW,CP(4,iRegion1,iRegion2),1,size(Six));
            C.W(Nix_inh{iRegion1},Nix_ex{iRegion2})  = W* wt_scale_factor;
        end
    end
end
C.Nix_ex = Nix_ex;
C.Nix_inh = Nix_inh;
C.Nix = Nix;
C.nNeurons = nNeurons;
C.nRegions = nRegions;
C.nEx = nExc;
C.nInh = nInh;

