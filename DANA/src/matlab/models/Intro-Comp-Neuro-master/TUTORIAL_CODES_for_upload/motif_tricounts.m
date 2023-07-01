function [ uniP, biP, Nconns, Nbi, Ntriples, N3expect ] = motif_tricounts( conn_matrix)
% motif_tricounts counts numbers of bidirectional and triplet connections 
% vs chance
%
%   [uniP, biP, Nconns, Nbi, Ntriples, N3expect ] = motif_tricounts(conn_matrix) 
%
%   The function requires as inputs:
%   conn_matrix is the input square connectivity matrix
%
%   It returns the following set of statistics:
%   uniP is mean connection probability (i.e. fraction of connections that
%   are present)
%
%   biP is empirical probability of a bidirectional connection (i.e. the
%   fraction of cell-pairs connected bidirectionally)
%
%   Nconns is total number of connections
%
%   Nbi is double the total number of bidirectionally connected pairs (i.e.
%   the number of connections 
%
%   This code is required for use with Tutorial 5.2 in Chapter 5 of the
%   textbook
%   An Introductory Course in Computational Neuroscience 
%   by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N1,N2] = size(conn_matrix);    % N1 and N2 should both equal no. of neurons
if ( N1 ~= N2 )                 % If they are not the same (i.e. not a square matrix)
    disp(['Error! N1', N1, 'N2', N2])   % Display an error message
end

 % Total number of connections is sum of non-zero entries in the
 % connectivity matrix
Nconns = sum(sum(conn_matrix>0));  
% Estimate of connection probability is no. of connections divided by number of 
% potential connections (including self-connections) 
uniP = Nconns/(N1*N2);  

Nself = sum(diag(conn_matrix>0));   % No. of self-connections on the diagonal
selfP = Nself/N1;                   % Probability of a self-connection

%% Now count number of bidirectionally connected pairs
Nbi = 0;                            % Initialize counter at zero
for i = 1:N1;                       % Loop through all potential neurons         
    for j = 1:i-1;                  % Loop through all potential partners not 
        if ( conn_matrix(i,j) )     % If there is a connection pre to post
           if ( conn_matrix(j,i) )  % and if there is a connection post to pre
               Nbi = Nbi + 1;       % Accumulate a bidirectional correction
           end
        end
    end
end

biP = 2*Nbi/Nconns;     % Probability of a connection being in a bidirectionally connected pair                

Npairs = N1*(N1-1)/2;   % Number of potential pairs for bidirectional connections
biP2 = Nbi/Npairs;      % Fraction of pairs with bidirectional connections
uniP2 = (Nconns-Nself-2*Nbi)/Npairs;    % Fraction of pairs with unidirectional connection
zeroP2 = 1-biP2-uniP2;                  % Fraction of pairs with no connection

% The next line converts the probability that a pair has only one
% connection (uniP2) to the probability that a connection is present in a
% pair while the other connection is absent. The factor of two arises
% because there are two ways of producing each singly connected pair
uniP_conn = uniP2/2;
%% Now calculate the expected numbers of each motif for the connection 
% patterns between three different neurons, comprising 3 pairs where each
% pair can be unconnected, singly connected or bidirectionally connected
% For some combinations, the order and direction of singly connected pairs
% matters. Combining the ordering and directionality, with the combinatorial 
% ways of arranging the group of 3 cells leads to different multiplicity 
% factors. 
% 
N3expect = zeros(1,16);                    
N3expect(1) = 1*zeroP2*zeroP2*zeroP2;           % No connections between any pair (multiplicity of 1)
N3expect(2) = 6*zeroP2*zeroP2*uniP_conn;        % One unidirectional connection (multiplicity of 6)
N3expect(3) = 3*zeroP2*zeroP2*biP2;             % One bidirectional connected pair (multiplicity of 3)
N3expect(4) = 3*zeroP2*uniP_conn*uniP_conn;     % 2 unidirectional connections (multiplicity of 3)
N3expect(5) = 3*zeroP2*uniP_conn*uniP_conn;     % 2 unidirectional connections (multiplicity of 3)
N3expect(6) = 6*zeroP2*uniP_conn*uniP_conn;     % 2 unidirectional connections (multiplicity of 6)
N3expect(7) = 6*zeroP2*uniP_conn*biP2;          % 1 uni- and 1 bidirectionally connected pair (multiplicity of 6)
N3expect(8) = 6*zeroP2*uniP_conn*biP2;          % 1 uni- and 1 bidirectionally connected pair (multiplicity of 6)
N3expect(9) = 3*zeroP2*biP2*biP2;               % 2 bidirectionally connected pairs (multiplictiy of 3)
N3expect(10) = 6*uniP_conn*uniP_conn*uniP_conn; % 3 unidirectional connections (multiplicity of 6)
N3expect(11) = 2*uniP_conn*uniP_conn*uniP_conn; % 3 unidirectional connections (multiplicity of 2)
N3expect(12) = 3*uniP_conn*uniP_conn*biP2;      % 2 uni- and 1 bidirectionally connected pairs (multiplicity of 3)
N3expect(13) = 6*uniP_conn*uniP_conn*biP2;      % 2 uni- and 1 bidirectionally connected pairs (multiplicity of 6_
N3expect(14) = 3*uniP_conn*uniP_conn*biP2;      % 2 uni- and 1 bidirectionally connected pairs (multiplicity of 3)
N3expect(15) = 6*uniP_conn*biP2*biP2;           % 1 uni- and 2 bidirectionally connected pairs (multiplicity of 6)
N3expect(16) = 1*biP2*biP2*biP2;                % All 3 pairs bidirectionally connected (multiplicity of 1)

N3expect=N3expect*N1*(N1-1)*(N1-2)/6;           % Multiply probabilties by total number of triplets


%% Now count the actual number of different motifs for triplets of cells
Ntriples = zeros(1,16);                 % Initialize as zero

% The variable 'type' will indicate the types of connections between the 3
% pairs of cells in a triplet i,j,k. Each pair (eg i,j) has 4
% possibilities, indicated with a number 0-3 in an element of 'type'.
type = zeros(1,3,'int8');           

for i = 1:N1;                           % First cell labeled i
    
    for j = i+1:N1;                     % Second cell labeled j
        type(1) = conn_matrix(i,j) + 2*conn_matrix(j,i);    % A number 0-3
        for k = j+1:N1
            type(2) = conn_matrix(j,k) + 2*conn_matrix(k,j);    % A number 0-3
            type(3) = conn_matrix(k,i) + 2*conn_matrix(i,k);    % A number 0-3
%            
            num = sum(type > 0 );       % Number of connected pairs
            num2 = sum( type == 3);     % Number of bidirectionally connected pairs
            if ( num == 0 ) 
                motif = 1;              % No connections in the triplet
            else 
                if ( num == 1 )         % One connected pair
                    if ( num2 == 0 )
                        motif = 2;      % One unidirectional connection 
                    else
                        motif = 3;      % One bidirectional connection
                    end
                else
                    if ( num == 2 )         % Two connected pairs
                        if ( num2 == 2 ) 
                            motif = 9;      % Both bidirectionally connected
                        else
                            c2 = find(type>0);      % Finds which pairs are connected
                            if ( num2 == 1 )        % One bidirectionally connected pair
                                if ( sum(type) == 4 )  % must be a "3" and a "1" 
                                    type2 = min(type,2);
                                    if ( mod(type2(2)-type2(1),3 ) == 1 )
                                        motif = 7;
                                    else
                                        motif = 8;
                                    end
                                else % sum(type) = 5, must be a "3" and a "2"
                                    type2 = max(type,1);
                                    if ( mod(type2(1)-type2(2),3 ) == 1 )
                                        motif = 7;  % Unidirectional connection to one of bidirectionally connected cells
                                    else
                                        motif = 8;  % Unidirectional connection away from one of bidirectionally connected cells
                                    end
                                end
                            else
                                if( type(c2(1)) == type(c2(2)) ) 
                                    motif = 6;  % Two single connections in same direction around triangle
                                else
                                    if ( mod(type(2)-type(1),3) == 1 ) 
                                        motif = 5;  % Two unidirectional connections with same endpoint
                                    else
                                        motif = 4;  % Two unidirectional connections with same start point
                                    end
                                end
                            end
                        end
                    else  % num = 3, so all 3 pairs of cells have connections
                        if ( num2 == 3 ) 
                            motif = 16;         % all cells bidirectionally connected
                        else
                            if ( num2 == 2 ) 
                                motif = 15;     % 2 bidirectionally and 1 unidirectionally connected pairs
                            else
                                if ( num2 == 1 )    % 1 bidirectionally connected and 2 unidirectionally connected pairs
                                    if ( sum(type) == 6 )
                                       if ( mod(type(2)-type(1),3) == 1 )
                                           motif = 14;  % Both unidirectional connections point away from bidirectionally connected pair
                                       else
                                           motif = 12;  % Both unidirectional conenctions point toward the bidirectionally connected pair
                                       end                                        
                                    else
                                        motif = 13; % One unidirectional points away, the other toward the bidirectionally connected pair
                                    end
                                else % num2 = 0 so all three connections are unidirectional
                                    if ( mod(type(1)+type(2)+type(3),3) == 0 )
                                        motif = 11; % Connections point the same direction around the triangle
                                    else
                                        motif = 10; % One connection points in the reverse direction of the other two
                                    end
                                end
                            end
                        end
                    end
                end
            end    
            Ntriples(motif) = Ntriples(motif) + 1;  % Accumulate 1 to the total number of the motif now determined
        end
    end
end

end

