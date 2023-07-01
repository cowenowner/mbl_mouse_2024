% Tutorial_5_2.m
%
% A code to produce a connectivity matrix with hidden structure, then
% extract that structure from the connectivity matrix and visualize it by
% reordering the connections. 
%
% The code also calls the function motif_tri_counts to extract the number
% of bidirectionally connected cell-pairs and the numbers of each distinct
% motif of connections within triplets of cells.
%
% This code is a solution to tutorial 5.2 of the textbook:
%
% "An Introductory Course in Computational Neuroscience" 
%
% by Paul Miller, Brandeis University.
%
% March 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                              
close all;

%% Set up initial connectivity matrix
Ngroups = 9;                    % No. of distinct groups or assemblies
Ncells = 200;                    % Total number of cells to consider

pconn_in_group = 0.5;           % Within-group connection probability
pconn_out_group = 0.02;         % Between-group connection probability
pconn_diff = pconn_in_group-pconn_out_group;    % Difference within-between

group_vals = ceil(rand(1,Ncells)*Ngroups);      % Randomly assign group IDs
conns = zeros(Ncells);                          % Set up connectivity matrix
for i = 1:Ncells;                               % Loop through each row
    % The following line produces a row of connections, each with 
    % probability that depends on whether the group_vals label of the 
    % column is the same as the group_vals label of the row
    conns(i,:) = rand(1,Ncells) < (pconn_out_group ...
        + pconn_diff*(group_vals == group_vals(i)) );
end

mean_pconn = sum(sum(conns))/(Ncells*(Ncells-1));   %  mean connection probability

% Now plot the initial randomly ordered connectivity matrix
imagesc(conns)
title('Unordered connectivity')
xlabel('Postsynaptic Cell No.')
ylabel('Presynaptic Cell No.')

%% Now group according to the prescribed groups of cells
% by sorting the numbers, all the group_vals entries with the same number
% will be placed successively in the vector "neworder1"
[groupnumber, neworder1] = sort(group_vals);

% A new matrix with rows and columns in the new order
newconns = conns(neworder1,neworder1);      
figure(2)
imagesc(newconns);              % View the connectivity ordered by group
title('Connectivity ordered by design')
xlabel('Postsynaptic Cell No.')
ylabel('Presynaptic Cell No.')

%% Now, without knowing group labels, attempt to reorder the connectivity
%  matrix to obtain the clustered structure.
%  In this method, we calculate the correlation between the set of
%  connected partners of one cell with the set of connected partners of
%  another cell. If two cells share many connected partners they are likely
%  to belong to the same group.
%  Since two cells in the same group are likely to be connected to each 
%  other, we insert diagonal elements into the connectivity matrix before 
%  calculating the correlation.

corr_conns1 = corr(conns);         % Column to column correlations
corr_conns2 = corr(conns');        % Row to row correlations
total_corr = corr_conns1+corr_conns2;   % Sum of column and row correlations

figure(4)
hist(total_corr(:))                     % Plot a histogram of these correlations
title('Distribution of connectivity correlations')
xlabel('Correlation of sets of connections')
ylabel('No. of cell-pairs')

% Select a threshold for same-group membership, based partly on the
% distribution observed in the histogram. This can be adjusted. A higher
% threshold results in more distinct groups being extracted.
corr_threshold = 0.5;                   

groupnum = 0; % This will count through the groups. Initialize it to zero                      
group = zeros(1,Ncells);    % This will label the group of each cell

for i = 1:Ncells;                           % Loop through all the cells
    if ( group(i) == 0 );                   % If cell i not already assigned a group
        groupnum = groupnum + 1;            % It must be in a new group
        cells_remain_list = i;              % List of cells in group to find partners
        cells_used_list = [];               % List of cells in group with partners found
        
        % The while loop will keep adding cells to the group with an
        % above-threshold correlation with a cell already in the group. 
        % All cells added to the group must be considered--and these cells
        % enter the "cells_remain_list".
        % Once they are considered -- and all their in-group partners found
        % -- they join the cells_used_list
        while (cells_remain_list)               % while cells remain in the list
            cell = cells_remain_list(1);        % keep taking the first cell
            % Find all cells with above-threshold correlation with the cell
            % under consideration
            groupcells = find(total_corr(cell,:)>corr_threshold);
            % Assign the same group label to all cells found this way 
            % (This will include the cell under consideration, since its
            % correlation with itself is maximum).
            group(groupcells) = groupnum;
            
            % Add the cell to the list of used cells
            cells_used_list = [cells_used_list, cell];
            
            % Add the ones just found to the list of remaining cells to
            % consider
            cells_remain_list = [cells_remain_list, groupcells];
            
            % Now remove cells fom cell_remain_list that are repeated
            % The list without repetitions is "new_cells_remain_list"
            new_cells_remain_list = [];
            
            % Loop through all cells in the "remain" list
            for k = 1:length(cells_remain_list)
                testcell = cells_remain_list(k);    % ID of a given cell
                
                % Then find if that cell is already in the 
                % "new_cells_remain_list"
                repeat = length(find(new_cells_remain_list == testcell));
                
                if ( ~repeat )  % If the cell is not present in the new list
                                % then add it to the new list
                    new_cells_remain_list = [new_cells_remain_list, testcell];
                end
            end
            
            % Now remove cells that are already on the cells_used_list
            cells_remain_list = [];     % To store cells not already used
            
            % Loop through all cells that potentially need checking
            for k = 1:length(new_cells_remain_list)
                testcell = new_cells_remain_list(k);    % ID of a given cell
                % Then find if this cell is in the list of used cells
                repeat = length(find(cells_used_list == testcell));
                if ( ~repeat )          % If it is not already used
                                        % add it to the list of remaining
                    cells_remain_list = [cells_remain_list, testcell];
                end
            end
            
        end;    % End of the "while loop" for remaining cells in this group
    end;        % End of "if" to test if the next cell needs a new group
    
end;   % End of "for loop" to go to next cell in the original list to look for new groups


%% Now group according to the extracted groups of cells
[groupnumber2, neworder2] = sort(group);

newconns2 = conns(neworder2,neworder2);
figure(3)
imagesc(newconns2)
title('Connectivity by extracted order')
xlabel('Postsynaptic Cell No.')
ylabel('Presynaptic Cell No.')
%% Finally, use the function "motif_tricounts" to generate the abundances 
% relative to chance of each pair and triplet motif.
[ uniP, biP, Nconns, Nbi, Ntriples, N3expect ] = motif_tricounts(conns);

% Plot the ratio of actual to expected abundance of each triplet motif, 
% where the expected abundance is calculated assuming any connection can
% exist with a constant probability calculated as no. of connections
% divided by number of possible connections.
figure(5)
plot(Ntriples./N3expect)
title('Motif counts versus expected no.')
xlabel('Motif No.')
ylabel('Ratio observed/random')






