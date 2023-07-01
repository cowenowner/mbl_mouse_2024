function Connectivity_plot(C,labels)
%function Connectivity_plot(C,labels)
% INPUT: The C structure from Connectivity_matrix
%        labels for each region (e.g. {'CA1', 'CA2'}
% OUTPUT: a plot.
%% Cowen
colors = {'b' 'g' 'r' 'k' 'c' 'm' 'b' 'g' 'r' 'k' 'c' 'm' 'b' 'g' 'r' 'k' 'c' 'm' 'b' 'g' 'r' 'k' 'c' 'm' 'b' 'g' 'r' 'k' 'c' 'm' };
clf
imagesc(C.W); % or spy
a_color = 'y';
box off
axis xy
hold on
[r,c] = size(C.W);
for iRegion = 1:C.nRegions
    h = rectangle('position',[1,C.Nix_ex{iRegion}(1), c,C.nNeurons(iRegion)-2],...
        'LineWidth',2,'LineStyle','-');
    set(h,'EdgeColor',colors{iRegion})
    h = rectangle('position',[C.Nix_ex{iRegion}(1),1, C.nNeurons(iRegion)-2,c],...
        'LineWidth',2,'LineStyle','-');
    set(h,'EdgeColor',colors{iRegion})
    midpos = C.Nix_ex{iRegion}(1) + C.nNeurons(iRegion)/2;
    text(1,midpos,labels{iRegion},'Color',a_color)
    text(sum(C.nNeurons)*.96,midpos,labels{iRegion},'Color',a_color)
    text(midpos,c*.97,labels{iRegion},'Color',a_color)
end
ylabel('Neuron Source')
xlabel('Neuron Target')
title('Connectivity')