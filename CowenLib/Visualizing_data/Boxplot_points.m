function Boxplot_points(V,G, varargin)
% PLot a boxplot and overlay points. In format of data (V) and categorical
% group information (G). 
% Cowen 2019
MarkerSize = 8;
MarkerFaceAlpha = 0.6;
LineWidth = 1.5;

Extract_varargin;
if nargin == 1 || isempty(G)
    D = [];
    for iC = 1:Cols(V)
        D{iC} = V(:,iC);
    end
    [V,G] = group_data(D);
end



hold on
u = unique(G);
for ii = 1:length(u)
    IX = G == u(ii);
    VV = V(IX);
    x = (rand(length(VV),1)-0.5)/6;
    % scatter(x + ii, VV,MarkerSize,'k','filled','MarkerFaceAlpha',MarkerFaceAlpha); % Alpha alwas screws stuff up in ppt.
    plot(x + ii, VV,'.','Color',[.5 .5 .5],'MarkerSize',7); 
end
% pubify_figure_axis
if length(u) > 2
    p = anova1(V,G,'off');
    str = 'ANOVA';
else
    p = ranksum(V(G == u(1)),V(G == u(2)));
    % [~,pp] = ttest2(Firing_Rate_Hz(G =='pDA'),Firing_Rate_Hz(G =='Narrow'))
    % [~,pp] = ttest2(log10(Firing_Rate_Hz(G =='pDA')),log10(Firing_Rate_Hz(G =='Narrow')))
    str = 'rank';
end

boxplot(V,G,'Notch','on','Whisker',1,'Symbol','')
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
lines = findobj(gcf, 'type', 'line');
set(lines, 'LineWidth', LineWidth);

title(sprintf('%s p=%0.8f',str, p),'FontSize',8 )
