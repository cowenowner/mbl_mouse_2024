function h= Plot_Q_cowen(t, Q, title_str, x_label, y_label, colormap_name)
% 
if nargin < 3
    title_str = '';
end
if nargin < 4
    x_label = '';
end
if nargin < 5
    y_label = '';
end
if nargin < 6
    colormap_name = 'turbo';
end

figure
h(1) = subplot(2,1,1);
imagesc(t,[],Q)
plot_vert_line_at_zero
pubify_figure_axis
clim(prctile(Q(:),[2.5 97.5]))
% ylabel('Neuron')
sgtitle(title_str)
colormap(colormap_name)
colorbar_label(y_label)

h(2) = subplot(2,1,2);
plot_confidence_intervals(t,Q)
xlabel(x_label); ylabel(y_label)
plot_vert_line_at_zero
pubify_figure_axis

set(gcf,'Position',[112         214        1072         678])
