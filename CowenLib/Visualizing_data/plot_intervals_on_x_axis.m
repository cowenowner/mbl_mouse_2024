function h = plot_intervals_on_x_axis(sted, labels, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h = plot_markers_simple(offsets, axis_to_use, linewidth, colors,linestyle))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth = 4;
colormap_to_use = lines(20);
a = axis;
min_y = a(3);
min_y_text = min_y*1.2;

Extract_varargin

if isempty(labels)
    labels = cell(Rows(sted),1);
end


for ii = 1:Rows(sted)
    if sted(ii,2) < sted(ii,1)
        sted(end+1,:) = [a(1) sted(ii,2)];
        labels{end+1} = labels{ii};
        sted(ii,:) = [sted(ii,1) a(2)];
        colormap_to_use(Rows(sted),:) = colormap_to_use(ii,:);
    end
end
hold on
for ii = 1:Rows(sted)
    plot(sted(ii,1:2), min_y([1 1]),'Color',colormap_to_use(ii,:),'LineWidth',linewidth)
    if ~isempty(labels{ii})
        text(sted(ii,1),min_y_text(1),labels{ii})
    end
end

