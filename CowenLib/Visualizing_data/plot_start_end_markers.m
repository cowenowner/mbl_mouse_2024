function plot_start_end_markers(start_end_intervals,y)
if nargin < 2
    y = 1;
end

hold on
plot(start_end_intervals(:,1),y*ones(size(start_end_intervals(:,1))),'g>')
if Cols(start_end_intervals) > 1
    plot(start_end_intervals(:,2),y*ones(size(start_end_intervals(:,2))),'r<')
end