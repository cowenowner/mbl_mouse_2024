function DANA_visualize_Events(EVT)
fn = fieldnames(EVT);
clrs = lines(length(fn));
figure
for ii = 1:length(fn)
    t = EVT.(fn{ii}).time_sec;
    if ~isempty(t)
        plot(t, ii*ones(size(t)),'o','Color',clrs(ii,:))
        hold on
        str = sprintf('%s %1.3f s',fn{ii},1/median(diff(t)));
        text(0,ii,str,'FontSize',7)
    end
end
xlabel('sec')
% legend(fn)
