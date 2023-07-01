function stop = DA_outfun(x, optimValues, state)

stop = false;
subplot(1,5,1:2)
bar(x)
ylabel('sec')
xlabel('ISI ID')
title('Inter-stim intervals')
box off

subplot(1,5,3:4)
hold on;
t_sec = cumsum(x);
% plot(t_sec,ones(size(t_sec))*optimValues.iteration,'.')
clrs = lines(length(t_sec));
for jj = 1:length(t_sec)
    plot(t_sec(jj),ones(size(t_sec))*optimValues.iteration,'.','Color',clrs(jj,:))
    hold on
end

ylabel('iter')
xlabel('sec')
title('Stim sequence')
% The following adds a lot of time.

subplot(1,5,5)
hold on
plot(1-optimValues.fval,optimValues.iteration, 'k.')
title('[DA]')
xlabel('[DA]')
ylabel('iter')
drawnow
