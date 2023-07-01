function stop = DA_outfun_lv(x, optimValues, state)
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
plot(t_sec,ones(size(t_sec))*optimValues.iteration,'b.')
ylabel('iter')
xlabel('sec')
title('Stim sequence')
% The following adds a lot of time.

lv = LocalVariance(x);
subplot(1,5,5)
hold on
plot(lv,ones(size(t_sec))*optimValues.iteration, 'k.')
title('LV')
xlabel('LV')
ylabel('iter')
drawnow
