% Generate ISIs given a continuum between fixed to exponential
% distributions.
%
%%
% clearvars
close all
wts = [0 1 0];
mean_freq = 20;
min_isi_s = 2/1000; % 2ms min interval 
n_sec_stim = 1;
pw = 2.8; % 1.8 is approximately LV = 1.6.
n_trains = 30;
edges_s = 0:.005:.2;
t_sec = []; lv = []; C = [];
for iT = 1:n_trains
     wts = rand(1,3);
     pw = pw + randn(1)*.2;
    % once and a while, just make some fixed freq like things
    if rand(1,1) < .1
       wts = [1 rand(1,2)];
    end
    if rand(1,1) < .02
       wts = [1 0 0];
    end
    [t_sec{iT},lv(iT)] = Generate_spikes_with_different_burst_stats(wts, pw, mean_freq, n_sec_stim,min_isi_s);
    d = diff(t_sec{iT});
    C(iT,:) = histcounts(d,edges_s);
    cv(iT) = std(d)/mean(d);
end
%%%%%%%%%%%%%%%%%%
[~,six] = sort(lv);
% [~,six] = sort(cv);
% Resort the t_sec to be in the order of burstiness.
tmp = [];
for ii = 1:length(t_sec)
    tmp{ii} = t_sec{six(ii)};
end
t_sec = tmp;
lv = lv(six);
cv = cv(six);
C = C(six,:);

figure
subplot(1,3,1)
imagesc(edges_s,[],C)
title('ISI histogram')
axis xy
pubify_figure_axis

figure
subplot(1,3,1:2);
plot_raster(t_sec)
axis tight
title('Pulse trains')
pubify_figure_axis
ylabel('Stim Sequence')
xlabel('sec')
subplot(1,3,3);
plot(lv,1:length(lv),'LineWidth',2)
pubify_figure_axis
xlabel('LV')
title('Burst Measure')
set(gca,'YTickLabel',[])
axis tight

figure
mid = round (length(t_sec)/2);
plot_raster(t_sec{1},1)
hold on
plot_raster(t_sec{mid},2)
plot_raster(t_sec{end},3)
axis tight
title(sprintf('Pulses %1.2f %1.2f %1.2f ',lv(1),lv(mid),lv(end)))
pubify_figure_axis
ylabel('Stim Sequence')
xlabel('sec')

save('spike_stats','t_sec','lv')
