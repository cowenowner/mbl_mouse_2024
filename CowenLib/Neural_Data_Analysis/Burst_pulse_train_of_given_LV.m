function [t] = Burst_pulse_train_of_given_LV(mean_freq, loc_var_tgt, n_pulses, duration_sec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: For 1) a given LV, 2) a desired number of repeating spikes (say 3
% spikes - so two inter-spike-intevals), and 3) a target frequency, generate a train for a given target
% local variance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% clearvars
% [t] = Burst_pulse_train_of_given_LV(10, 0:.2:1.4, 20, 10)
PLOT_IT = false;
if nargin == 0
    mean_freq = 20;
    loc_var_tgt = 1.1;
    n_pulses = 3; % number of pulses. Must be 3 or more. 2 will just give you a tonic firing.
    duration_sec = 10;
end
if length(loc_var_tgt) > 1
    t = cell(length(loc_var_tgt),1);
    for ii = 1:length(loc_var_tgt)
        t{ii} = Burst_pulse_train_of_given_LV(mean_freq, loc_var_tgt(ii), n_pulses, duration_sec);
    end
    if 1
        figure
        subplot(1,3,1:2)
        plot_raster(t)
        subplot(1,3,3)
        barh(loc_var_tgt)
        xlabel('LV');axis tight; grid on
    end
    return
end

max_iter = 2360;
% Starting point.
rnd_var = .4*randn(n_pulses-1,1)*1/mean_freq;
% optimize

STIM = normrnd(1/mean_freq,1/mean_freq,n_pulses-1,1) + rnd_var; % Starting with equal spacing seems best.
STIM = abs(STIM);

fun = @(x)LV_opt(x, loc_var_tgt); % the function to optimze (in this directory) - it determines how much dopamine is released following each stim sequence and returns -1* this.
if PLOT_IT
    options = optimset('OutputFcn',{@Outfun_lv},'MaxIter', max_iter); % DA_outfun does the nice plot.
    figure(1010)
else
    options = optimset('MaxIter', max_iter); % DA_outfun does the nice plot.
end
[x,fval,exitflag,output] = fminsearch(fun,STIM,options); % note - with a custom outfun you can record the history. You can then tell
% make sure the sum of x is equal to our target frequency.
% add some time so that the interval lasts
% target_dur_sec = 1;
% x(1) = x(1)+target_dur_sec-sum(x);
%   x(end+1) = target_dur_sec-sum(x);
if PLOT_IT
    title(['Target: ' num2str(loc_var_tgt)])
end
isis = repmat(x,ceil(duration_sec/sum(x))+5,1);
% rescale to the proper sampling rate.
v = (1/mean_freq)/mean(isis);
isis = isis*v;
t = cumsum(isis,1);
% n_spikes = len

t = t(t<= duration_sec);
% length(t)/(t(end)-0)
% Optimize such that the trains do not overlap with the scan pulses.

% figure
% plot(t,ones(size(t)),'o')
end

function stop = Outfun_lv(x, optimValues, state)
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
end
function [objective]  = LV_opt(ISI, lv_tgt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize so that local variance is varied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isi_thresh = 0.003; % could have this as an input.
v = (ISI(1:(end-1))-ISI(2:end))./(ISI(1:(end-1))+ISI(2:end));
LV = 3*sum(v.^2)/(length(ISI)-1);
% penalize any negative ISIs or ISIs shorter than some interval.
negpen = double(sum(ISI<isi_thresh))*1.9;
objective = (LV-lv_tgt).^2 + negpen;
end