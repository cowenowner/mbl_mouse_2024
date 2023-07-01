% One approach is to generate a random train of ISI that has a desired LV is to use optimization
% We use optimization to change the intervals in such a way to optimize the LV
% measure. There are probably better approaches, like explicity changing
% the distribution and drawing from it - say from exponential to a point
% gaussian. The issue is that the mean neads to stay the same.
%
%%
mean_freq = 20;
n_pulses = 25;
loc_var_tgt = 2;
max_iter = 2360;
% Starting point.
rnd_var = .4*randn(n_pulses-1,1)*1/mean_freq;

% STIM = exprnd(1/mean_freq,n_pulses-1,1);
% STIM = repmat(1/mean_freq,n_pulses-1,1); % Starting with equal spacing seems best.
STIM = repmat(1/mean_freq,n_pulses-1,1) + rnd_var; % Starting with equal spacing seems best.
STIM = normrnd(1/mean_freq,1/mean_freq,n_pulses-1,1) + rnd_var; % Starting with equal spacing seems best.
STIM = abs(STIM);
v = LocalVariance(STIM)
% optimize
fun = @(x)DA_LV_opt(x, loc_var_tgt); % the function to optimze (in this directory) - it determines how much dopamine is released following each stim sequence and returns -1* this.
options = optimset('OutputFcn',{@DA_outfun_lv},'MaxIter', max_iter); % DA_outfun does the nice plot.
% options = optimset('Display','iter','PlotFcns',@optimplotx,'TolX',1e-2,'TolFun',1e-4,'MaxIter',1e3); % optimplotx is awesome - spits out the actaul data.
% options = optimset('Display','iter','PlotFcns',@optimplotx,'MaxIter', 500); % optimplotx is awesome - spits out the actaul data.
figure
[x,fval,exitflag,output] = fminsearch(fun,STIM,options); % note - with a custom outfun you can record the history. You can then tell
title(['Target: ' num2str(loc_var_tgt)])