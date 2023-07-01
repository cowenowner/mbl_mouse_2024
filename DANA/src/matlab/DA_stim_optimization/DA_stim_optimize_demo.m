% A demonstration of how we could, in real time, optimize a sequence
% of pulses to trigger maximal dopamine release.
%
% Assumption: A specific sequence of stim pulses will be optimal for
% maximizing dopamine release.
%
% Approach: Create a sequence of target stim sequences (TGT) stims defined as an ordered sequence of
% inter-stim intervals of length n_pulses. Dopamine release will be maximal
% when the actual stim (STIM) is closest (measured as MSE or some other
% distance measure). This routine searches the high-dimensional inter-stim
% interval space to find the optimal sequences that evoke maximal dopamine
% release. It uses the non-linear fminsearch method.
%
% We could assume that there is only one optimal solution, but this may be
% a very sketchy assumption. As a result, mulitple TGT vectors (rows) can
% be included to see how fminsearch works when there are multiple minima.
%
% ISSUES: sometimes negative ISIs come out of optimization. How do we adapt
% to this? Just set negatives to zero during optimization? 
% May not converge when mulutiple pattenrs - most likely given the max
% operation in the optimization function. A batch optimization procedure
% might work better - run a couple random iterations.
% Another constraint could also limit the ISIs during optimization.
% 

%% Cowen 2020

% Each element of TGT is the delay, in sec, between adjacten stims (sec).
close all
mean_freq = 40; % mean stim frequency.
n_pulses = 10; % number of pulses allowed in each block of stim. Tried 40 pulses and never really converges. 
max_iter = 140;
rng(4); % For replication.
TGT = []; % the targets that will produce the most dopamine
% By having >1 TGT, you make the problem much harder, but probably more
% realistic. If it does fluctuate, then how do we get it to better converge
% on one solution instead of some half-assed middle solution?
%
% since rand averages the .5, that means that the last pulse will 
mx = .5*n_pulses;
TGT(1,:) = rand(1,n_pulses - 1) ; % Random sequences but make last one such that it adds to 1s.
if sum(TGT(1,:)) < mx
    TGT(1,end) = TGT(1,end) + mx-sum(TGT(1,:));
end
TGT(1,:) = TGT(1,:)/mean_freq;

% TGT(2,:) = rand(1,n_pulses - 1) /mean_freq; % Random sequences
% TGT(2,:) = sort(rand(1,n_pulses - 1) /mean_freq); % Random sequences fast to slow.

figure
plot(TGT')
title('ISIs that evoke max DA')
xlabel('ISI ID')
ylabel('sec')


% Starting point.
STIM = rand(1,n_pulses-1)*1/mean_freq;
% optimize
fun = @(x)DA_opt1(x, TGT); % the function to optimze (in this directory) - it determines how much dopamine is released following each stim sequence and returns -1* this.
options = optimset('OutputFcn',{@DA_outfun},'MaxIter', max_iter); % DA_outfun does the nice plot.
% options = optimset('Display','iter','PlotFcns',@optimplotx,'TolX',1e-2,'TolFun',1e-4,'MaxIter',1e3); % optimplotx is awesome - spits out the actaul data.
% options = optimset('Display','iter','PlotFcns',@optimplotx,'MaxIter', 500); % optimplotx is awesome - spits out the actaul data.
figure
[x,fval,exitflag,output] = fminsearch(fun,STIM,options); % note - with a custom outfun you can record the history. You can then tell
figure
clrs = lines(20);

for ii = 1:size(TGT,1)
    subplot(2,size(TGT,1),ii)
    stem(TGT(ii,:),'k')
    ylabel('sec')
    yyaxis right
    stem(x,'r')
    axis tight
    xlabel('ISI ID')
    title(['Pattern ' num2str(ii)])
    subplot(2,Rows(TGT),ii+size(TGT,1))
    t_sec = cumsum(TGT(ii,:));
    act_t_sec = [0 cumsum(x)];

    for jj = 1:length(t_sec)
    plot(t_sec(jj),ii,'.','Color',clrs(jj,:))
    hold on
    plot(act_t_sec(jj),ii,'o','Color',clrs(jj,:))

    end
    xlabel('sec')
    title('Stim sequence')
end
