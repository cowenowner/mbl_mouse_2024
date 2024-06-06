% I have struggled mightily with cross corrs since grad school so hopefully
% this will help those embarking on cross-corr adventures...
% 
% Cross corrs are wonderful but they can yield surprisingly unintuitive results and also are
% subject to disturbing corruptions due to non-stationarities in the data.
% The later 
%
% For a list of some of the corruptions of cross-corrs and the solutions
% (via the generation of null distributions via suffling and randomizing in
% clever ways) see... (these are for spike trains, but the general
% principles and confounds are mostly the same for cross corss between
% continuous signals).
%
% Brody, C.D., 1999. Correlations without synchrony. Neural Comput 11, 1537–1551. https://doi.org/10.1162/089976699300016133
% Palm, G., Aertsen, A.M., Gerstein, G.L., 1988. On the significance of correlations among neuronal spike trains. Biol Cybern 59, 1–11. https://doi.org/10.1007/BF00336885
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(10); % fixes the random seed so that we get the same results.
n_pts = 250; % for voltammetry, I imagine we're concatenating all [DA] vectors under a single condition or do we calculate the cross corr per trial and then average the cross corrs across trials within that condition (5 trials per condition)?
n_lags = 100; % size of the xcorr
n_shuffle = 1000; % for generating a null distribution
%% Lets do the cross corr between 2 random number sequences.
% Intuition (including mine) would suggest a flat cross corr with some
% jitter and the mean somewhere near the average of the two signals.
v1 = rand(1,n_pts); % make a random sequence
v2 = rand(1,n_pts); % make another random sequence

[xc_orig,lags] = xcorr(v1, v2, n_lags,'coeff'); % NOTE: the normalized flag does not fix it.
figure;
subplot(2,1,1); hold on; plot(v1);plot(v2)
title('Data going into the cross corr'); legend('v1','v2')
subplot(2,1,2);plot(lags,xc_orig);grid on; ylabel('r')
title('The CrossCorr (xcorr(v1,v2)))')
% What happened? This rediculous triangle! Definitely not flat.
% This is a problem. Why?

%% Could it be because we did not subtract the mean?
[xc_orig,lags] = xcorr(v1 - mean(v1), v2 - mean(v2), n_lags,'coeff');
figure;
subplot(2,1,1); hold on; plot(v1 - mean(v1));plot(v2 - mean(v2))
title('Subtract mean: Data going into the cross corr'); legend('v1-mn','v2-mn')
subplot(2,1,2);plot(lags,xc_orig);grid on; ylabel('r')
title('The CrossCorr (xcorr(v1,v2)))')
% Well, at least that solved this problem. Note to self: always subtract the mean before running xcorr. 
% There are many other probems though as indicated in Brody - more
% nefarious problems due to non-stationarities and other things. cross
% corrs essentially presume the mean in both v1 or v2 stays the same.

%% Before getting to this point, let's look at what happens when we cross
% corr a random signal with a vector of ones...
v1b = v1 - mean(v1); v2b =  ones(size(v1));
[xc_orig,lags] = xcorr(v1b,v2b, n_lags,'coeff');
figure;
subplot(2,1,1); hold on; plot(v1b);plot(v2b)
title('data with ones: Data going into the cross corr'); legend('v1-mn','ones only')
subplot(2,1,2);plot(lags,xc_orig);grid on; ylabel('r')
title('The CrossCorr (xcorr(v1,v2)))')
% WTF! That gives us a very assymetric xcorr. I honestly do not know why.
% Someone can chime in on this. It's interesting but clearly 100% artifact
% and would lead to considerable confusion and false attribution of a
% temporal relationship between these two signals.
% You also can't replace the ones(size(v1)) with zeros(size(v1)) as that
% will alwyas yield an xcorr of zoro.

%% So what can one do (beyond subtracting a mean).
% First, let's create an oscillating signal.
x = linspace(0,10*pi,n_pts);
v1 = sin(x) + randn(size(x))*.1;
v2 = sin(x) + randn(size(x))*.1;
v2 = circshift(v2,10); % this shifts the vector nicely forward. The lage in the xcorr should be this amount.

[xc_orig,lags] = xcorr(v1 - mean(v1), v2 - mean(v2), n_lags,'coeff');
figure;
subplot(2,1,1); hold on; plot(v1 - mean(v1));plot(v2 - mean(v2))
title('Data going into the cross corr'); legend('v1-mn','v2-mn')
subplot(2,1,2);plot(lags,xc_orig);grid on; ylabel('r')
title('Cross corr(v1,v2)')
% Sweet- this actually looks correct in my mind. the peak is offset by
% about 10 lags. This looks like an xcorr that I would expect.
%% Let's re-evaluate however what not subtracting the mean would do... I'll
% just add 10 to one of the inputs (this is a scalar - will not change the
% temporal relationship AT ALL)
[xc_orig,lags] = xcorr(10 + (v1 - mean(v1)), v2 - mean(v2), n_lags,'coeff');
figure;
subplot(2,1,1); hold on; plot(v1);plot(v2)
title('Data going into the cross corr'); legend('v1','v2')
subplot(2,1,2);plot(lags,xc_orig);grid on; ylabel('r')
title('The CrossCorr (xcorr(v1,v2)))')
% Yep - causes crazyness. Terrible. Strange cross-corr. Yes, there is still
% an offest of about 10 but there is this very strange assymetry.


%% Non-stationarities and other nastyness.
% Slow changes in one or both of the vectors also corrupt the xcorr and
% this will be clearly an issue with voltammetry and extracellular ephys.
% Both signals are VERY drifty. The general solution for this is to
% subtract a 'null' cross corr from the original cross corr where the null
% cross corr is generated from surrogate input vectors with the same
% general statistics as the original (for example, by reordering the timing
% of the signal randomly in one of the vectors) and then repeating this
% multiple times. This does not entirely solve the problem as we shoul
% ideally randomize at local intervals in the input vector so that
% non-stationarities are preserved so that they can be subtracted away.
% This sounds confusing, so let's do the simpler thing of just subtracting
% a null distribution generated by generating random input...
%%
x = linspace(0,10*pi,n_pts);
v1 = sin(x) + randn(size(x));
v2 = sin(x) + randn(size(x));
v2 = circshift(v2,10); % this shifts the vector nicely forward.
[xc_orig,lags] = xcorr(v1 - mean(v1), v2 - mean(v2), n_lags,'coeff');

shuffle_xc = nan(n_shuffle,2*n_lags+1);
for iShuff = 1:n_shuffle
    v1sh = v1(randperm(length(v1)));
    % We could shuffle v2 as well but that is probably not necessary.
    shuffle_xc(iShuff,:) = xcorr(v1sh - mean(v1sh), v2 - mean(v2), n_lags,'coeff');
end
figure
subplot(4,1,1); hold on; plot(v1);plot(v2)
title('Data going into the cross corr'); legend('v1','v2')
subplot(4,1,2);plot(lags,xc_orig);grid on; ylabel('r')
title('The CrossCorr (xcorr(v1,v2)))')
subplot(4,1,3);imagesc(lags,[],shuffle_xc)
title('cross corrs produced by randomizing the original data.')
subplot(4,1,4);hold on;
plot(lags,mean(shuffle_xc),'k') % mean of the shuffle.
plot(lags,mean(shuffle_xc) + std(shuffle_xc)*2,'r') % 97% confidence interval
plot(lags,mean(shuffle_xc) - std(shuffle_xc)*2,'r') % 97% confidence interval
plot(lags,xc_orig,'b');grid on; ylabel('r')
legend('shuf xc', '2std', 'original')

%% one can then make a 'final' cross corr that is now in z scores from the
% null distribution...
new_xc = (xc_orig - mean(shuffle_xc))./std(shuffle_xc);
figure
plot(lags,new_xc,'m');grid on; ylabel('z')
title('Normalized Cross Corr')

% Again, the above is not really going to eliminate artifacts when there
% are more subtle non-stationarities in the data, but it's a good start and
% maybe enough. (note: this method for example is still corrupted if you
% forget to subtract the mean from one of the signals.
%
% Feel free to play around with the parameters and the data!
%