%% LFP demo
sFreq = 500;
sig_freq = 10;
sig_freq2 = 30;
dur_sec = 10;
t = 0:1/sFreq:dur_sec;
% create radiance at the freq you want
v = linspace(0,2*pi*sig_freq*dur_sec,sFreq*dur_sec); % 2 pi is the ful cycle 360 degrees
v2 = linspace(0,2*pi*sig_freq2*dur_sec, sFreq*dur_sec); % 2 pi is the ful cycle 360 degrees
LFP = sin(v);
LFP_withnoise = sin(v) + randn(size(v));
LFP_twofreqs = sin(v) + cos(v2) + randn(size(v));

figure
plot(LFP_twofreqs)
figure
pwelch(LFP_twofreqs,[],[],[],sFreq)
hold on
pwelch(LFP_twofreqs,sFreq,sFreq/2,[],sFreq)