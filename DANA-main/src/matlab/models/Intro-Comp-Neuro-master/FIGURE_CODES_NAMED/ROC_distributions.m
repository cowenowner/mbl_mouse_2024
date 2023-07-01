% ROC_distributions.m
% 
% This code was used to produce Figure 3.10 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
dx = 0.1;
x = 0:0.001:50;       % Vector of possible response values

%% First produce base probability distribution (Distribution-a)
% Base process (0) has mean response of 15Hz and standard deviation of 5Hz
mu_0 = 15;          % mean
sigma_0 = 5;        % standard deviation
% Now produce a Gaussian (Distribution-a) for the base process
P_a = exp(-(x-mu_0).*(x-mu_0)./(2*sigma_0*sigma_0));
P_a = P_a / (dx*sum(P_a));  % Integral of probability distribution must be 1.

%% Distribution-b is a mixture of the base process and a higher rate process
% Process (1) is a higher rate process
mu_1 = 40;
sigma_1 = 5;

p1 = 0.3;           % Contribution of process 1 to Distribution-b
p0 = 1-p1;          % Contribution of base process to Distribution-b.

% Now produce the mixture of two Gaussians for Distribution-b.
P_b = (p0/sigma_0)*exp(-(x-mu_0).*(x-mu_0)./(2*sigma_0*sigma_0)) ...
    + (p1/sigma_1)*exp(-(x-mu_1).*(x-mu_1)./(2*sigma_1*sigma_1));
P_b = P_b / (dx*sum(P_b));  % Integral of probability distribution must be 1.

%% Distribution-c is a mixture of two different higher-rate Gaussians
mu_2 = 20;          % Recognition (process-2) mean-rate
sigma_2 = 5;
mu_3 = 40;          % Recall (process-3) mean-rate
sigma_3 = 5;

p3 = 0.25;          % Fraction of contribution from process-3
p2 = 1-p3;          % Fraction of contribution from process-2

% Now produce the mixture of two Gaussians for Distribution-c.
P_c = (p2/sigma_2)*exp(-(x-mu_2).*(x-mu_2)./(2*sigma_2*sigma_2)) ...
    + (p3/sigma_3)*exp(-(x-mu_3).*(x-mu_3)./(2*sigma_3*sigma_3));
P_c = P_c / (dx*sum(P_c));  % Integral of probability distribution must be 1.


%% Distribution-d is a single Gaussian (one process) with greater standard deviation
mu_4 = 23;          % Mean rate from process-4
sigma_4 = 8;        % Greater standard deviation in rate for process-4

% Produce Distribution-d as a Gaussian entirely from process-r
P_d = exp(-(x-mu_4).*(x-mu_4)./(2*sigma_4*sigma_4));
P_d = P_d / (dx*sum(P_d));  % Integral of probability distribution must be 1.

% Calculate the probability as a function of threshold (x) that a draw from
% each of the distributions is greater than the threshold. This is
% equivalent to a false positive when drawn from the baseline distribution
% (no stimulus) or a true positive when drawn from each of the three
% stimulus-present distributions.
false_p = 1-cumsum(P_a)/sum(P_a);
true_p_b = 1 - cumsum(P_b)/sum(P_b);
true_p_c = 1 - cumsum(P_c)/sum(P_c);
true_p_d = 1 - cumsum(P_d)/sum(P_d);

%% Set up the plotting parameters and plot the distibutions and ROC curves
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf

%% First plot the individual probability distributions, each of the 
%  stimulus-present ones (P_b, P_c, P_d) with the baseline,
%  stimulus-absent distribution (P_a). 
subplot('Position',[0.08 0.72 0.23 0.23])
plot(x,P_a,'k--')
hold on
plot(x,P_b,'k')
xlabel('Response')
ylabel('Probability')
set(gca,'YTick',[])
set(gca,'XTick',[])
legend('No stim', 'Stim')

subplot('Position',[0.4 0.72 0.23 0.23])
plot(x,P_a,'k--')
hold on
plot(x,P_c,'k')
ylabel('Probability')
xlabel('Response')
set(gca,'YTick',[])
set(gca,'XTick',[])
legend('No stim', 'Stim')

subplot('Position',[0.72 0.72 0.23 0.23])
plot(x,P_a,'k--')
hold on
plot(x,P_d,'k')
ylabel('Probability')
xlabel('Response')
set(gca,'YTick',[])
set(gca,'XTick',[])
legend('No stim', 'Stim')

%% Next plot the ROC curve, true-positives versus false-positives
subplot('Position',[0.08 0.4 0.23 0.23])
plot(false_p,true_p_b,'k')
hold on
plot([0 1], [0 1],'k:')
xlabel('P(False +ve)')
ylabel('P(True +ve)')
set(gca,'XTick',[0 1])
set(gca,'YTick',[0 1])

subplot('Position',[0.4 0.4 0.23 0.23])
plot(false_p,true_p_c,'k')
hold on
plot([0 1], [0 1],'k:')
xlabel('P(False +ve)')
ylabel('P(True +ve)')
set(gca,'XTick',[0 1])
set(gca,'YTick',[0 1])

subplot('Position',[0.72 0.4 0.23 0.23])
plot(false_p,true_p_d,'k')
hold on
plot([0 1], [0 1],'k:')
xlabel('P(False +ve)')
ylabel('P(True +ve)')
set(gca,'XTick',[0 1])
set(gca,'YTick',[0 1])

%% Finally plot the z-scores of true positives versus z-scores of false positives
subplot('Position',[0.08 0.08 0.23 0.23])
plot(sqrt(2)*erfcinv(2*(1-false_p)),sqrt(2)*erfcinv(2*(1-true_p_b)),'k')
axis([-3 3 -3 3])
hold on
plot([-3 3], [-3 3],'k:')
xlabel('Z(False +ve)')
ylabel('Z(True +ve)')
set(gca,'YTick',[-3 0 3])
set(gca,'XTick',[-3 0 3])

subplot('Position',[0.4 0.08 0.23 0.23])
plot(sqrt(2)*erfcinv(2*(1-false_p)),sqrt(2)*erfcinv(2*(1-true_p_c)),'k')
axis([-3 3 -3 3])
hold on
plot([-3 3], [-3 3],'k:')
xlabel('Z(False +ve)')
ylabel('Z(True +ve)')
set(gca,'YTick',[-3 0 3])
set(gca,'XTick',[-3 0 3])

subplot('Position',[0.72 0.08 0.23 0.23])
plot(sqrt(2)*erfcinv(2*(1-false_p)),sqrt(2)*erfcinv(2*(1-true_p_d)),'k')
axis([-3 3 -3 3])
hold on
plot([-3 3], [-3 3],'k:')
xlabel('Z(False +ve)')
ylabel('Z(True +ve)')
set(gca,'YTick',[-3 0 3])
set(gca,'XTick',[-3 0 3])


annotation('textbox',[0 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A1')
annotation('textbox',[0.32 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A2')
annotation('textbox',[0.64 0.98 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A3')
annotation('textbox',[0 0.66 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B1')
annotation('textbox',[0.32 0.66 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B2')
annotation('textbox',[0.64 0.66 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B3')
annotation('textbox',[0 0.34 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C1')
annotation('textbox',[0.32 0.34 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C2')
annotation('textbox',[0.64 0.34 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C3')
