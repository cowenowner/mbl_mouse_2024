% Figure_8_2.m
%
% Code to plot parabola r.(r-r_T) for a threshold r_T.
%
% This code produces Figure 8.2 of 
% An Introductory Course in Computational Neuroscience 
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_pre = 10;                 % presynaptic rate
r_post = 0:0.1:25;          % set of postsynaptic rates
r_T = 15;                   % threshold rate
epsilon = 0.001;            % scales down plasticity rate
% dWdt is rate of plasticity as a function of r_post
dWdt = epsilon*r_pre*r_post.*(r_post-r_T);    

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf
fill([0 0 r_T r_T], [-5 5 5 -5],[0.75, 0.75, 0.75], ...
    'EdgeColor','none')
hold on

% Plot plasticity rate
plot(r_post,dWdt,'k')
hold on
%Plot dotted line to indicate zero crossings
plot(r_post,zeros(size(r_post)),'k:')

annotation('textbox',[0.22 0.85 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','Depression')
annotation('textbox',[0.58 0.85 0.1 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','Potentiation')

set(gca,'YTick',[0 1 2])
set(gca,'Layer','top')
xlabel('Postsynaptic activity, r_{j} (Hz)')
ylabel('Plasticity rate, dW_{ij}/dt')
axis([0 25 -0.7 2.5])
