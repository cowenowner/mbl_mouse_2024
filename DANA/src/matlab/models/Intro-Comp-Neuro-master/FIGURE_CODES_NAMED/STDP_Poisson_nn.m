% STDP_Poisson_nn.m
% Plots the analytical result for two Poisson spike trains, presynaptic at
% a fixed rate, r_i, postsynaptic at a variable rate, r_j.
% The result is for an STDP window comprising two exponential decays, with
% amplitudes A_plus and A_minus respectively for potentiation
% (pre-before-post) and depression (post-before-pre) with corresponding
% time constants of tau_plus and tau_minus respectively.
%
% The resulting formula for the expected (mean) rate of change of synaptic
% strength is:
% <dW_ij/dt> = r_i.r_j.[ A_plus.tau_plus/(1+r_i.tau_plus)
%                       - A_minus.tau_minus/(1+r_j.tau_minus) ]
%
% The code produces Figure 8.13 of Chapter 8, Appendix A, in the book
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_i = 20;               % presynaptic rate
r_j = 0:0.1:40;         % range of postsynaptic rates
tau_plus = 0.020;       % time constant for STDP potentiation window
tau_minus = 0.020;      % time constant for STDP depression window
A_plus = 0.01;          % maximum amplitude of potentiation from spike-pair
A_minus = 0.011;        % maximum amplotude of depression from spike-pair 

% Produce curve of dW/dt from formula (see Eq. 8.16 in textbook)
dW_by_dt = r_i*r_j.*(A_plus*tau_plus/(1 + r_i*tau_plus) ...
    - A_minus*tau_minus./(1+r_j*tau_minus) );

%% Now plot the figure
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf
plot(r_j,dW_by_dt,'k')
hold on
plot(r_j,zeros(size(r_j)),'k:')
xlabel('Postsynaptic rate (Hz)')
ylabel('Net plasticity dW_{ij}/dt')
