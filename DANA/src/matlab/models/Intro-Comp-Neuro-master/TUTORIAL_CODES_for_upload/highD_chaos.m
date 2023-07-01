% highD_chaos.m
%
% A firing-rate model network with coupled excitatory and inhibitory units
% whose connections are strong and sparsely random, designed to be in the 
% regime that produces chaotic activity.
%
%   The circuit is simulated twice, with the only difference in the initial
%   conditions of a single unit shifted by 1 in 10^14 of a Hz.
%
%   This code was used to produce Figures 7.13 and 7.14 of Chapter 7 and is 
%   required for Tutorial 8.4 of Chapter 8 in the textbook:
%   An Introductory Course in Computational Neuroscience,
%   by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;                      % clear all prior variables and parameters from memory

%% Set up the network parameters

% The next two lines set the random number generator to reproducibly use
% the same seed. This allows for regeneration of the identical circuit
% which is designated with random connectivity.
s = RandStream('mt19937ar','Seed',1)
RandStream.setGlobalStream(s);

NE = 200;           % Number of excitatory (E) units              
NI = 200;           % Number of inhibitory (I) units
gEE = 90;           % mean sum of E connection strengths to an E-unit
gEI = 80;           % mean sum of E connection strengths to an I-unit
gIE = 90;           % mean sum of I connection strengths to an E-unit
gII = 80;           % mean sum of I connection strengths to an I-unit

% define matrices of connections
WEE = double(2*rand(NE));       % E-to-E random from 0-2
WEI = double(2*rand(NE,NI));    % E-to-I random from 0-2
WIE = double(2*rand(NI,NE));    % I-to-E random from 0-2
WII = double(2*rand(NI));       % I-to-I random from 0-2

% Sparse connection probabilities (all identical at 1/20)
connEE_prob = 0.05;
connIE_prob = 0.05;
connEI_prob = 0.05;
connII_prob = 0.05;
% The matrices connEE to connII are binary, the same size as WEE to WII
connEE = rand(NE)<connEE_prob;
connIE = rand(NI,NE)<connIE_prob;
connEI = rand(NE,NI)<connEI_prob;
connII = rand(NI)<connII_prob;

% Now render the connectivity matrices sparse
WEE = WEE.*connEE;
WIE = WIE.*connIE;
WEI = WEI.*connEI;
WII = WII.*connII;

% Finally normalize them so that the mean connection strengths are as
% desired
for i = 1:NE
    WEE(:,i) = WEE(:,i)/mean(WEE(:,i));
    WIE(:,i) = WIE(:,i)/mean(WIE(:,i));
end
for i = 1:NI
    WEI(:,i) = WEI(:,i)/mean(WEI(:,i));
    WII(:,i) = WII(:,i)/mean(WII(:,i));
end
% The following lines ensure the mean excitatory and inhibitory inputs to a
% unit scale with gEE gEI gIE or gII.
WEE = WEE*gEE/NE;
WEI = WEI*gEI/NE;
WIE = WIE*gIE/NI;
WII = WII*gII/NI;

%% Single unit firing-rate curve parameters
rmax = double(100);     % maximum rate
Ith = double(10);       % input needed for half-maximum rate
Isigma = double(1);     % change in input needed for significant rate change

% firing rate curve is a sigmoid inline function
rate_fn = @(x) double(rmax./(1+exp(-(x-Ith)/Isigma)));

%% Simulation parameters
dt = 0.00001;                   % very small time-step for accuracy
tmax = 1;                       % maximum simulation time
t = 0:dt:tmax;                  % vector of time points
Nt = length(t);                 % number of time points
tauE = 0.010;                   % time constant for E-units
exp_dt_tauE = exp(-dt/tauE);    % used in the exponential Euler method
tauI = 0.005;                   % time constnat for I-units
exp_dt_tauI = exp(-dt/tauI);    % used in the exponential Euler method

rE1 = double(zeros(NE,Nt));     % matrix of rates of E-units
rI1 = double(zeros(NI,Nt));     % matrix of rates of I-units

rE1(:,1) = [[0:NE-1]*rmax/NE]'; % initialize E-units with uniform spread of rates
rI1(:,1) = [[0:NE-1]*rmax/NE]'; % initialize I-units with uniform spread of rates

%% Now simulate through time (first simulation)
for i = 2:Nt
    
    % Use the exponential Euler method for integration
    I_E = WEE'*rE1(:,i-1) - WIE'*rI1(:,i-1);    % input current to E-units
    rEinf = rate_fn(I_E);                       % send input through f-I curve
    rE1(:,i) = rEinf + (rE1(:,i-1) - rEinf)*exp_dt_tauE;    % update E-unit rates
    
    I_I = WEI'*rE1(:,i-1)  - WII'*rI1(:,i-1);   % input current to I-units
    rIinf = rate_fn(I_I);                       % send input through f-I curve
    rI1(:,i) = rIinf + (rI1(:,i-1) - rIinf)*exp_dt_tauI;    % update I-unit rates
    
end

%% Plot the results of the first simulation, selecting three units for visualization
figure(1)
clf
subplot(2,3,1)
plot(t,rE1(1,:),'k')
hold on
ylabel('Activity (Hz)')
title('E-unit #1')
subplot(2,3,2)
plot(t,rE1(50,:),'k')
title('E-unit #50')
hold on
subplot(2,3,3)
plot(t,rE1(100,:),'k')
hold on
title('E-unit #100')
drawnow

%% Now repeat the process with a miniscule difference in initial conditions 
rE2 = double(zeros(NE,Nt));
rI2 = double(zeros(NI,Nt));

rE2(:,1) = [[0:NE-1]*rmax/NE]';
rI2(:,1) = [[0:NE-1]*rmax/NE]';

rE2(1,1) = 1e-14;       % This is the change in initial conditions from zero

%% Now simulate through time (second simulation)
for i = 2:Nt
    
    I_E = double(WEE'*rE2(:,i-1) - WIE'*rI2(:,i-1));
    rEinf = rate_fn(I_E);
    rE2(:,i) = rEinf + (rE2(:,i-1) - rEinf)*exp_dt_tauE;
    
    I_I = double(WEI'*rE2(:,i-1)  - WII'*rI2(:,i-1));
    rIinf = rate_fn(I_I);
    rI2(:,i) = rIinf + (rI2(:,i-1) - rIinf)*exp_dt_tauI;
    
end

%% Now plot results from the second simulation as black
figure(1)
subplot(2,3,4)
plot(t,rE2(1,:),'k')
hold on
xlabel('Time (sec)')
ylabel('Activity (Hz)')

subplot(2,3,5)
plot(t,rE2(50,:),'k')
hold on
xlabel('Time (sec)')
subplot(2,3,6)
plot(t,rE2(100,:),'k')
hold on
xlabel('Time (sec)')


annotation('textbox',[0.06 0.96 0.05 0.05],'String','A1','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.36 0.96 0.05 0.05],'String','A2','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.64 0.96 0.05 0.05],'String','A3','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.06 0.48 0.05 0.05],'String','B1','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.36 0.48 0.05 0.05],'String','B2','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.64 0.48 0.05 0.05],'String','B3','LineStyle','none','FontSize',16,'FontWeight','bold')

%% Finally produce Figure 7.14 in the textbook, by measuring the 
%  time-evolution of the difference in firing rates due to the tiny change 
%  in initial conditions.
figure(2)
clf
subplot('Position',[0.2 0.6 0.75 0.32])
plot(t,mean(abs(rE2-rE1)),'k')      % absolute difference in rates
axis([0 0.5 0 60])
ylabel('Mean abs(\Deltar)')

subplot('Position',[0.2 0.12 0.75 0.32])
plot(t,log10(mean(abs(rE2-rE1))),'k')   % plot on a log scale to see exponenital rise
ylabel('Log_{10}[Mean abs(\Deltar)]')
xlabel('Time (sec)')
axis([0 0.5 -15 2.5])

annotation('textbox',[0.00 0.98 0.02 0.02],'String','A','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.00 0.52 0.02 0.02],'String','B','LineStyle','none','FontSize',16,'FontWeight','bold')
