% Figure_5_11C.m
% This model is based on the two-compartment model of Pinsky and Rinzel (1994)
% The dendritic compartment produces calcium spikes, which couple to the
% somatic compartment to produce regular bursts of action potentials.
%
% The code requires the two functions,
% PR_soma_gating and PR_dend_gating
% to be in the path.
%
% This code is used to produce Figure 5.11C of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
dt = 2e-6;
tmax=2;

%% Set up parameters for the neuron
E_L = -0.060;   % leak reversal potential
E_Na = 0.060;   % reversal for sodium channels
E_K = -0.075;   % reversal for potassium channels
E_Ca = 0.080;   % reversal for calcium channels
E_H = -0.020;   % reversal of hyperpolarization-activated current

S_frac = 1/3;  % fraction of total membrane area that is soma
D_frac = 1-S_frac; % rest of area is dendritic

% Conductance values for somatic channels follow

% Conductance values for somatic channels follow
G_LS = 1e-9*S_frac;     % somatic leak conductance in Siemens
G_Na = 3e-6*S_frac;     % sodium conductance (Soma)
G_K = 2e-6*S_frac;      % potassium conductance (Soma)

% Conductance values for dendritic channels follow
G_LD = 1e-9*D_frac;         % dendritic leak conductance in Siemens
G_Ca = 2.5e-6*D_frac;         % calcium conductance (Dendrite)
G_KAHP = 0.06e-6*D_frac;     % Potassium conductance to generate after-hyperpolarization
G_KCa = 5e-6*D_frac;        % calcium-dependent Potassium conductance
G_H = 15e-9*D_frac;          % hyperpolarization-activated conductance

G_Link = 25e-9; % conductance linking dendrite and soma

tau_Ca = 50e-3;             % time constant for buffering of calcium
convert_Ca = 1e6/D_frac;  % conversion changing calcium charge entry per unit area into concentration

CmS = 100e-12*S_frac;     % somatic membrane capacitance in Farads
CmD = 100e-12*D_frac;     % dendritic membrane capacitance in Farads

t=0:dt:tmax;        % time vector
VS=zeros(size(t));  % somatic voltage vector
VD=zeros(size(t));  % dendritic voltage vector
VS(1) = E_L;    % set the inititial value of somatic voltage
VD(1) = E_L;    % set the inititial value of dendritic voltage


Ca=zeros(size(t));  % dendritic calcium level (extra Ca above base level)
Ca(1) = 0;          % initialize with no (extra) Ca in cell.

I_LD= zeros(size(t));   % leak current in dendrite
I_LS= zeros(size(t));   % leak current in soma
I_Na = zeros(size(t));  % sodium current (soma)
I_K = zeros(size(t));   % potassium current (soma)
I_Ca = zeros(size(t));  % calcium current (dendrite)
I_KAHP = zeros(size(t)); % after-hyperpolarization current (dendrite)
I_KCa = zeros(size(t)); % calcium-dependent potassium current (dendrite)
I_H = zeros(size(t));   % hyperpolarization-activated current (dendrite)

n=zeros(size(t));   % n: potassium activation gating variable
m=zeros(size(t));   % m: sodium activation gating variable
h=zeros(size(t));   % h: sodim inactivation gating variplot(t,V)able
n(1) = 0.4;         % initialize near steady state at resting potential
h(1) = 0.5;         % initialize near steady state at resting potential

mca=zeros(size(t));     % Ca current activation gating variable
mkca=zeros(size(t));    % K_Ca current activation gating variable
mkahp = zeros(size(t)); % K_AHP current activation gating variable
mH = zeros(size(t));    % I_H activation gating variable
mkahp(1) = 0.2;         % initialize near steady state at resting potential
mkca(1) = 0.2;          % initialize near steady state at resting potential
Ca(1) = 1e-6;           % initialize near steady state at resting potential
mH(1) = 0.5;            % initialize near steady state at resting potential

Itot=zeros(size(t)); % in case we want to plot and look at the total current

%% Now add current pulses at different points on the cycle and analyze the
%  change in response due to the pulse.

Ipulse_amp = 10e-12;

Npulses = 100;
Ipulse = 0;
shift = zeros(1,Npulses);

for trial = 0:Npulses;  % INitially trial 0 is used to get the baseline
    
    if ( trial == 1 )
        i_tpulse = i_startphase + floor((i_stopphase-i_startphase)*[1:Npulses]/Npulses);
    end
    
    for i = 2:length(t); % now see how things change through time
        
        if ( trial > 0 )
            if ( ( i >= i_tpulse(trial) ) && ( i < i_tpulse(trial)+0.005/dt ) )
                Ipulse = Ipulse_amp;
            else
                Ipulse = 0;
            end
        end
        I_LS(i) = G_LS*(E_L-VS(i-1));
        I_LD(i) = G_LD*(E_L-VD(i-1));
        
        % Take variables from last time-point to update all variables in
        % this time-point ("tmp" stands for temporary)
        Vm = VS(i-1);
        VmD = VD(i-1);
        Catmp = Ca(i-1);
        mtmp = m(i-1);
        htmp = h(i-1);
        ntmp = n(i-1);
        mcatmp = mca(i-1);
        mkcatmp = mkca(i-1);
        mkahptmp = mkahp(i-1);
        mHtmp = mH(i-1);
        
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        [ alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n ] = PR_soma_gating(Vm);
        [ alpha_mca, beta_mca, alpha_mkca, beta_mkca, alpha_mkahp, beta_mkahp ] = PR_dend_gating(VmD, Catmp);
        
        % Update somatic gating variables based on rate constants
        m(i) = mtmp + dt*( alpha_m*(1-mtmp) - beta_m*mtmp );
        h(i) = htmp + dt*( alpha_h*(1-htmp) - beta_h*htmp );
        n(i) = ntmp + dt*( alpha_n*(1-ntmp) - beta_n*ntmp );
        
        % Update dendritic gating variables based on rate constants
        mca(i) = mcatmp + dt*( alpha_mca*(1-mcatmp) - beta_mca*mcatmp );
        mkca(i) = mkcatmp + dt*( alpha_mkca*(1-mkcatmp) - beta_mkca*mkcatmp );
        mkahp(i) = mkahptmp + dt*( alpha_mkahp*(1-mkahptmp) - beta_mkahp*mkahptmp );
        
        % Gating variables for hyperpolarization-activated conductance, G_H
        mH_inf = 1/(1+exp((VmD+0.070)/0.006));
        tau_mH = 0.272 + 1.499/(1 + exp(-(VmD+0.0422)/0.00873));
        mH(i) = mH(i-1) + (mH_inf-mH(i-1))*dt/tau_mH;    % Update mH
        
        % Now update all conductances and currents
        G_Na_now = G_Na*m(i)*m(i)*h(i);     % instantaneous sodium conductance
        I_Na(i) = G_Na_now*(E_Na-VS(i-1));  % sodium current in soma
        
        G_K_now = G_K*n(i)*n(i);            % instantaneous potassium conductance
        I_K(i) = G_K_now*(E_K-VS(i-1));     % potassium delayed rectifier current, soma
        
        G_Ca_now = G_Ca*mca(i)*mca(i);      % instantaneous calcium conductance
        I_Ca(i) = G_Ca_now*(E_Ca-VD(i-1));  % caclium current in dendrite
        
        % G_KCa depends on both [Ca] and V
        if ( Ca(i-1) > 250e-6 )
            G_KCa_now = G_KCa*mkca(i);
        else
            G_KCa_now = G_KCa*mkca(i)*Ca(i-1)/250e-6;
        end
        I_KCa(i) = G_KCa_now*(E_K-VD(i-1)); % calcium-dependent potassium current in dendrite
        
        G_KAHP_now = G_KAHP*mkahp(i);           % K_AHP instantaneous conductance
        I_KAHP(i) = G_KAHP_now*(E_K-VD(i-1));   % after-hyperoloarization potassium current in dendrite
        
        G_H_now = G_H*mH(i);            % Ih instantaneous conductance
        I_H(i) = G_H_now*(E_H-VD(i-1)); % hyperpolarization-activated mixed cation current in dendrite
        
        I_Link(i) = G_Link*(VD(i-1)-VS(i-1));   % current from dendrite to soma
        
        IS(i) = I_LS(i)+I_Na(i)+I_K(i)+I_Link(i); % total current in soma
        ID(i) = I_LD(i)+I_Ca(i)+I_KCa(i)+I_KAHP(i)+I_H(i)-I_Link(i)+Ipulse; % total current in dendrite
        
        gS_Tot = G_LS+G_Na_now+G_K_now+G_Link;
        VS_inf = (G_LS*E_L + G_Na_now*E_Na + G_K_now*E_K ...
            + VD(i-1)*G_Link )/gS_Tot;
        
        gD_Tot = G_LD+G_Ca_now+G_KCa_now+G_KAHP_now+G_H_now+G_Link;
        VD_inf = (G_LD*E_L + G_Ca_now*E_Ca + G_KCa_now*E_K + G_KAHP_now*E_K ...
            + G_H_now*E_H + VS(i-1)*G_Link + Ipulse)/gD_Tot;
        
        VS(i) = VS_inf - (VS_inf-VS(i-1))*exp(-dt*gS_Tot/CmS);  % Update the membrane potential, V.
        VD(i) = VD_inf - (VD_inf-VD(i-1))*exp(-dt*gD_Tot/CmD);  % Update the membrane potential, V.
        Ca_inf = tau_Ca*convert_Ca*I_Ca(i);
        Ca(i) = Ca_inf - (Ca_inf-Ca(i-1))*exp(-dt/tau_Ca);  % update Ca level
        
    end
    
    
    %% Detect the initiation time of individual bursts
    inburst = 0;
    tburst = [];
    Nbursts = 0;
    for i = 1:length(t)
        if ( inburst == 0 ) && ( VD(i) > 0 )
            inburst = 1;
            tburst(end+1) = t(i);
            Nbursts = Nbursts + 1;
        end
        if (inburst == 1 ) && ( VD(i) < -0.050 )
            inburst = 0;
        end
    end
    %% Now decide where a phase of "zero" corresponds to in the oscillation
    %  and calculate the time where this occurs -- time should be well after any
    %  initial transients.
    if ( trial == 0 )
        i_startphase = round(tburst(3)/dt);
        i_stopphase = round(tburst(4)/dt);
        period = tburst(4) - tburst(3);
    else
        shift(trial) = 2*pi*(1 - (tburst(4) - tburst(3))/period );
        if ( shift(trial) < -pi )
            shift(trial) = shift(trial) + 2*pi;
        end
        if ( shift(trial) > pi )
            shift(trial) = shift(trial) - 2*pi;
        end
        
        data = shift(trial)
    end
    %% Now plot the graphs
    if ( mod(trial,10) == 0 )
        set(0,'DefaultLineLineWidth',2,...
            'DefaultLineMarkerSize',8, ...
            'DefaultAxesLineWidth',2, ...
            'DefaultAxesFontSize',14,...
            'DefaultAxesFontWeight','Bold');
        
        % First plot somatic membrane potential and dendritic membrane potential
        % one above the other for entire time window of 2 sec.
        figure(trial+1)
        clf
%         subplot('Position',[0.16 0.57 0.8 0.39])
        plot(t,VS,'k')
        hold on
        axis([tburst(3)-0.05 tburst(4)+0.05 -0.075 0.050])
        ylabel('Membrane Potential (V)')
        
%         subplot('Position',[0.16 0.13 0.8 0.39])
%         plot(t,VD*1000,'k')
%         
%         axis([1 2 -85 50])
%         ylabel('V_D (mV)')
        xlabel('Time (sec)')
        drawnow
        
    end
end


figure()
plot([1:Npulses]*2*pi/Npulses,shift,'k')
hold on
plot([1:Npulses]*2*pi/Npulses,[1:Npulses]*0,'k:')
xlabel('Phase of pulse')
ylabel('Phase shift')
axis([0 2*pi -0.005 0.06])
set(gca,'YTick',[ 0 0.02 0.04])
set(gca,'XTick',[0 pi 2*pi])
set(gca,'XTickLabel',{'0' '\pi' '2\pi'})


