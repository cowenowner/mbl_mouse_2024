% Figure_4_17.m
% This model is based on the two-compartment model of Pinsky and Rinzel (1994)
% The dendritic compartment produces calcium spikes, which couple to the
% somatic compartment to produce regular bursts of action potentials.
%
% Hyperpolarization-activated conductance is added to show its impact on
% the standard PR model.
%
% The code requires the two functions,
% PR_soma_gating and PR_dend_gating
% to be in the path.
%
% This code is used to generate figures 4.17 and 4.18 of Chapter 3 in the
% textbook:
% "An Introductory Course in Computational Neuroscience"
%
% by Paul Miller, Brandeis University, June 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
dt = 10e-6;
tmax=12;

%% Set up default parameters for plotting
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf

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
G_H = 0e-6*D_frac;          % hyperpolarization-activated conductance

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

G_H_vec = [0 5 10 15]*1e-9; % Values of Ih-conductance to use

%% Carry out the simulations through many trials
Ntrials = length(G_H_vec);
burstrate_vec = zeros(1,Ntrials);

for trial = 1:Ntrials
    G_H = G_H_vec(trial)
    
    for i = 2:length(t); % now see how things change through time
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
        ID(i) = I_LD(i)+I_Ca(i)+I_KCa(i)+I_KAHP(i)+I_H(i)-I_Link(i); % total current in dendrite
        
        gS_Tot = G_LS+G_Na_now+G_K_now+G_Link;
        VS_inf = (G_LS*E_L + G_Na_now*E_Na + G_K_now*E_K ...
            + VD(i-1)*G_Link )/gS_Tot;
        
        gD_Tot = G_LD+G_Ca_now+G_KCa_now+G_KAHP_now+G_H_now+G_Link;
        VD_inf = (G_LD*E_L + G_Ca_now*E_Ca + G_KCa_now*E_K + G_KAHP_now*E_K ...
            + G_H_now*E_H + VS(i-1)*G_Link )/gD_Tot;
        
        VS(i) = VS_inf - (VS_inf-VS(i-1))*exp(-dt*gS_Tot/CmS);  % Update the membrane potential, V.
        VD(i) = VD_inf - (VD_inf-VD(i-1))*exp(-dt*gD_Tot/CmD);  % Update the membrane potential, V.
        Ca_inf = tau_Ca*convert_Ca*I_Ca(i);
        Ca(i) = Ca_inf - (Ca_inf-Ca(i-1))*exp(-dt/tau_Ca);  % update Ca level
        
    end
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
    if ( Nbursts > 1 )
        burstrate_vec(trial) = 1/(tburst(Nbursts)-tburst(Nbursts-1));
    end
    
    %% Next section is for plotting the figures
    % Plot somatic membrane potential and dendritic membrane potential
    % one above the other for entire time window of 2 sec.
    figure(1)
    
    if ( trial == 1 )
        subplot('Position',[0.15 0.8 0.44 0.185])
        plot(t,VS*1000,'k')
        ylabel('V_S (mV)')
        axis([tmax-4 tmax -85 55])
        set(gca,'XTick',[8 9 10 11 12])
        
        subplot('Position',[0.65 0.8 0.31 0.185])
        plot(t,VS*1000,'k')
        axis([tburst(end-1)-0.025 tburst(end-1)+0.025 -85 55])
        set(gca,'XTick',...
            [tburst(end-1)-0.025, tburst(end-1), tburst(end-1)+0.025]);
        set(gca,'XTickLabel',{'-25', '0', '25'})
    end
    if ( trial == 2 )
        subplot('Position',[0.15 0.56 0.44 0.185])
        plot(t,VS*1000,'k')
        ylabel('V_S (mV)')
        axis([tmax-4 tmax -85 55])
        set(gca,'XTick',[8 9 10 11 12])
        
        subplot('Position',[0.65 0.56 0.31 0.185])
        plot(t,VS*1000,'k')
        axis([tburst(end-1)-0.025 tburst(end-1)+0.025 -85 55])
        set(gca,'XTick',...
            [tburst(end-1)-0.025, tburst(end-1), tburst(end-1)+0.025]);
        set(gca,'XTickLabel',{'-25', '0', '25'})
        
    end
    if ( trial == 3 )
        subplot('Position',[0.15 0.32 0.44 0.185])
        plot(t,VS*1000,'k')
        ylabel('V_S (mV)')
        axis([tmax-4 tmax -85 55])
        set(gca,'XTick',[8 9 10 11 12])
        
        subplot('Position',[0.65 0.32 0.31 0.185])
        plot(t,VS*1000,'k')
        axis([tburst(end-1)-0.025 tburst(end-1)+0.025 -85 55])
        set(gca,'XTick',...
            [tburst(end-1)-0.025, tburst(end-1), tburst(end-1)+0.025]);
        set(gca,'XTickLabel',{'-25', '0', '25'})
        
    end
    if ( trial == Ntrials )
        subplot('Position',[0.15 0.08 0.44 0.185])
        plot(t,VS*1000,'k')
        ylabel('V_S (mV)')
        xlabel('Time (sec)')
        axis([tmax-4 tmax -85 55])
        set(gca,'XTick',[8 9 10 11 12])
        
        subplot('Position',[0.65 0.08 0.31 0.185])
        plot(t,VS*1000,'k')
        xlabel('Time (msec)')
        
        set(gca,'YTickLabel',[])
        axis([tburst(end-1)-0.025 tburst(end-1)+0.025 -85 55])
        set(gca,'XTick',...
            [tburst(end-1)-0.025, tburst(end-1), tburst(end-1)+0.025]);
        set(gca,'XTickLabel',{'-25', '0', '25'})
    end
    drawnow
end;        % Go to next trial with new G_h
%% End of simulations
% Add labels to the graph (Fig 4.17)
annotation('textbox',[0 0.99 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.0 0.75 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.0 0.51 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.0 0.27 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','D')
annotation('textbox',[0.795 0.972 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','G_{h}=0')
annotation('textbox',[0.795 0.732 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','G_{h}=5nS')
annotation('textbox',[0.795 0.492 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','G_{h}=10nS')
annotation('textbox',[0.795 0.252 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','G_{h}=15nS')

%% Finally show the V-dependence of the gating variables for Ih
%  on a separate graph (Fig 3.17)
V = -0.100:0.001:0;     % Range of membrane potential

mH_inf = 1./(1+exp((V+0.070)/0.006));       % To plot m_h vs V
tau_mH = 0.272 + 1.499./(1 + exp(-(V+0.0422)/0.00873)); % to plot tau_m_h vs V

figure(2)
subplot('Position',[0.1 0.18 0.38 0.78])
plot(V*1000,mH_inf,'k')
xlabel('Membrane potential (mV)')
ylabel('m_{h} (steady state)')

subplot('Position',[0.6 0.18 0.38 0.78])
plot(V*1000,tau_mH,'k')
xlabel('Membrane potential (mV)')
ylabel('\tau_{h} (sec)')

