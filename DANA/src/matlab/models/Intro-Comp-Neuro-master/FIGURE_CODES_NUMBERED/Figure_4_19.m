% Figure_4_19.m
%
% This model is a variant with 3 dendritic and one somatic compartment, 
% based on the two-compartment model of Pinsky and Rinzel (1994)
%
% The simulation has a sequence of inputs arriving simultaneously either on
% the three different dendritic branches or on the same dendritic branch.
% 
% This code was used to produce Figure 4.19 in the textbook 
% An Introductory Course in Computational Neuroscience 
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set up default plotting parameters 
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
%% Set up parameters for the simulation
clear;              % Clear all variables and parameters from memory
dt = 10e-6;         % Time-step in sec (so 10 microsec)
tmax=2;             % Total simulation time in sec

I0 = 0;             % Baseline current
Itotal = 0.15e-9;   % Sets the current scale
IappSfrac = 0;      % Fraction of current to enter soma is zero
ND = 3;             % Number of dendritic branches

istart = 0.25;      % time applied current starts (sec)
ilength= 1.5;       % length of applied current pulse
isine_flag = 1;     % used to produce a series of regular pulses
freq = 3;           % frequency of pulses (Hz)

E_L = -0.075;   % leak reversal potential
E_Na = 0.060;   % reversal for sodium channels
E_K = -0.075;   % reversal for potassium channels
E_Ca = 0.080;   % reversal for calcium channels

S_frac = 1/4;           % fraction of total membrane area that is soma
D_frac = 1-S_frac;      % rest of area is dendritic
D_frac = D_frac/ND;     % fraction for each dendritic compartment

G_LD = 1e-7*D_frac;     % leak conductance in Siemens for each dendritic compartment
G_LS = 1e-7*S_frac;     % leak conductance in Siemens for Soma
G_Na = 2e-6*S_frac;     % sodium conductance in soma
G_K = 5e-6*S_frac;      % potassium conductance in soma

% Dendritic conductance values for a single dendritic compartment are
% calculated below
G_Ca = 4.5e-6*D_frac;       % calcium conductance in each dendritic compartment
G_KAHP = 0.00e-6*D_frac;    % potassium conductance to generate after-hyperpolarization
G_KCa = 8e-6*D_frac;        % calcium-dependent Potassium conductance

G_Link = 6e-9*ones(ND,1);   % conductance linking dendrite and soma

tau_Ca = 50e-3;             % time constant for removal of calcium from cell
convert_Ca = 5e5/D_frac;    % conversion changing calcium charge entry per unit area into concentration

CmS = 100e-12*S_frac;     % somatic membrane capacitance in Farads
CmD = 100e-12*D_frac;     % dendritic compartment membrane capacitance in Farads 

% Loop through two trials.
% Trial 1 has applied current separate across dendritic compartments
% Trial 2 has applied current entering a single compartment
for trial = 1:2
    if ( trial == 1 )
        IappDfrac = [0.5; 0.5; 0.5];
    else
        IappDfrac = [1; 0; 0];
    end
    t=0:dt:tmax;        % time vector
    Nt = length(t);     % Number of time points
    VS=zeros(size(t));  % somatic voltage vector
    VD=zeros(ND,Nt);    % dendritic voltage vector (one row per compartment)
    
    Ca=zeros(ND,Nt);    % dendritic calcium level (extra Ca above base level)
    Ca(:,1) = 0;        % initialize with no (extra) Ca in cell.
    
    I_LD= zeros(ND,Nt);     % leak current in dendritic compartments
    I_LS= zeros(size(t));   % leak current in soma
    I_Na = zeros(size(t));  % sodium current in soma
    I_K = zeros(size(t));   % potassium current in soma
    I_Ca = zeros(ND,Nt);    % calcium current (one row per dendritic compartment)
    I_KAHP = zeros(ND,Nt);  % AHP current (one row per dendritic compartment)
    I_KCa = zeros(ND,Nt);   % K_Ca current (one row per dendritic compartment)
    
    VS(1) = E_L;            % set the inititial value of somatic voltage
    VD(:,1) = -0.050;       % set the inititial value of all dendritic voltages
    
    n=zeros(size(t));   % n: potassium activation gating variable
    m=zeros(size(t));   % m: sodium activation gating variable
    h=zeros(size(t));   % h: sodim inactivation gating variplot(t,V)able
    n(1) = 0.4;
    h(1) = 0.5;
    
    % set up gating variables for dendrites as matrices (one row per
    % compartment)
    mca=zeros(ND,Nt);       % Ca current activation gating variable
    mkca=zeros(ND,Nt);      % K_Ca current activation gating variable
    mkahp = zeros(ND,Nt);   % AHP current activation gating variable
    mkahp(:,1) = 0.2;
    mca(:,1) = 0.02;
    mkca(:,1) = 0.02;
    Ca(:,1) = 1e-4;         % initial calcium concentration (Molar)
    
    %% The following section sets up the applied current as a series of pulses
    % baseline current is I0
    Iapp = zeros(size(t));
    for i = 1:round(istart/dt)
        Iapp(i) = I0;
    end
    % The following section adds a pulse whenever the cosine is above 0.5
    for i=round(istart/dt)+1:round((istart+ilength)/dt) % make non-zero for duration of current pulse
        if ( isine_flag ) % if the pulse is oscillatory
            Iapp(i) = I0 + 0.5*Itotal* ...
                (1 + sign(cos((t(i)-istart/dt)*freq*2*pi)-0.5));
        else
            Iapp(i) = Itotal;
        end
    end
    % Baseline current after the applied pulses
    for i = round((istart+ilength)/dt):length(Iapp)
        Iapp(i) = I0;
    end
    
    Itot=zeros(size(t)); % to store the total current
    
    for i = 2:length(t); % now see how things change through time
        I_LS(i) = G_LS*(E_L-VS(i-1));
        I_LD(:,i) = G_LD*(E_L-VD(:,i-1));
        
        % Values of variables at this time-point are extracted from arrays
        % for ease of viewing the code
        Vm = VS(i-1);       % allows for conversion of units 
        VmD = VD(:,i-1);    % allows for conversion of units 
        Catmp = Ca(:,i-1);
        mtmp = m(i-1);
        htmp = h(i-1);
        ntmp = n(i-1);
        mcatmp = mca(:,i-1);
        mkcatmp = mkca(:,i-1);
        mkahptmp = mkahp(:,i-1);
        
        % Find the rate constants using the functions PR_soma_gating and
        % PR_dend_gating
        [ alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n ] = PR_soma_gating(Vm);
        [ alpha_mca, beta_mca, alpha_mkca, beta_mkca, alpha_mkahp, beta_mkahp ] = PR_dend_gating(VmD, Catmp);
        
        % Update all the gating variables using the rate constants
        m(i) = mtmp + dt*( alpha_m.*(1-mtmp) - beta_m.*mtmp );
        h(i) = htmp + dt*( alpha_h.*(1-htmp) - beta_h.*htmp );
        n(i) = ntmp + dt*( alpha_n.*(1-ntmp) - beta_n.*ntmp );
        mca(:,i) = mcatmp + dt*( alpha_mca.*(1-mcatmp) - beta_mca.*mcatmp );
        mkca(:,i) = mkcatmp + dt*( alpha_mkca.*(1-mkcatmp) - beta_mkca.*mkcatmp );
        mkahp(:,i) = mkahptmp + dt*( alpha_mkahp.*(1-mkahptmp) - beta_mkahp.*mkahptmp );
        
        % Update each conductance and associated current using the gating
        % variables
        G_Na_now = G_Na*m(i)*m(i)*h(i);
        I_Na(i) = G_Na_now*(E_Na-VS(i-1)); % sodium current in soma
        
        G_K_now = G_K*n(i)*n(i);
        I_K(i) = G_K_now*(E_K-VS(i-1)); % potassium delayed rectifier current, soma
        
        G_Ca_now = G_Ca*mca(:,i).*mca(:,i);
        I_Ca(:,i) = G_Ca_now.*(E_Ca-VD(:,i-1)); % persistent sodium current in dendrite
        
        G_KCa_now = min(Ca(:,i-1),250e-6).*G_KCa.*mkca(:,i)/250e-6;
        I_KCa(:,i) = G_KCa_now.*(E_K-VD(:,i-1)); % calcium-dependent potassium current in dendrite
        
        G_KAHP_now = G_KAHP*mkahp(:,i);
        I_KAHP(:,i) = G_KAHP_now.*(E_K-VD(:,i-1)); % calcium-dependent potassium current in dendrite
        I_Link(:,i) = G_Link.*(VD(:,i-1)-VS(i-1));
        
        % Update applied currents
        IappD = Iapp(i)*IappDfrac;
        IappS = Iapp(i)*IappSfrac;
        
        % Total somatic current
        IS(i) = I_LS(i)+I_Na(i)+I_K(i)+sum(I_Link(:,i))+IappS; % total current in soma
        % Total dendritic current in each compartment
        ID(:,i) = I_LD(:,i)+I_Ca(:,i)+I_KCa(:,i)+I_KAHP(:,i)-I_Link(:,i)+IappD; % total current in dendrite
           
        % We will use the exponential method for integrating as it is
        % stable.
        % This requires calculating total conductance and the temporary
        % "steady state" membrane potential in each compartment
        gS_Tot = G_LS+G_Na_now+G_K_now+sum(G_Link);
        VS_inf = (G_LS*E_L + G_Na_now*E_Na + G_K_now*E_K ...
            + sum(VD(:,i-1).*G_Link) + IappS)/gS_Tot;
        
        gD_Tot = G_LD+G_Ca_now+G_KCa_now+G_KAHP_now+sum(G_Link);
        VD_inf = (G_LD*E_L + G_Ca_now*E_Ca + G_KCa_now*E_K + G_KAHP_now*E_K ...
            + VS(i-1)*G_Link + IappD )./gD_Tot;
        
        % Use a stable/robuse exponential method to step toward the current
        % "steady state" membrane potential in each compartment
        VS(i) = VS_inf - (VS_inf-VS(i-1))*exp(-dt*gS_Tot/CmS);  % Update the membrane potential, V.
        VD(:,i) = VD_inf - (VD_inf-VD(:,i-1)).*exp(-dt*gD_Tot/CmD);  % Update the membrane potential, V.

        % Update calcium with the same method
        Ca_inf = tau_Ca*convert_Ca*I_Ca(:,i);
        Ca(:,i) = Ca_inf - (Ca_inf-Ca(:,i-1)).*exp(-dt/tau_Ca);  % update Ca level
        
    end
    
    % figure(1) is Figure 4.19 in the textbook
    figure(1)
    if ( trial == 1 )
        clf
        % top-left panel is somatic potential, Trial 1
        subplot('Position',[0.11 0.56 0.35 0.33])
        plot(t,VS*1000,'k')
        hold on
        axis([0 2 -85 50])
        ylabel('V_S (mV)')
        title('I_{App1} = I_{App2} = I_{App3} = I_{0}/2')
        
        % bottom-left panel is first dendritic compartment, Trial 1
        subplot('Position',[0.11 0.14 0.35 0.33])
        plot(t,VD(1,:)*1000,'k')
        
        axis([0 2 -85 50])
        ylabel('V_{D1} (mV)')
        xlabel('Time (sec)')
        
    else
        % top-right panel is somatic potential, Trial 2
        subplot('Position',[0.62 0.56 0.35 0.33])
        plot(t,VS*1000,'k')
        hold on
        axis([0 2 -85 50])
        title('I_{App1} = I_{0}; I_{App2} = I_{App3} = 0')
        ylabel('V_S (mV)')

        % bottom-right panel is first dendritic compartment, Trial 2
        subplot('Position',[0.62 0.14 0.35 0.33])
        plot(t,VD(1,:)*1000,'k')       
        axis([0 2 -85 50])
        xlabel('Time (sec)')
        ylabel('V_{D1} (mV)')
        
    end
    
end

% Finally label each panel in the figure
annotation('textbox',[0 0.95 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B1')
annotation('textbox',[0.0 0.52 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','B2')
annotation('textbox',[0.5 0.95 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C1')
annotation('textbox',[0.5 0.52 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','C2')
