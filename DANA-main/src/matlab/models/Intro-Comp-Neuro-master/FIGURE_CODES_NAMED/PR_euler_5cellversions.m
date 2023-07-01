% PR_euler_5cellversions.m
%
% This model is the two-compartment model of Pinsky and Rinzel (1994)
%
% Five different levels of Ca conductance are used to demonstrate its
% effect on bursting.
%
% This code was used to produce Figure 8.9 in the textbook:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

figure(1)
clf

dt = 10e-6;         % time-step (sec)
tmax=2;             % maximum time to simulate (sec)

% 5 different types of cell with different calcium conductance will be
% simulated
for cell_version = 1:5
    
    E_L = -0.060;   % leak reversal potential
    E_Na = 0.060;   % reversal for sodium channels
    E_K = -0.075;   % reversal for potassium channels
    E_Ca = 0.080;   % reversal for calcium channels
    
    S_frac = 1/3;  % fraction of total membrane area that is soma
    D_frac = 1-S_frac; % rest of area is dendritic
    
    G_LD = 5e-9*D_frac;     % dendritic leak conductance in Siemens 
    G_LS = 5e-9*S_frac;     % somatic leak conductance in Siemens
    G_Na = 3e-6*S_frac;     % sodium conductance in Siemens
    G_K = 2e-6*S_frac;      % potassium conductance in Siemens
    % 5 different levels of calcium conductance are used as given below
    switch cell_version
        case 1
            G_Ca = 1e-6*D_frac;     % calcium conductance
        case 2
            G_Ca = 1.5e-6*D_frac;   % calcium conductance
        case 3
            G_Ca = 2e-6*D_frac;     % calcium conductance
        case 4
            G_Ca = 2.5e-6*D_frac;   % calcium conductance
        case 5
            G_Ca = 3e-6*D_frac;     % calcium conductance
    end
    G_KAHP = 0.04e-6*D_frac;        % Potassium conductance to generate after-hyperpolarization
    G_KCa = 2.5e-6*D_frac;          % calcium-dependent Potassium conductance
    
    G_Link = 20e-9; % conductance linking dendrite and soma divided by total membrane area
    
    tau_Ca = 50e-3;       % time constant for removal of calcium from cell
    convert_Ca = 0.25e7/D_frac;  % conversion changing calcium charge entry per unit area into concentration
    
    CmS = 100e-12*S_frac;     % specific membrane capacitance in Farads per mm-square
    CmD = 100e-12*D_frac;     % specific membrane capacitance in Farads per mm-square
    
    t=0:dt:tmax;        % time vector
    VS=zeros(size(t));  % somatic voltage vector
    VD=zeros(size(t));  % dendritic voltage vector
    VD(1) = E_L;
    VS(1) = E_L;
    
    Ca=zeros(size(t));  % dendritic calcium level (extra Ca above base level)
    Ca(1) = 0;          % initialize with no (extra) Ca in cell.
    
    I_LD= zeros(size(t));       % leak current in dendrite
    I_LS= zeros(size(t));       % leak current in soma
    I_Na = zeros(size(t));      % sodium current
    I_K = zeros(size(t));       % potassium current
    I_Ca = zeros(size(t));      % calcium current
    I_KAHP = zeros(size(t));    % AHP current
    I_KCa = zeros(size(t));     % calcium-dependent potassium current
    
    VS(1) = E_L;    % set the inititial value of somatic voltage
    VD(1) = E_L;    % set the inititial value of dendritic voltage
    
    n=zeros(size(t));   % n: potassium activation gating variable
    m=zeros(size(t));   % m: sodium activation gating variable
    h=zeros(size(t));   % h: sodim inactivation gating variplot(t,V)able
    n(1) = 0.4;
    h(1) = 0.5;
    
    mca=zeros(size(t)); % CaT current activation gating variable
    mkca=zeros(size(t)); % CaT current inactivation gating variable
    mkahp = zeros(size(t));
    mkahp(1) = 0.2;
    mkca(1) = 0.2;
    Ca(1) = 1e-6;
        
    Itot=zeros(size(t)); % in case we want to plot and look at the total current    
    
    %% Finished initialization, now simulate through time
    for i = 2:length(t); % now see how things change through time
        I_LS(i) = G_LS*(E_L-VS(i-1));
        I_LD(i) = G_LD*(E_L-VD(i-1));
        
        % values from last time-step are given as scalar variables for ease
        % of reading code
        Vm = VS(i-1); 
        VmD = VD(i-1); 
        Catmp = Ca(i-1);
        mtmp = m(i-1);
        htmp = h(i-1);
        ntmp = n(i-1);
        mcatmp = mca(i-1);
        mkcatmp = mkca(i-1);
        mkahptmp = mkahp(i-1);
        
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        % PR_soma_gating and PR_dend_gating are functions in separate files
        [ alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n ] = PR_soma_gating(Vm);
        [ alpha_mca, beta_mca, alpha_mkca, beta_mkca, alpha_mkahp, beta_mkahp ] = PR_dend_gating(VmD, Catmp);
        
        m(i) = mtmp + dt*( alpha_m*(1-mtmp) - beta_m*mtmp );
        h(i) = htmp + dt*( alpha_h*(1-htmp) - beta_h*htmp );
        n(i) = ntmp + dt*( alpha_n*(1-ntmp) - beta_n*ntmp );
        
        mca(i) = mcatmp + dt*( alpha_mca*(1-mcatmp) - beta_mca*mcatmp );
        mkca(i) = mkcatmp + dt*( alpha_mkca*(1-mkcatmp) - beta_mkca*mkcatmp );
        mkahp(i) = mkahptmp + dt*( alpha_mkahp*(1-mkahptmp) - beta_mkahp*mkahptmp );
        
        % Update all conductances and currents using gating variables
        G_Na_now = G_Na*m(i)*m(i)*h(i);
        I_Na(i) = G_Na_now*(E_Na-VS(i-1)); % sodium current in soma
        
        G_K_now = G_K*n(i)*n(i);
        I_K(i) = G_K_now*(E_K-VS(i-1)); % potassium delayed rectifier current, soma
        
        G_Ca_now = G_Ca*mca(i)*mca(i);
        I_Ca(i) = G_Ca_now*(E_Ca-VD(i-1)); % persistent sodium current in dendrite
        
        if ( Ca(i-1) > 250e-6 )
            G_KCa_now = G_KCa*mkca(i);
        else
            G_KCa_now = G_KCa*mkca(i)*Ca(i-1)/250e-6;
        end
        I_KCa(i) = G_KCa_now*(E_K-VD(i-1)); % calcium-dependent potassium current in dendrite
        
        G_KAHP_now = G_KAHP*mkahp(i);
        I_KAHP(i) = G_KAHP_now*(E_K-VD(i-1)); % calcium-dependent potassium current in dendrite
        I_Link(i) = G_Link*(VD(i-1)-VS(i-1));
                
        IS(i) = I_LS(i)+I_Na(i)+I_K(i)+I_Link(i); % total current in soma
        ID(i) = I_LD(i)+I_Ca(i)+I_KCa(i)+I_KAHP(i)-I_Link(i); % total current in dendrite
        
        gS_Tot = G_LS+G_Na_now+G_K_now+G_Link;
        VS_inf = (G_LS*E_L + G_Na_now*E_Na + G_K_now*E_K ...
            + VD(i-1)*G_Link )/gS_Tot;
        
        gD_Tot = G_LD+G_Ca_now+G_KCa_now+G_KAHP_now+G_Link;
        VD_inf = (G_LD*E_L + G_Ca_now*E_Ca + G_KCa_now*E_K + G_KAHP_now*E_K ...
            + VS(i-1)*G_Link )/gD_Tot;
        
        VS(i) = VS_inf - (VS_inf-VS(i-1))*exp(-dt*gS_Tot/CmS);  % Update the membrane potential, V.
        VD(i) = VD_inf - (VD_inf-VD(i-1))*exp(-dt*gD_Tot/CmD);  % Update the membrane potential, V.
        Ca_inf = tau_Ca*convert_Ca*I_Ca(i);
        Ca(i) = Ca_inf - (Ca_inf-Ca(i-1))*exp(-dt/tau_Ca);  % update Ca level
                
    end
    %% Plot the figures
    figure(1)
    subplot('Position',[0.18 1.035-0.19*cell_version 0.8 0.13])
    plot(t,VS*1000,'k')
    hold on
    axis([0 2 -85 50])
    ylabel('V_S (mV)')
    if ( cell_version == 5 )
        xlabel(' Time (sec)')
    end
    
    figure(2)
    subplot('Position',[0.18 1.035-0.19*cell_version 0.8 0.13])
    plot(t,Ca,'k')
    hold on
    axis([0 2 0 0.03])
    ylabel('[Ca] ')
    if ( cell_version == 5 )
        xlabel(' Time (sec)')
    end
    mean(VS(ceil(end/2):end))
end

%% Finally add labels to figure panels
figure(1)

annotation('textbox',[0.00 0.99 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.00 0.8 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.00 0.61 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.00 0.42 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','D')
annotation('textbox',[0.00 0.23 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','E')
