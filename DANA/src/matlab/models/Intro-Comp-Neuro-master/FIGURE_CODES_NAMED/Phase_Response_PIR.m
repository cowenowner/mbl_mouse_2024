% Phase_Response_PIR.m
% This model contains a T-type Calcium current to generate a
% post-inhibitory rebound as a models of thalamic relay cells.
%
% This code is used to produce Figure 5.11D of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dt = 1e-6;
Ibase = 0;
Istep = 0;
%
tmax = 1;
tvec = 0:dt:tmax;

E_L = -0.070;       % leak reversal potential
E_Na = 0.055;       % reversal for sodium channels
E_K = -0.090;       % reversal for potassium channels
E_Ca = 0.120;        % reversal potential for Ca current

G_L = 10e-9;        % leak conductance in Siemens
G_Na = 1e-6;      % sodium conductance in Siemens
G_K = 1e-6;       % potassium conductance in Siemens
G_CaT = 5e-6;    % T-type calcium conductance in Siemens

Cm = 0.1e-9;        % membrane capacitance in Farads

V = zeros(size(tvec));
V(1) = E_L;

m = zeros(size(tvec));
h = zeros(size(tvec));
n = zeros(size(tvec));
mca = zeros(size(tvec));
hca = zeros(size(tvec));
%% Now add current pulses at different points on the cycle and analyze the
%  change in response due to the pulse.

Ipulse_amp = 10e-12;

Npulses = 200;
Ipulse = 0;
shift = zeros(1,Npulses);

for trial = 0:Npulses;  % INitially trial 0 is used to get the baseline
    
    if ( trial == 1 )
        i_tpulse = i_startphase + floor((i_stopphase-i_startphase)*[1:Npulses]/Npulses);
    end
    
    I = Ibase*ones(size(tvec));
    if ( trial > 0 )
        I(i_tpulse(trial):i_tpulse(trial)+round(0.005/dt)) = ...
            I(i_tpulse(trial):i_tpulse(trial)+round(0.005/dt)) + Ipulse_amp;
    end
    for i = 2:length(tvec); % now see how things change through time
        
        Vm = V(i-1); % For ease of repeated use in the equations
        
        % alpha_m is sodium activation rate constant
        if ( Vm == -0.035 )                            % to prevent division by zero
            alpha_m = 1e3;
        else
            alpha_m = 1e5*(Vm+0.035)./(1-exp(-100*(Vm+0.035)));
        end
        beta_m = 4000*exp(-(Vm+0.060)./0.018);      % sodium deactivation rate
        
        alpha_h = 350*exp(-50*(Vm+0.058));          % sodium inactivation rate
        beta_h = 5000./(1+exp(-100*(Vm+0.028)));    % sodium deinactivation rate
        
        if ( Vm == -0.034 )                            % prevent division by zero
            alpha_n = 500;                          % potassium activation rate constant
        else
            alpha_n = 5e4*(Vm+0.034)./(1-exp(-100*(Vm+0.034)));
        end
        beta_n = 625*exp(-12.5*(Vm+0.044));     % potassium inactivation rate
        
        % From the alpha and beta for each gating variable we find the steady
        % state values (_inf) and the time constants (tau_) for each m,h and n.
        
        m_inf = alpha_m./(alpha_m+beta_m);      % sodium activation variable
        
        tau_h = 1./(alpha_h+beta_h);            % sodium inactivation time constant
        h_inf = alpha_h./(alpha_h+beta_h);      % sodium inactivation variable
        
        tau_n = 1./(alpha_n+beta_n);            % potassium activation time constant
        n_inf = alpha_n./(alpha_n+beta_n);      % potassium activation variable
        
        mca_inf = 1./(1+exp(-(Vm+0.052)/0.0074));   % Ca T-current activation variable
        
        hca_inf = 1./(1+exp(500*(Vm+0.076)));       % Ca T-current inactivation variable
        
        % Below the Ca T-current time constant has two terms, one for Vm<-80mV,
        % the other for Vm > -80mV. The expressions in parenthesis, such as
        % (Vm<-0.080) yield 1 or 0 depending on the value of Vm, so either the
        % first line is used or the second term line used.
        tau_hca = 1e-3*exp(15*(Vm+0.467)).*(Vm<-0.080) ...
            + 1e-3*(28+exp(-(Vm+0.022)/0.0105)).*(Vm>=-0.080);
        
        m(i) = m_inf;    % Update m, assuming time constant is neglible.
        
        h(i) = h(i-1) + dt*(h_inf - h(i-1))/tau_h;    % Update h
        
        n(i) = n(i-1) + dt*(n_inf - n(i-1))/tau_n;    % Update n
        
        mca(i) = mca_inf;    % Update mca, assuming time constant is negligible
        hca(i) = hca(i-1) + dt*(hca_inf - hca(i-1))/tau_hca;    % Update hca
        
        % Now calculate instantaneous values of the conductances
        G_Na_now = G_Na*m(i)*m(i)*m(i)*h(i);
        
        G_K_now = G_K*n(i)*n(i)*n(i)*n(i);
        
        G_CaT_now = G_CaT*mca(i)*mca(i)*hca(i);
        
        % Now update the membrane potential
        V(i) = Vm + dt*(G_L*(E_L-Vm) + G_Na_now*(E_Na-Vm) ...
            +G_K_now*(E_K-Vm) + G_CaT_now*(E_Ca-Vm) + I(i-1) )/Cm;
    end
    
    %% Detect the initiation time of individual bursts
    inburst = 0;
    tburst = [];
    Nbursts = 0;
    for i = 1:length(tvec)
        if ( inburst == 0 ) && ( V(i) > -0.010 )
            inburst = 1;
            tburst(end+1) = tvec(i);
            Nbursts = Nbursts + 1;
        end
        if (inburst == 1 ) && ( V(i) < -0.065 )
            inburst = 0;
        end
    end
    %% Now decide where a phase of "zero" corresponds to in the oscillation
    %  and calculate the time where this occurs -- time should be well after any
    %  initial transients.
    if ( trial == 0 )
        i_startphase = round(tburst(3)/dt);
        i_stopphase = round(tburst(4)/dt);
        period = (tburst(6) - tburst(3))/3;
    else
        shift(trial) = 2*pi*(3 - (tburst(6) - tburst(3))/period );
        if ( shift(trial) < -pi )
            shift(trial) = shift(trial) + 2*pi;
        end
        if ( shift(trial) > pi )
            shift(trial) = shift(trial) - 2*pi;
        end
        
        data = [shift(trial) Nbursts tburst(4)]
    end
    
    %% Now plot the graphs
    if ( mod(trial,10) == 0 )
        set(0,'DefaultLineLineWidth',2,...
            'DefaultLineMarkerSize',8, ...
            'DefaultAxesLineWidth',2, ...
            'DefaultAxesFontSize',14,...
            'DefaultAxesFontWeight','Bold');
        
        figure(trial+1)
        plot(tvec,V,'k')
        axis([tburst(3)-0.05 tburst(4)+0.05 -0.075 0.005])
        xlabel('Time (sec)')
        ylabel('Membrane potential (V)')
        drawnow
    
    end
end

figure()
plot([1:Npulses]*2*pi/Npulses,shift,'k')
hold on
plot([1:Npulses]*2*pi/Npulses,[1:Npulses]*0,'k:')
xlabel('Phase of pulse')
ylabel('Phase shift')
axis([0 2*pi -0.1 0.2])
set(gca,'YTick',[-0.1 0 0.1 0.2])
set(gca,'XTick',[0 pi 2*pi])
set(gca,'XTickLabel',{'0' '\pi' '2\pi'})


