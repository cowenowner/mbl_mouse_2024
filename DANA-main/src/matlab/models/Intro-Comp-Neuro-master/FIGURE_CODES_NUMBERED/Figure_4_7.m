% Figure_4_7.m
% This code runs through multiple trials of applied current, commencing
% each trial from the final state of the prior trial.
%
% This model is the Hodgkin-Huxley model in new units.
%
% This code is used to produce Figure 4.7 of the textbook
% An Introductory Course in Computational Neuroscience
% by Paul Miller
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dt = 2e-8;          % time-step for integration (ms)


%% Neuron parameters

V_L = -0.060;       % leak reversal potential (V)
E_Na = 0.045;       % reversal for sodium channels (V)
E_K = -0.082;       % reversal for potassium channels (V)
V0 = -0.065;

G_L = 30e-9;        % specific leak conductance (S)
G_Na = 12e-6;       % specific sodium conductance (S)
G_K = 3.6e-6;         % specific potassium conductance (S)

Cm = 100e-12;       % specific membrane capacitance (F)


Imin1 = 0.2e-9;
Imax1 = 1.0e-9;
Imin2 = 0.4e-9;
Imax2 = 1.2e-9;
Istep = 0.01e-9;

Ievec = [Imin1:Istep:Imax1; Imin2:Istep:Imax2];  % series of applied currents (A)

Ntrials = length(Ievec(1,:));    % Number of trials to loop through

Nspikes = zeros(2,Ntrials);
sustained = zeros(2,Ntrials);
rate = zeros(2,Ntrials);

for fig_number = 2:2
    if ( fig_number == 1 )
        Nramps = 1;
        istart = 0.2;        % time applied current starts (sec)
        ilength= 1;        % length of applied current pulse (sec)
        Ibase = 0e-7;
        
    else
        Nramps = 2;
        istart = 0;        % time applied current starts (sec)
        ilength=1;        % length of applied current pulse (sec)
        Ibase = 0e-7;
    end
    tmax=istart+ilength;             % maximum time of simulation (s)
    t=0:dt:tmax;        % time vector
    V=zeros(size(t));           % membrane potential vector
    n=zeros(size(t));       % n: potassium activation gating variable
    m=zeros(size(t));       % m: sodium activation gating variable
    h=zeros(size(t));       % h: sodim inactivation gating variable
    
    for ramp = 1:Nramps
        
        for trial = 1:Ntrials
            if ( fig_number == 1 )
                Ie = Ievec(1,trial)
            else
                if ( ramp == 1 )
                    Ie= Ievec(2,trial);           % New applied current each trial
                else
                    Ie = Ievec(2,Ntrials+1-trial);
                end
            end
            spike = zeros(size(t));     % records spike times
            %% current clamp initialization
            Iapp=Ibase*ones(size(t));   % Applied current, relevant in current-clamp mode
            for i=round(istart/dt)+1:round((istart+ilength)/dt+1) % make non-zero for duration of current pulse
                Iapp(i) = Ie;
            end
            
            
            %% General initialization
            if ( fig_number == 1 )
                V(1) = V_L;
                m(1) = 0;
                h(1) = 0;
                n(1) = 0;
            else
                V(1) = V(end);          % set the inititial value of voltage
                m(1) = m(end);
                h(1) = h(end);
                n(1) = n(end);
            end
            Itot=zeros(size(t));    % in case we want to plot and look at the total current
            I_Na=zeros(size(t));    % record sodium curret
            I_K=zeros(size(t));     % record potassium current
            I_L=zeros(size(t));     % record leak current
            lastspiketime = 0;
            inspike = 0;
            for i = 2:length(t); % now see how things change through time
                
                Vm = V(i-1);          % membrane potential for calculations
                
                % Sodium and potassium gating variables are defined by the
                % voltage-dependent transition rates between states, labeled alpha and
                % beta.
                
                % First, sodium activation rate
                if ( Vm == -0.045 )     % to avoid dividing zero by zero
                    alpha_m = 1e3;      % value calculated analytically
                else
                    alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
                end
                beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
                alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
                beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
                
                if ( Vm == -0.060)      % to avoid dividing by zero
                    alpha_n = 100;      % value calculated analytically
                else;                   % potassium activation rate
                    alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
                end
                beta_n = 125*exp((-Vm-0.070)/0.08);     % potassium deactivation rate
                
                % From the alpha and beta for each gating variable we find the steady
                % state values (_inf) and the time constants (tau_) for each m,h and n.
                
                tau_m = 1/(alpha_m+beta_m);
                m_inf = alpha_m/(alpha_m+beta_m);
                
                tau_h = 1/(alpha_h+beta_h);
                h_inf = alpha_h/(alpha_h+beta_h);
                
                tau_n = 1/(alpha_n+beta_n);
                n_inf = alpha_n/(alpha_n+beta_n);
                
                
                m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
                
                h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
                
                n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
                
                I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-V(i-1)); % total sodium current
                
                I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-V(i-1)); % total potassium current
                
                I_L(i) = G_L*(V_L-V(i-1));    % Leak current is straightforward
                
                Itot(i) = I_L(i)+I_Na(i)+I_K(i)+Iapp(i); % total current is sum of leak + active channels + applied current
                
                V(i) = V(i-1) + Itot(i)*dt/Cm;        % Update the membrane potential, V.
                
                if ( ( inspike == 0 ) && ( V(i) > 0 ) )
                    if ( t(i) > istart )
                        spike(i+1) = 1;
                    end
                    inspike = 1;
                    lastspiketime = t(i);
                else
                    if ( inspike == 1 )
                        if ( V(i) < -0.020 )
                            inspike = 0;
                        end
                    end
                end
                
            end
            
            spiketimes = dt*find(spike)
            Nspikes(ramp,trial) = length(spiketimes);
            if ( Nspikes(ramp,trial) > 1 )
                ISIs = diff(spiketimes);
                if ( lastspiketime + max(ISIs)*1.1 > tmax )
                    sustained(ramp,trial) = 1;
                    rate(ramp,trial) = 1/ISIs(end);
                end
            end
            %         figure(10*ramp+trial)
            %         plot(t,-70-V)
            %         drawnow
        end
        
    end
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    figure(fig_number)
    clf
    plot(Ievec(fig_number,:)*1e9,rate(1,:),'k');
    hold on
    if ( fig_number == 2 )
        hold on
        plot(Ievec(fig_number,end:-1:1)*1e9,rate(2,:),'k-');
    end
    single_spikes1 = find(Nspikes(1,:) == 1);
    plot(Ievec(fig_number,single_spikes1)*1e9,zeros(1,length(single_spikes1)),'.k')
    
    multi_spikes1 = find((Nspikes(1,:) > 1).*(sustained(1,:) == 0));
    plot(Ievec(fig_number,multi_spikes1)*1e9,zeros(1,length(multi_spikes1)),'*k')

    if ( ramp == 2 )
        single_spikes2 = find(Nspikes(2,:) == 1);
        plot(Ievec(fig_number,single_spikes2)*1e9,zeros(1,length(single_spikes2)),'ok')
        multi_spikes2 = find((Nspikes(2,:) > 1).*(sustained(2,:) == 0));
        plot(Ievec(fig_number,multi_spikes2)*1e9,zeros(1,length(multi_spikes2)),'dk')
    end
    
    xlabel('Applied current (nA)')
    ylabel('Firing rate (Hz)')
    
    if ( fig_number == 2 )
        x = [0.74 0.74];
        y = [0.3 0.5];
        annotation('arrow',x, y)
        
        x = [0.34 0.34];
        y = [0.5 0.3];
        annotation('arrow',x, y)
        
        annotation('textbox',[0.43 0.4 0.15 0.1],'String','Bistable Range','LineStyle','none','FontSize',16,'FontWeight','Bold')
    end
end
save('HH_f_I_bi.mat')
