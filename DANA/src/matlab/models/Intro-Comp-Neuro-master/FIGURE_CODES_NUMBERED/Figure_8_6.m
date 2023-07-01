% Figure_8_6.m
% 
%  Models a recurrent circuit with 4 excitatory cells and 1 inhibitory
%  cell, with recurrent excitatory connections undergoing synaptic
%  plasticity via the triplet rule for STDP.
%  see Pfister and Gerstner, J Neurosci 26:9673 (2006)
%
% The circuit undergoes plasticity with two types of input. 
% In one case, the cells are excited sequentially to produce sequence
% recall after learning.
% In the second case, pairs of excitatory cells are associated by co-stimulation.
% see Clopath et al, Nat Neurosci 13:344 (2010)
%
%  This code was used to produce Figure 8.6 in the book:
%  An Introductory Course in Computational Neuroscience,
%  by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

tmax = 5;               % time for a run (sec)
dt = 0.0001;            % timestep for integration (sec)
t = 0:dt:tmax;          % time vector (time in sec)

%% Set up the inputs to the neuron
NE = 4;                 % No. of excitatory cells
NI = 1;                 % No. of inhibitory cells
N_Groupinputs = 20;           % total number of input cells
Ninputs = NE*N_Groupinputs;

phase_offset = zeros(1,Ninputs);  % phase of cycle for each input, default 0

freq = 5;              % Default frequency of oscillating input (Hz)

tausyn = 0.002;         % synaptic timeconstant (sec)
itau = ceil(tausyn/dt); % number of timesteps in tau

rapp = zeros(length(t),Ninputs);  % input spike rate as function of time for each of two groups
rmax = 200;                  % peak input spike rate

%% Set up initial synaptic input strengths
Win = 4e-9;
W0 = 0.5e-9;                 % initial average synaptic strengths
Wrand = 0e-9;

%% Set up rules for changes in synaptic strengths
dW_stdpp_3 = 0.03e-9;   % maximum increase in strength for each triplet (S)
dW_stdpp_2 = 0.02e-9;   % maximum increase in strength for each pairing (S)
dW_stdpm = 0.060e-9;    % maximum decrease in strength for each pairing (S)
tau_post_2 = 0.030;     % time window of 2 spikes for depression (sec)
tau_pre_2 = 0.020;      % time window of 2 spikes for potentiation (sec)
tau_post_3 = 0.050;     % time window of 2 post-spikes for potentiation (sec)
tau_s = 0.050;          % excitatory synaptic time constant
tau_sI = 0.005;         % inhibitory synaptic time constant 
Wmax = 15*W0;           % maximum allowed synaptic strength


%% Parameters for the LIF neuron being simulated
Vth = -0.050;           % threshold voltage (V)
Vpeak = 0.050;          % Peak deflection added for visualization (V)
E_L = -0.070;           % leak voltage (V)
E_synE = 0.000;         % Reversal potential of synaptic inputs (V)
E_synI = -0.080;        % Reversal potential of synaptic inputs (V)
Vreset = -0.090;        % reset voltage (V)
Cm = 200e-12;           % specific membrane capacitance (F)
G_L = 5e-9;             % specific leak conductance (S)

%% Set up matrices for storing initial and final voltage traces and 
%  final connectivity matrix
V_init = zeros(NE+NI,length(t),2);
V_final = zeros(NE+NI,length(t),2);
WEE_final = zeros(NE,NE,2);

%% First do a sequence of four inputs, then 2 groups of 2 inputs.
for part = 1:2
    
    WEE = W0*ones(NE) + Wrand*randn(NE);
    WEE = WEE - diag(diag(WEE));  % Remove diagonal entries
    WEE = min(WEE,Wmax);        % ensure W is no higher than Wmax
    WEE = max(WEE,0);           % and no lower than zero
    WEI = 8*W0*ones(NE,NI);
    WIE = 25*W0*ones(NI,NE);
    
    % in the following section "thresh" determines the length of each input
    % pulse (the fraction of a cycle at which the sine is above thresh is
    % the fraction of a cycle with input).
    % phase_offset determines whether different input cells spike at the
    % same time or at different phase-relationships to each other.
    % the input cells are split into equal-sized groups, one group for each
    % excitatory LIF neuron being simulated.
    if ( part == 1 ) % Sequence of four inputs
        Ntrials = 6;
        for i = 1:NE;
            phase_offset((i-1)*N_Groupinputs+1:i*N_Groupinputs) = i*2*pi/NE;
            thresh = 0.98;
        end
    else
        if ( part == 2 ) 
            Ntrials = 2;
            phase_offset(1:Ninputs/2) = 0;
            phase_offset(Ninputs/2+1:Ninputs) = pi;
            thresh = 0.9;
        end
    end
    %% Begin the loop through all trials with plasticity implemented by the
    %  the continuous update method
    
    for trial = 1:Ntrials   % begin trials
        
        %  Inputs are excitatory conductances due to periodic pulses of 
        %  Poisson inputs. 
        for i = 1:Ninputs
            rapp(:,i) = rmax*(sin(2*pi*freq*t - phase_offset(i) ) > thresh );   % input rate for group 1
        end
        
        % In the final trial just simulate one input group for one cycle
        if ( ( trial == Ntrials ) && ( Ntrials > 1 ) )
            rapp(:,N_Groupinputs+1:end) = 0;
            rapp(round(1/(freq*dt))+1:end,:) = 0;
        end
        spikesin = zeros(Ninputs,length(t));    % series of spikes for each input
        spikesout = zeros(NE,length(t));        % series of spikes generated by LIF neuron
        spikesoutI = zeros(NI,length(t));       % series of spikes generated by LIF neuron
        G_synE = zeros(NE,length(t));           % input current to LIF neuron
        G_synI = zeros(NI,length(t));           % input current to LIF neuron
        
        for cell = 1:NE
            for i = (cell-1)*N_Groupinputs+1:cell*N_Groupinputs                     % first group of inputs
                for j = 1:length(t)
                    if rapp(j,i)*dt > rand(1)       % random Poisson process at rate Iapp
                        spikesin(i,j) = 1;          % to generate spikes at inputs
                        
                        kmax = min(j+5*itau,length(t));     % how long to add current
                        for k = j+1:kmax                    % from the spike
                            G_synE(cell,k) = G_synE(cell,k) + ...         % add current from input i
                                Win*exp(-dt*(k-j)/tausyn);
                        end
                    end
                end
            end
        end
        
        G_synI = ones(NI,1)*mean(G_synE);
        % Completed generation of set of input current from input spike trains
        %% Now see how the LIF neuron responds to these inputs
        
        V = zeros(NE,length(t));    % membrane potential of LIF neuron
        V(:,1) = E_L;               % initialize at leakage potential
        VI = zeros(NI,length(t));   % membrane potential of LIF neuron
        VI(:,1) = E_L;              % initialize at leakage potential
        
        post_2 = zeros(NE,1);       % depends on postsynaptic spike for LTD
        post_3 = zeros(NE,1);       % depends on prior postsynaptic spike for LTP
        pre_2 = zeros(NE,1);        % depends on presynaptic spike for LTP

        % synaptic conductances indexed by presynaptic neuron
        sE = zeros(NE,1);           % excitatory synaptic conductances
        sI = zeros(NI,1);           % inhibitory synaptic conductances
        
        for i = 2:length(t)     % integrate through time
            % all synaptic conductances and spike-window terms decay
            % exponentially between updates:
            post_2 = post_2*exp(-dt/tau_post_2);
            post_3 = post_3*exp(-dt/tau_post_3);
            pre_2 = pre_2*exp(-dt/tau_pre_2);
            sE = sE*exp(-dt/tau_s);
            sI = sI*exp(-dt/tau_sI);
            
            %% First update membrane potential of E cells
            GtotE = G_L + G_synE(:,i-1) + WEE'*sE + WIE'*sI;
            Vinf = (E_L*G_L + E_synE*(G_synE(:,i-1)+ WEE'*sE) ...
                + E_synI*WIE'*sI)./GtotE;           % calculate steady state voltage
            V(:,i) = Vinf + (V(:,i-1)-Vinf) ...
                .*exp(-dt*GtotE/Cm);  % and integrate towards it
            
            cellsfired = find(V(:,i) >= Vth);     % cells whose voltage passes threshold
            
            spikesout(cellsfired,i) = 1;           % generate spike in LIF neuron
            V(cellsfired,i) = Vreset;              % and reset the membrane potential
            V(cellsfired,i-1) = Vpeak;
            
            % Now update synaptic strengths with the continuous update rule
            % First when there is a postsynaptic spike, possuble
            % potentiation
            WEE(:,cellsfired) = WEE(:,cellsfired) + pre_2(:) ...
                *(dW_stdpp_2*ones(1,length(cellsfired)) ...
                + dW_stdpp_3*post_3(cellsfired)');
            % Then when there is a presynaptic spike, possible depression
            WEE(cellsfired,:) = WEE(cellsfired,:) - dW_stdpm ...
                *(ones(length(cellsfired),1)*post_2(:)');
            WEE = min(WEE,Wmax);                        % now ensure no W is above maximum allowed
            WEE = max(WEE,0);                           % and no W is less than zero
            WEE = WEE.*(1-eye(NE));                     % remove self-connections
            
            % Update all terms due to the cells that fired a spike
            post_2(cellsfired) = 1 + post_2(cellsfired);
            post_3(cellsfired) = 1 + post_3(cellsfired);
            pre_2(cellsfired) = 1 + pre_2(cellsfired);
            sE(cellsfired) = 1;
            
            %% Then update membrabe potential of I cells
            GtotI = G_L + G_synI(:,i-1) + WEI'*sE;
            VinfI = (E_L*G_L + E_synE*(G_synI(:,i-1) + WEI'*sE) ) ...
                ./GtotI;           % calculate steady state voltage
            VI(:,i) = VinfI + (VI(:,i-1)-VinfI) ...
                .*exp(-dt*GtotI/Cm);  % and integrate towards it
            
            cellsfiredI = find(VI(:,i) >= Vth);     % cells whose voltage passes threshold
            
            spikesoutI(cellsfiredI,i) = 1;           % generate spike in LIF neuron
            VI(cellsfiredI,i) = Vreset;              % and reset the membrane potential
            VI(cellsfiredI,i-1) = Vpeak;
            sI(cellsfiredI) = 1;
            
        end
        
        % Plot matrix of connections and voltage trace each trial
        if ( mod(trial,1) == 0 )
            WEE./W0
            
            figure(1)
            imagesc(WEE)
            caxis([0 Wmax])
            drawnow
            figure(2)
            
            if ( trial < Ntrials )
                plot(t,V)
                axis([0 1 Vreset Vth])
                drawnow
            end
        end
        
        if ( trial == 1 )
            figure(100*part)
            subplot(3,1,1)
            for cell = 1:NE
                plot(t,V(cell,:)-0.2*cell,'k')
                hold on
            end
            axis([0 0.5 -0.1-NE*0.2 0.05-0.2])
            xlabel('Time (sec)')
            ylabel('Membrane potential')
            set(gca,'YTick',[])
            V_init(1:NE,:,part) = V;
            V_init(NE+1:end,:,part) = VI;
        end
        if ( trial == Ntrials )
            figure(100*part)
            subplot(3,1,2)
            for cell = 1:NE
                plot(t,V(cell,:)-0.2*cell,'k')
                hold on
            end
            axis([0 0.5 -0.1-NE*0.2 0.05-0.2])
            xlabel('Time (sec)')
            ylabel('Membrane potential')
            set(gca,'YTick',[])
            V_final(1:NE,:,part) = V;
            V_final(NE+1:end,:,part) = VI;
            
        end
        drawnow
        
    end
    figure(100*part)
    subplot(3,1,3)
    imagesc(WEE)
    colormap gray
    WEE_final(:,:,part) = WEE;
    
    figure(3)
    
    plot(t,V)
    axis([0 1 Vreset Vth])
    drawnow
    
end

%% Now plot the figure for the textbook (Figure 8.6)
clf
for part = 1:2
    figure(300)
    
    subplot('Position',[0.08 (1.02-0.52*part)+0.32 0.35 0.145])
    for cell = 1:NE
        plot(t,squeeze(V_init(cell,:,part))-0.2*cell,'k')
        hold on
    end
    axis([0 0.5 -0.1-NE*0.2 0.05-0.2])
    ylabel('V_{m}')
    set(gca,'YTick',[])
    if (part == 1 )
        title('Early Sequence Training')
    end
    if (part == 2 )
        title('Early Association Training')
    end
    
    subplot('Position',[0.08 (1.02-0.52*part)+0.09 0.35 0.145])
    for cell = 1:NE
        plot(t,squeeze(V_final(cell,:,part))-0.2*cell,'k')
        hold on
    end
    axis([0 0.5 -0.1-NE*0.2 0.05-0.2])
    xlabel('Time (sec)')
    ylabel('V_{m}')
    set(gca,'YTick',[])
    title('Test after Training')
    
    subplot('Position',[0.55 (1.01-0.51*part)+0.09 0.4 0.37])
    imagesc(1e9*squeeze(WEE_final(:,:,part)))
    colormap gray
    caxis([0 1e9*Wmax])
    colorbar
    title('Final Connection Strengths (nS)')
    ylabel('Presynaptic Cell')
    xlabel('Postsynaptic Cell')
    set(gca,'XTick',[1:4])
    set(gca,'YTick',[1:4])
    
    %% Finally label the figures A1-A3, B1-B3  at the appropriate points
    if ( part == 1 )
        annotation('textbox',[0.00 0.98 0.02 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','A1')
        annotation('textbox',[0.00 0.76 0.02 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','A2')
        annotation('textbox',[0.48 0.98 0.02 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','A3')
    end
    if ( part == 2 )
        annotation('textbox',[0.00 0.47 0.02 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','B1')
        annotation('textbox',[0.00 0.25 0.02 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','B2')
        annotation('textbox',[0.48 0.47 0.02 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','B3')
    end
    
    
end

save STDP_recurrent
