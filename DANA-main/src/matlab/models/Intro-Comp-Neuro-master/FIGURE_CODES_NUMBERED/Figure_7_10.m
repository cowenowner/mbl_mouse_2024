% Figure_7_10.m
% Fitzhugh-Nagumo model with only time-dependence.
% A 2-variable model of a spiking neuron, with a fast variable "V"
% and a slow adaptation variable, "w" representing inactivation and
% K-channels.
%
% This code is used to produce Figures 7.9 and 7.10 in the book:
% An Introductory Course in Computational Neuroscience,
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clf

% time for integration
dt = 0.001;         % time-step in ms
tmax = 100;         % max time in ms
t = 0:dt:tmax;      % vector of time points

% initialize V and w to zero
V= zeros(size(t));
w = zeros(size(t));

%% Basic parameters for the model (see notes)
V1 = -0.07;         % approximate fixed point without current
V2 = 0.050;         % approximate peak of oscillating voltage
V0 = -0.050;        % level above which adaptation becomes positive

beta = 8;           % sets importance of V self-feedback over w
tau_V = 0.01;       % time constant (ms) for rapid change in V
tau_w = 200;        % time constant (ms) for change in w


%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf

%% Loop through three different values of applied current
Iappvec = [0.009 0.01 0.1];
Nloops = length(Iappvec);

for loop = 1:Nloops
    Iapp = Iappvec(loop);
    % Initial conditions
    Vinit = -0.045;
    winit = 0.008;
    if ( loop == 3 )
        winit = 0.1;
    end
    % Now start time integration for V and w.
    V(1) = Vinit;
    w(1) = winit;
    for i = 2:length(t)
        V(i) = V(i-1) + dt*(-beta*V(i-1)*(V(i-1)-V1)*(V(i-1)-V2)-w(i-1)+Iapp)/tau_V;
        w(i) = w(i-1) + dt*(V(i-1)-V0-w(i-1) )/tau_w;
    end
    
    %% Now plot in the V-w phase plane the nullclines,
    %  which means the curve dV/dt = 0 then the line
    % dw/dt = 0.
    subplot('Position',[-0.23+0.33*loop 0.74 0.21 0.22])
    Vplot = -0.2:0.0001:0.2;
    plot(1000*Vplot,-1000*beta*Vplot.*(Vplot-V1).*(Vplot-V2)+1000*Iapp,'k', ...
        'LineWidth',2);
    hold on
    plot(1000*Vplot,1000*(Vplot-V0),'k:','LineWidth',2);
    
    % Next how V and w covary during the simulation
    plot(1000*V(1:round(end/3)),1000*w(1:round(end/3)),'k--','LineWidth',2)
    
    if ( loop ~= 2 )
        plot(1000*V(end),1000*w(end),'.k','MarkerSize',30)
    end
    
    if ( loop < 3 )
        axis([-100 100 4 12])
    else
        axis([-100 100 96 104])
    end
    title(strcat(['I^{(app)} = ',num2str(Iapp)]));
    ylabel('Adaptation, w (mV)')
    xlabel('Potential, V (mV)')
    
    %% Separate figure (for Fig 7.9) using one value of current
    if ( loop == 2 )
        figure(2)
        Vplot = -0.2:0.0001:0.2;
        plot(1000*Vplot,-1000*beta*Vplot.*(Vplot-V1).*(Vplot-V2)+1000*Iapp,'k', ...
            'LineWidth',2);
        hold on
        plot(1000*Vplot,1000*(Vplot-V0),'k:','LineWidth',2);
        
        % Next how V and w covary during the simulation 
        % (skip initial transient and only one cycle)
        plot(1000*V(round(end/6):round(end/3)),1000*w(round(end/6):round(end/3)),'k--','LineWidth',2)
        ylabel('Adaptation, w (mV)')
        xlabel('Membrane Potential, V (mV)')
        axis([-100 100 8.5 11.5])
        
        figure(1)
    end
    %% Then multiple plots of V and w against time for Fig. 7.10
    subplot('Position',[-0.23+0.33*loop 0.41 0.21 0.22])
    plot(t,V*1000,'k','LineWidth',2)
    ylabel('Potential, V (mV)')
    xlabel('Time (ms)')
    
    axis([0 tmax -90 90])
    
    subplot('Position',[-0.23+0.33*loop 0.08 0.21 0.22])
    plot(t,w*1000,'k','LineWidth',2)
    ylabel('Adaptation, w (mV)')
    if ( loop < 3 )
        axis([0 tmax 4 12])
    else
        axis([0 tmax 95 105])
    end
    xlabel('Time (ms)')
end

% Label all panels in the graph
annotation('textbox',[0.0 0.99 0.02 0.02],'String','A1','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.0 0.65 0.02 0.02],'String','A2','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.0 0.35 0.02 0.02],'String','A3','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.33 0.99 0.02 0.02],'String','B1','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.33 0.65 0.02 0.02],'String','B2','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.33 0.35 0.02 0.02],'String','B3','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.66 0.99 0.02 0.02],'String','C1','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.66 0.65 0.02 0.02],'String','C2','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.66 0.35 0.02 0.02],'String','C3','LineStyle','none','FontSize',16,'FontWeight','bold')
