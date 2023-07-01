% Figure_7_6.m
% This code is an extended version of nullcline_bistable.m. It includes the
% set of arrows for the trajectories in the phase-plane.
% It also adds an oscillating network.
%
% The fixed points for two coupled equations are shown as the intersection
% of two nullclines, each nullcline representing the fixed points of one
% equation if the other variable were held constant. The equations are:
%
% tau.dr1/dr = -r1 + W11.r1 + W21.r2 - theta1 + I01
% tau.dr2/dt = -r2 + W22.r2 + W12.r1 - theta2 + I02
%
% The time constant does not impact the nullclines, which are obtained by
% setting the left-hand side of each equation to zero.
%
% Two additional conditions cause the nullclines to have a "Z" shape,
% rather than each be a perfect straight line. These conditions are;
% 1) r1 and r2 cannot be less than 0.
% 2) r1 and r2 cannot be more than rmax.
%
% This code is used to produce Figures 7.5 and 7.6 in the book:
% An Introductory Course in Computational Neuroscience
% by Paul Miller (Brandeis University, 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = 0.010;
I01 = 0;                    % Default is no applied current
I02 = 0;                    % Default is no applied current
rmax = 60;                  % Maximum firing rate of each unit
clf
% part 1 is the integrator, parts 2&3 are used for Figure 7.5 in the
% textbook, parts 4&5 are used for Figure 7.6 in the textbook
for part = 1:5
    
    switch part
        
        case 1
            % see lines 10-11 for meaning of parameters
            % connection strengths are W
            W11 = 0.975;
            W12 = -0.025;
            W22 = 0.975;
            W21 = -0.025;
            % thresholds are theta, negative values are spontaneous rates
            theta1 = -0.5;
            theta2 = -0.5;
        case 2
            % see lines 10-11 for meaning of parameters
            % connection strengths are W
            W11 = 1.05;
            W12 = -0.1;
            W22 = 1.05;
            W21 = -0.1;
            % thresholds are theta, negative values are spontaneous rates
            theta1 = 2;
            theta2 = 2;
        case 3
            % see lines 10-11 for meaning of parameters
            % connection strengths are W
            W11 = 2.25;
            W12 = 1.5;
            W21 = -2.25;
            W22 = 0;
            % thresholds are theta, negative values are spontaneous rates
            theta1 = -20;
            theta2 = 20;
        case 4              % part 4 produced Figure 7.6A in the textbook
            % see lines 10-11 for meaning of parameters
            % connection strengths are W
            W11 = 2.25;
            W12 = 1.5;
            W21 = -2.25;
            W22 = -1;
            % thresholds are theta, negative values are spontaneous rates
            theta1 = -20;
            theta2 = -15;
        case 5              % part 5 produced Figure 7.6B in the textbook
            % see lines 10-11 for meaning of parameters
            % connection strengths are W
            W11 = 2.25;
            W12 = 1.5;
            W21 = -2.25;
            W22 = -1;
            % thresholds are theta, negative values are spontaneous rates
            theta1 = -20;
            theta2 = -15;
            I02 = -15;
    end
    
    % set up appropriate figure for each network
    if ( part == 1 )
        figure(part)
        clf
    else
        if ( part < 4 )
            figure(2)
            subplot('Position',[0.1+0.5*(part-2) 0.22 0.36 0.72])
        else
            figure(3)
            subplot('Position',[0.1+0.5*(part-4) 0.62 0.36 0.36])
        end
    end
    % r2-Nullcline, plot of r2 as a function of r1 with dr2/dt = 0
    r1vals = -100:0.1:100;
    r2null = (W12*r1vals - theta2+I02)/(1-W22);
    
    % remove values where rates are < 0 or > rmax
    allowed_vals = find((r2null >= 0 ).*( r2null <= rmax ));
    plot(r1vals(allowed_vals),r2null(allowed_vals),'k:');
    hold on
    
    % r1-Nullcline, plot of r1 as a function of r2 with dr1/dt = 0
    r2vals = -100:0.1:100;
    r1null = (W21*r2vals - theta1+I01)/(1-W11);
    
    % remove values where rates are < 0 or > rmax
    allowed_vals = find((r1null >= 0 ).*( r1null <= rmax ));
    plot(r1null(allowed_vals),r2vals(allowed_vals),'k');
    hold on
    
    % legend('dr_{2}/dt = 0', 'dr_{1}/dt = 0')
    
    % Now remainder of r2 nullcline:
    % vals of r1 where r2 = 0 produces negative dr2/dt
    allowed_vals = find(W12*r1vals - theta2 +I02 < 0 );
    plot(r1vals(allowed_vals),zeros(1,length(allowed_vals)),'k:');
    hold on
    
    % Now vals of r1 where r2 = rmax produces positive dr2/dt
    allowed_vals = find( (W22-1)*rmax + W12*r1vals - theta2 +I02> 0 );
    plot(r1vals(allowed_vals),rmax*ones(1,length(allowed_vals)),'k:');
    hold on
    
    % Now remainder of r1 nullcline:
    % vals of r2 where r1 = 0 produces negative dr1/dt
    allowed_vals = find(W21*r2vals - theta1 + I01 < 0 );
    plot(zeros(1,length(allowed_vals)),r2vals(allowed_vals),'k');
    hold on
    
    % Now vals of r2 where r1 = rmax produces positive dr1/dt
    allowed_vals = find( (W11-1)*rmax + W21*r2vals - theta1 + I01 > 0 );
    plot(rmax*ones(1,length(allowed_vals)),r2vals(allowed_vals),'k');
    hold on
    
    xvals = 0:5:60;
    yvals = 0:5:60;
    [xgrid, ygrid] = meshgrid(xvals,yvals);     % Grid of (x,y) values
    Nvals = length(xgrid(:));                   % Number of points in grid
    W = [W11 W21; W12 W22];                     % Connection matrix
    Ithresh = [theta1; theta2];                 % Threshold of each cell
    
    % The vector [dr(1)/dt; dr(2)/dt] is calculated at each grid point
    drdtvals = (W-eye(2))*[xgrid(:)' ; ygrid(:)'] - Ithresh*ones(1,Nvals);
    drdtvals = drdtvals/tau;    % Divide by time constant for rate of change
    drdtvals = drdtvals';       % We need two columns, not two rows, for the plot
    
    % The next section enforces that dr/dt can not be positive at maximum
    % rates nor negative when rate is zero. These lines force dr/dt to be
    % zero if the dynamical equation would send the rates out of bounds.
    xzeros = find(xgrid(:) == 0 );
    xmax = find(xgrid(:) == rmax );
    yzeros = find(ygrid(:) == 0 );
    ymax = find(ygrid(:) == rmax );
    drdtvals(xzeros,1) = max(drdtvals(xzeros,1),0);
    drdtvals(xmax,1) = min(drdtvals(xmax,1),0);
    drdtvals(yzeros,2) = max(drdtvals(yzeros,2),0);
    drdtvals(ymax,2) = min(drdtvals(ymax,2),0);
    
    % The command "quiver" produces arrows at the grid points specified by
    % each entry in the first two vectors with the x- and y-components of
    % the arrow length specified by each entry of the 3rd and 4th vectors,
    % respectively.
    quiver(xgrid(:), ygrid(:), drdtvals(:,1), drdtvals(:,2) ,'k'); % Plot as arrows
    
    axis([-2 62 -2 62])
    xlabel('r_{1} (Hz)')
    ylabel('r_{2} (Hz)')
    
    if ( part == 2 )
        plot(60,0,'.k','MarkerSize',40)
        plot(0,0,'.k','MarkerSize',40)
        plot(0,60,'.k','MarkerSize',40)
        plot(0,40,'ok','MarkerSize',10)
        plot(40,0,'ok','MarkerSize',10)
    end
    if ( part == 3 )
        plot(30.5885,25.8824,'ok','MarkerSize',10)
    end
    if ( part == 4 )
        plot(7.15,12.85,'.k','MarkerSize',40)
    end
    if ( part == 5 )
        plot(45.72,34.28,'.k','MarkerSize',40)
    end
    
    %% Set up the time vector
    dt = 0.0001;                % time-step
    tmax = 1;                   % maximum time to simulate
    tvec = 0:dt:tmax;           % vector of time-points
    Nt = length(tvec);          % number of time-points
    Nunits = 2;                 
    Iapp = [I01*ones(1,Nt);I02*ones(1,Nt)];
    
    r = zeros(Nunits,Nt);   % array of rate of each cell as a function of time
    %% Now simulate through time
    for i=2:Nt
        I = W*r(:,i-1) + Iapp(:,i-1);                   % total current to each unit
        newr = r(:,i-1) + dt/tau*(I-Ithresh-r(:,i-1));  % Euler-Mayamara update of rates
        r(:,i) = max(newr,0);                           % rates are not negative
        r(:,i) = min(r(:,i),rmax);                      % rates can not be above rmax
    end
    %% Now plot figures rate vs time for each color-coded unit in figure 1,
    %  with a new panel for each trial updating every 1/10 of the
    %  simulation
    
    if ( part < 4 )
        figure(part+10)
        clf
        plot(tvec(1:i),r(1,1:i),'k')
        hold on
        plot(tvec(1:i),r(2,1:i),'r')
    else
        figure(3)
        subplot('Position',[0.1+0.5*(part-4) 0.12 0.36 0.36])
        plot(tvec(1:i),r(1,1:i),'k')
        hold on
        plot(tvec(1:i),r(2,1:i),'k:')
        axis([0 tmax 0 rmax])
        xlabel('Time (sec)')
        ylabel('Rate (Hz)')
        legend('r_{1}', 'r_{2}')
    end
    
end

% Label all panels for Figure 7.6 in the textbook
figure(3)
annotation('textbox',[0.0 0.98 0.02 0.02],'String','C','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.5 0.98 0.02 0.02],'String','D','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.0 0.48 0.02 0.02],'String','E','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.5 0.48 0.02 0.02],'String','F','LineStyle','none','FontSize',16,'FontWeight','bold')
