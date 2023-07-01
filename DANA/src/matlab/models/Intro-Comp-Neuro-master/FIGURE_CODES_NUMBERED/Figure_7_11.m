% Figure_7_11.m
%
% This code is an extended version of nullcline_bistable.m. It includes the
% set of arrows for the trajectories in the phase-plane.
% It also adds an oscillating network.
%
% The fixed points for two coupled equations are shown as the intersection
% of two nullclines, each nullcline representing the fixed points of one
% equation if the other variable were held constant. The equations are:
%
% tau.dr1/dr = -r1 + W11.r1 + W21.r2 - theta1
% tau.dr2/dt = -r2 + W22.r2 + W12.r1 - theta2
%
% The time constant does not impact the nullclines, which are obtained by
% setting the left-hand side of each equation to zero.
%
% Two additional conditions cause the nullclines to have a "Z" shape,
% rather than each be a perfect straight line. These conditions are;
% 1) r1 and r2 cannot be less than 0.
% 2) r1 and r2 cannot be more than rmax.
%
% This code is used to produce Figure 7.11 in the book:
% An Introductory Course in Computational Neuroscience
% by Paul Miller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed parameters
tau = 0.010;        % time constant
W11 = 0.95;         % connection strength 1-to-1
W12 = -1;           % connection strength 1-to-2
W21 = -0.5;         % connection strength 2-to-1
W22 = 0;            % connection strength 2-to-2
theta1 = -10;       % threshold of unit 1 (negative means spontaneous rate)
theta2 = -40;       % threshold of unit 2 (negative means spontaneous rate)
rmax = 60;          % maximum rate of each unit

%% Loop with two different initial conditions
for part = 1:2
    
    switch part
        case 1
            r1init = 15.5;      % Unit 1 initial rate
            r2init = 0;         % Unit 2 initial rate
        case 2
            r1init = 38.5;      % Unit 1 initial rate
            r2init = 60;        % Unit 2 initial rate
    end
    
    subplot('Position',[0.62-0.5*(part-1) 0.62 0.36 0.36])
    
    % r2-Nullcline, plot of r2 as a function of r1 with dr2/dt = 0
    r1vals = -100:0.1:100;
    r2null = (W12*r1vals - theta2)/(1-W22);
    
    % remove values where rates are < 0 or > rmax
    allowed_vals = find((r2null >= 0 ).*( r2null <= rmax ));
    plot(r1vals(allowed_vals),r2null(allowed_vals),'k:');
    hold on
    
    % r1-Nullcline, plot of r1 as a function of r2 with dr1/dt = 0
    r2vals = -100:0.1:100;
    r1null = (W21*r2vals - theta1)/(1-W11);
    
    % remove values where rates are < 0 or > rmax
    allowed_vals = find((r1null >= 0 ).*( r1null <= rmax ));
    plot(r1null(allowed_vals),r2vals(allowed_vals),'k');
    hold on
    
    % Now remainder of r2 nullcline:
    % vals of r1 where r2 = 0 produces negative dr2/dt
    allowed_vals = find(W12*r1vals - theta2 < 0 );
    plot(r1vals(allowed_vals),zeros(1,length(allowed_vals)),'k:');
    hold on
    
    % Now vals of r1 where r2 = rmax produces positive dr2/dt
    allowed_vals = find( (W22-1)*rmax + W12*r1vals - theta2 > 0 );
    plot(r1vals(allowed_vals),rmax*ones(1,length(allowed_vals)),'k:');
    hold on
    
    % Now remainder of r1 nullcline:
    % vals of r2 where r1 = 0 produces negative dr1/dt
    allowed_vals = find(W21*r2vals - theta1 < 0 );
    plot(zeros(1,length(allowed_vals)),r2vals(allowed_vals),'k');
    hold on
    
    % Now vals of r2 where r1 = rmax produces positive dr1/dt
    allowed_vals = find( (W11-1)*rmax + W21*r2vals - theta1 > 0 );
    plot(rmax*ones(1,length(allowed_vals)),r2vals(allowed_vals),'k');
    hold on
    
%% Grid of values to produce the vector field, which will show direction
%  of rates of change for different pairs of (r1,r2)
    xvals = 0:5:60;
    yvals = 0:5:60;
    [xgrid, ygrid] = meshgrid(xvals,yvals);     % Grid of (x,y) values
    Nvals = length(xgrid(:));                   % Number of points in grid
    W = [W11 W21; W12 W22];                     % Connection matrix
    Ithresh = [theta1; theta2];                   % Threshold of each cell
    
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
    
    plot(0,40,'.k','MarkerSize',40)
    plot(60,0,'.k','MarkerSize',40)
    plot(22.2,17.8,'ok','MarkerSize',10)
    
    %% Set up the time vector
    dt = 0.0001;
    tmax = 0.25;
    tvec = 0:dt:tmax;
    Nt = length(tvec);
    Nunits = 2;
    
    r = zeros(Nunits,Nt);   % array of rate of each cell as a function of time
    r(:,1) = [r1init; r2init];
    %% Now simulate through time
    for i=2:Nt
        I = W*r(:,i-1);                   % total current to each unit
        newr = r(:,i-1) + dt/tau*(I-Ithresh-r(:,i-1));  % Euler-Mayamara update of rates
        r(:,i) = max(newr,0);                           % rates are not negative
        r(:,i) = min(r(:,i),rmax);                      % rates can not be above rmax
    end
    %% Now plot figures rate vs time for each color-coded unit in figure 1,
    %  with a new panel for each trial updating every 1/10 of the
    %  simulation
    
    subplot('Position',[0.62-0.5*(part-1) 0.12 0.36 0.36])
    plot(tvec(1:i),r(1,1:i),'k')
    hold on
    plot(tvec(1:i),r(2,1:i),'k:')
    axis([0 tmax 0 rmax])
    xlabel('Time (sec)')
    ylabel('Rate (Hz)')
    legend('r1','r2')
    subplot('Position',[0.62-0.5*(part-1) 0.62 0.36 0.36])
    plot(r(1,1:i),r(2,1:i),'k--')
    
    
end


annotation('textbox',[0.0 0.98 0.02 0.02],'String','B','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.5 0.98 0.02 0.02],'String','C','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.0 0.48 0.02 0.02],'String','D','LineStyle','none','FontSize',16,'FontWeight','bold')
annotation('textbox',[0.5 0.48 0.02 0.02],'String','E','LineStyle','none','FontSize',16,'FontWeight','bold')
