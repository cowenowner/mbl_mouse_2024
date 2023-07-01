% Figure_7_4.m
%
% The fixed points for two coupled equations are shown as the intersection
% of two nullclines, each nullcline representing the fixed points of one
% equation if the other variable were held constant. The equations are:
%
% tau.dr1/dr = -r1 + Ws.r1 + Wx.r2 - theta
% tau.dr2/dt = -r2 + Ws.r2 + Wx.r1 - theta
%
% The time constant does not impact the nullclines, which are obtained by
% setting the left-hand side of each equation to zero.
%
% Two additional conditions cause the nullclines to have a "Z" shape,
% rather than each be a perfect straight line. These conditions are;
% 1) r1 and r2 cannot be less than 0.
% 2) r1 and r2 cannot be more than rmax.
%
% This code is used to produce Figure 7.4 in the book:
% An Introductory Course in Computational Neuroscience 
% by Paul Miller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for part = 1:2;         % part 1 is the integrator, part 2 used for Fig. 7.4
    % The connection strenghts and firing-rate threshold can vary
    if ( part == 1 )
        Ws = 0.975;     % recurrent connection strength
        Wx = -0.025;    % cross-connection strength
        theta = -0.5;   % threshold (negative means spontaneously active)
        rmax = 60;      % maximum rate
    else
        Ws = 1.05;      % recurrent connection strength
        Wx = -0.1;      % cross-connection strength
        theta = 2;      % threshold
        rmax = 60;      % maximum rate
    end
    
    figure(part)
    clf
    
    % r2-Nullcline, plot of r2 as a function of r1 with dr2/dt = 0
    r1vals = -100:0.1:100;
    r2null = (Wx*r1vals - theta)/(1-Ws);
    
    % Exclude values where r<0 or r>rmax
    allowed_vals = find((r2null >= 0 ).*( r2null <= rmax ));
    
    % Plot the nullcline where dr2/dt = 0
    plot(r1vals(allowed_vals),r2null(allowed_vals),'k:');
    hold on
    
    % r1-Nullcline, plot of r1 as a function of r2 with dr1/dt = 0
    r2vals = -100:0.1:100;
    r1null = (Wx*r2vals - theta)/(1-Ws);

    % Exclude values where r<0 or r>rmax
    allowed_vals = find((r1null >= 0 ).*( r1null <= rmax ));
    
    % Plot the nullcline where dr1/dt = 0
    plot(r1null(allowed_vals),r2vals(allowed_vals),'k');
    hold on
    
    legend('dr_{2}/dt = 0', 'dr_{1}/dt = 0')
    
    % Now remainder of r2 nullcline:
    % vals of r1 where r2 = 0 produces negative dr2/dt
    allowed_vals = find(Wx*r1vals - theta < 0 );
    plot(r1vals(allowed_vals),zeros(1,length(allowed_vals)),'k:');
    hold on
    
    % Now vals of r1 where r2 = rmax produces positive dr2/dt
    allowed_vals = find( (Ws-1)*rmax + Wx*r1vals - theta > 0 );
    plot(r1vals(allowed_vals),rmax*ones(1,length(allowed_vals)),'k:');
    
    hold on
    
    
    % Now remainder of r1 nullcline:
    % vals of r2 where r1 = 0 produces negative dr1/dt
    allowed_vals = find(Wx*r2vals - theta < 0 );
    plot(zeros(1,length(allowed_vals)),r2vals(allowed_vals),'k');
    hold on
    
    % Now vals of r2 where r1 = rmax produces positive dr1/dt
    allowed_vals = find( (Ws-1)*rmax + Wx*r2vals - theta > 0 );
    plot(rmax*ones(1,length(allowed_vals)),r2vals(allowed_vals),'k');
    
    % Finish up the plot with axis-labels
    hold on    
    axis([-20 70 -20 70])
    xlabel('r_{1} (Hz)')
    ylabel('r_{2} (Hz)')
    
    % Finally indicate the fixed points for Figure 7.4 with a solid circle
    % for stable ones and an open circle for unstable ones
    if ( part == 2 )
        plot(60,0,'.k','MarkerSize',40)
        plot(0,0,'.k','MarkerSize',40)
        plot(0,60,'.k','MarkerSize',40)
        plot(0,40,'ok','MarkerSize',10)
        plot(40,0,'ok','MarkerSize',10)
    end
end
