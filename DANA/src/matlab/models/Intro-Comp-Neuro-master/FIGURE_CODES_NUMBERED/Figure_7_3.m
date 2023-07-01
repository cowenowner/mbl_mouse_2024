% Figure_7_3.m
%
% Plots the fixed points of r as a function of recurrent
% feedback, W.
%
% For each value of W, the code steps through values of r and checks the
% value of dr/dt. If dr/dt passes through zero then the value of r at the
% fixed point (where dr/dt = 0) is interpolated.
%
% This code is used to produce Figure 7.3 in the book:
% An Introductory Course in Computational Neuroscience,
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 0.01;                  % Step-size for finding fixed points
r_in = 0:dr:100;            % Set of rates for finding fixed points
Npts = length(r_in);        % No. of points used to calculate fixed points

rmax = 100;                 % Maximum rate in f-I curve
I_sigma = 20;               % Range of inputs to change rate (inverse steepness)
Ith = 50;                   % Input for half-maximum rate
% Define the firing rate curve as a sigmoid function
f_of_r = @(x) rmax./(1+(exp(-(x-Ith)/I_sigma)));

Wvals = [0.5:0.0005:1.5];   % Set of values of feedback, W, to use one at a time
Nvals = length(Wvals);      % Number of subplots, one per value of W

fixed_points = zeros(3,Nvals);  % Up to 3 fixed points per W value allowed

%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf

%% Loop through different values of W, calculated dr/dt as a function of r
%   then plot the resulting curve for each value of W
for i = 1:Nvals;                % Loop through W
    Wrec = Wvals(i);            % Value of W for this plot
    
    input = Wrec*r_in;          % Input is recurrent feedback only
    
    % Calculate dr/dt as a vector from the dynamical equations as:
    % dr/dt = (-r + f(I))/tau    
    drdt = (-r_in + f_of_r(input))/tau;
       
    % Now loop through values of dr/dt (a vector with a value at each r) to
    % see where dr/dt passes through zero.
    for j = 2:Npts
        if ( ( drdt(j) <= 0 ) && ( drdt(j-1) > 0 ) )    % stable fixed point
            if ( j < Npts/2 )
                fixed_points(1,i) = dr*(j + drdt(j)/(drdt(j-1)-drdt(j)));
            else
                fixed_points(3,i) = dr*(j + drdt(j)/(drdt(j-1)-drdt(j)));                
            end
        end
        if ( ( drdt(j) >= 0 ) && ( drdt(j-1) < 0 ) )    % unstable fixed point
                fixed_points(2,i) = dr*(j + drdt(j)/(drdt(j-1)-drdt(j)));
        end
    end
end

% low_vals is set of stable fixed points of low rate
low_vals = find(fixed_points(1,:) > 0 );
plot(Wvals(low_vals),fixed_points(1,low_vals),'k')
hold on
% unstable_vals is set of unstable fixed points
unstable_vals = find(fixed_points(2,:) > 0 );
plot(Wvals(unstable_vals),fixed_points(2,unstable_vals),'k:')
% high_vals is set of stable fixed points of high rate
high_vals = find(fixed_points(3,:) > 0 );
plot(Wvals(high_vals),fixed_points(3,high_vals),'k')
legend('stable','unstable')
xlabel('Feedback strength, W')
ylabel('Firing rate, r (Hz)')
