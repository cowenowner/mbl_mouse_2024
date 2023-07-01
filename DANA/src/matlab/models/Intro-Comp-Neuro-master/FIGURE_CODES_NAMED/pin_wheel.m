% pin_wheel.m Cartoon of a small pin-wheel map in V1.
%
%  This code was used to produce Figure 6.18 in the book:
%  An Introductory Course in Computational Neuroscience
%  by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% For the figure in a black and white textbook set coloroff to 1.
% Otherwise set coloroff to 0 for a better figure.
coloroff = 1;

% Set up a rectangular array of neurons for two pinwheels.
Nx = 400;
Ny = 200;
% Define coordinates of the first pinwheel center.
% The second one will simply be a reflection of the first producing a
% mirror copy shifted in the x-direction.
center1x = Nx/4+0.5;
center1y = Ny/2+0.5;

% angle_pref will contained the preferred angle for an oriented bar
% stimulus that can range from -pi/2 to + pi/2.
angle_pref = zeros(Nx,Ny);
for i = 1:Nx/2;
    % For positive x then the arctangent of y/x produces an angle between
    % -pi/2 and pi/2. In this case we divide by two to obtain an angle
    % between -pi/4 and +pi/4.
    if ( i > center1x )
        for j = 1:Ny;
            angle_pref(i,j) = 0.5*atan((j-center1y)/(i-center1x));
        end
    else;   % For negative x
        for j = 1:center1y;     % For negative y, a negative angle
            angle_pref(i,j) = -pi/2+0.5*atan((j-center1y)/(i-center1x));
        end
        for j = ceil(center1y):Ny; % For positive y, a positive angle
            angle_pref(i,j) = pi/2+0.5*atan((j-center1y)/(i-center1x));
        end
    end
end

% Now produce a reflected pinwheel in the right half of x-coordinates
angle_pref(Nx:-1:Nx/2+1,:) = angle_pref(1:Nx/2,:);

% In the plot below the absolute value of the preferred angle is used to
% prevent a discontinuity between -pi/2 and +pi/2 when they represent the
% same orientation of a bar.
% If a color plot is used then a full rotation through distince colors in
% the RGB plane can be used to produce a more accurate figure where for
% example -pi/4 is distinct from +pi/4.
% When using image or imagesc, the transpose of the matrix is needed to 
% convert coordinates (x,y) into (column no., row no.) instead of 
% (row, column), which is the standard matrix notation.
figure(1)
clf
subplot('Position',[0.025 0.05 0.95 0.9])
set(gca,'YDir','normal')
if ( coloroff )
    imagesc(abs(angle_pref'))
    colormap(gray);         % Sets the colors the range of values map onto
else    
    imagesc(angle_pref');   
    npoints = 10000;     % Determines the gradation of the color scale
    vec = 2*pi/npoints*[0:npoints]';    % A vector from 0 to 2.pi
    reds = (0.5*(1+cos(vec)));            % Redness cycles along the vector
    greens = (0.5*(1+cos(vec+2*pi/3)));   % Greenness offset from redness
    blues = (0.5*(1+cos(vec-2*pi/3)));    % Blueness offset from redness and greenness.
    map = [reds greens blues];          % New colormap in RGB values
    colormap(map);          % Sets the colors the range of values map onto
end

% By default Matlab plots an image of a matrix with the rows ordered with
% row 1 at the top. If the rows represent y-coordinates, then we need to
% set the y-direction as 'normal' not 'reverse'.
set(gca,'YDir','normal')

% Finally add some labels of the angles at specific points on the plot
annotation('textbox',[0.08 0.54 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','Color','k','String','\pi/2')
annotation('textbox',[0.87 0.54 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','Color','k','String','\pi/2')
annotation('textbox',[0.5 0.54 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','Color','w','String','0')
annotation('textbox',[0.23 0.75 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','\pi/4')
annotation('textbox',[0.23 0.25 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','-\pi/4')
annotation('textbox',[0.7 0.75 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','\pi/4')
annotation('textbox',[0.7 0.25 0.02 0.02],'LineStyle','none','FontSize',16,'FontWeight','Bold','String','-\pi/4')
set(gca,'XTick',[]);
set(gca,'YTick',[]);
hold on
plot(center1x,center1y,'xk')
plot(3*Nx/4+0.5,Ny/2+0.5,'xk')

