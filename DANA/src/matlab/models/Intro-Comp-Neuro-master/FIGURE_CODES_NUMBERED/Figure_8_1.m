% Figure_8_1.m
%
% Code to plot the boundaries for long-term potentiation (LTP) or long-term
% depression (LTD) as a function of the firing rate of the presynaptic cell
% and the postsynaptic cell.
%
% This code produces Figure 8.1 of 
% An Introductory Course in Computational Neuroscience 
% by Paul Miller (Brandeis University, 2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_T = 15;                           % threshold rate for LTP
r_max = 35;                         % maximum rate for plots
%% Set up the plotting parameters
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf

%% Now loop through 4 rules, each in its own panel of the figure
for rule = 1:4
    col_num = mod(rule-1,2);        % Panel is in column 0 or column 1
    row_num = (rule-1-col_num)/2;   % Panel is in row 0 or row 1
    
    xpos = 0.11+0.47*col_num;       % x position of panel depends on column
    ypos = 0.575-0.47*row_num;      % y position of panel depends on row
    
    subplot('Position',[xpos ypos 0.36 0.36])   % Set position and size of panel
    
    % Now produce a shaded rectangle for LTP in light gray
    fill([r_T r_T r_max r_max], [r_T r_max r_max r_T],[0.8, 0.8, 0.8], ...
        'EdgeColor','none')
    hold on
    % Label this as LTP within the shaded rectangle
    annotation('textbox',[xpos+0.22 ypos+0.26 0.1 0.02],'LineStyle','none', ...
        'FontSize',16,'FontWeight','Bold','String','LTP')
    
    % Rules 2 and 3 include LTD for low presynaptic and high postsynaptic
    % firing rates.
    if ( ( rule == 2 ) || (rule == 3 ) )
        % Produce a dark gray shaded rectangle for LTD
        fill([r_T r_max r_max r_T],[0 0 r_T r_T],[0.4, 0.4, 0.4], ...
            'EdgeColor','none')
        % Label this as LTD within the shaded rectangle
        annotation('textbox',[xpos+0.22 ypos+0.08 0.1 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','LTD')
    end
    
    % Rules 2 and 4 include LTD for high presynaptic and low postsynatpic 
    % firing rates.
    if ( ( rule == 2 ) || (rule == 4 ) )
        % Produce a dark gray shaded rectangle for LTD
        fill([0 0 r_T r_T], [r_T r_max r_max r_T],[0.4, 0.4, 0.4], ...
            'EdgeColor','none')
        % Label this as LTD within the shaded rectangle
        annotation('textbox',[xpos+0.03 ypos+0.26 0.1 0.02],'LineStyle','none', ...
            'FontSize',16,'FontWeight','Bold','String','LTD')
        
    end
    
    % Now plot horizontal and vertical dotted lines at the threshold rate
    plot([r_T r_T],[0 r_max],':k')
    plot([0 r_max],[r_T r_T],':k')
    
    axis([0 r_max 0 r_max])             % Set the axes sizes
    set(gca,'XTick',[0 r_T])            % Only zero and threshold rates matter
    set(gca,'XTickLabel',{'0', 'r_T'})  % Label the x-tick marks
    set(gca,'YTick',[0 r_T])            % Only zero and threshold rates matter
    set(gca,'YTickLabel',{'0', 'r_T'})  % Label the y-tick marks

    % Only label the y-axis on the left column
    if ( col_num == 0 )
        ylabel('Postsynaptic Rate')
    end
    % Only label the x-axis on the bottom row
    if ( row_num == 1 )
        xlabel('Presynaptic Rate')
    end
    
    % Title each figure panel with the Rule number 
    title(strcat(['Rule ',num2str(rule)]))
    set(gca,'Layer','top')              % Ensures axes overlay the "fills"
end

%% Finally label the figures A-D at the appropriate points
annotation('textbox',[0.00 0.98 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','A')
annotation('textbox',[0.50 0.98 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','B')
annotation('textbox',[0.00 0.5 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','C')
annotation('textbox',[0.50 0.5 0.02 0.02],'LineStyle','none', ...
    'FontSize',16,'FontWeight','Bold','String','D')
