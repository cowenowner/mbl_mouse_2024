%Vector_Stats_SpecCalcsG.m 

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% GUI that obtains the controls that specify
% the calculations to be made

% Variables to be obtained and defaults:
%   Flag1A=0: Test for uniformity if =1
%   Flag1B=0: Test for vonMises distribution if =1
%   Flag2P=0: Generate plots if =1
%   Flag2A=0: Inference on vector mean, kappa known if =1
%   Flag2B=0: Inference on vector mean, kappa unknown if =1
%   Flag2C=0: Confidence interval on kappa if =1
%   Flag4A=0: Separate out two components of mixture if =1
%   alfa=0.05: Signficance level; must be in range (0.001, 0.25)
%   MeanKnown=0: Known or hypothesized vector mean direction
%   KappaKnown: Known concentration parameter
%   RtrnCode=0: Error return

% Set up the figure

Hfig_GUI = figure;       % Handle to figure window

set(Hfig_GUI, 'NumberTitle', 'off', ...
    'Name','Specify calculations to do', ...
    'Menubar','none', ...
    'Units','inches',  ...
    'Position', [0.3, 2.0, 5, 5.5]);

% Set up the entries

HC_1 = uicontrol(Hfig_GUI, 'Style','text', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 38, 65, 2], ...
    'String','Standard parameters will be estimated');

HC_FR = uicontrol(Hfig_GUI, 'Style','frame', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[8, 10, 82, 25]);

HC_2P = uicontrol(Hfig_GUI, 'Style','checkbox', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 35.5, 65, 2], ...
    'Value', Flag2P, ...
    'String','Generate plots (rose diagrams, Q-Q)', ...
    'Callback',['if get(HC_2P,''Value'') == 1,' ...
                  'Flag2P = str2num(''1'');,' ...
                'else,' ...
                  'Flag2P = str2num(''0'');,'...
                'end']);

HC_1A = uicontrol(Hfig_GUI, 'Style','checkbox', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 32, 65, 2], ...
    'Value', Flag1A, ...
    'String','Perform tests for uniform distribution', ...
    'Callback',['if get(HC_1A,''Value'') == 1,' ...
                  'Flag1A = str2num(''1'');,' ...
                'else,' ...
                  'Flag1A = str2num(''0'');,'...
                'end']);

HC_1B = uicontrol(Hfig_GUI, 'Style','checkbox', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 29, 65, 2], ...
    'Value', Flag1B, ...
    'String','Perform test for von Mises distribution', ...
    'Callback',['if get(HC_1B,''Value'') == 1,' ...
                  'Flag1B = str2num(''1'');,' ...
                'else,' ...
                  'Flag1B = str2num(''0'');,'...
                'end']);
    
HC_2A = uicontrol(Hfig_GUI, 'Style','checkbox', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 26, 65, 2], ...
    'Value', Flag2A, ...
    'String','Inference on vector mean, concentration known', ...
    'Callback',['if get(HC_2A,''Value'') == 1,' ...
                  'Flag2A = str2num(''1'');,' ...
                'else,' ...
                  'Flag2A = str2num(''0'');,'...
                'end']);
    
HC_2B = uicontrol(Hfig_GUI, 'Style','checkbox', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 23, 65, 2], ...
    'Value', Flag2B, ...
    'String','Inference on vector mean, concentration unknown', ...
    'Callback',['if get(HC_2B,''Value'') == 1,' ...
                  'Flag2B = str2num(''1'');,' ...
                'else,' ...
                  'Flag2B = str2num(''0'');,'...
                'end']);
    
HC_2C = uicontrol(Hfig_GUI, 'Style','checkbox', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 20, 68, 2], ...
    'Value', Flag2C, ...
    'String','Confidence interval on concentration (Kappa)', ...
    'Callback',['if get(HC_2C,''Value'') == 1,' ...
                  'Flag2C = str2num(''1'');,' ...
                'else,' ...
                  'Flag2C = str2num(''0'');,'...
                'end']);
    
HC_ALF1 = uicontrol(Hfig_GUI, 'Style','text', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 17, 65, 2], ...
    'String','Significance level (alfa) for inference');
HC_ALF2 = uicontrol(Hfig_GUI, 'Style','edit', ...
    'Units','characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[72, 17, 16, 2], ...
    'String',num2str(alfa), ...
    'Callback','alfa = str2num(get(HC_ALF2,''String''));');

HC_MN1 = uicontrol(Hfig_GUI, 'Style','text', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 14, 65, 2], ...
    'String','Vector mean (Theta, deg.), known or hypothesized');
HC_MN2 = uicontrol(Hfig_GUI, 'Style','edit', ...
    'Units','characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[80, 14, 8, 2], ...
    'String',num2str(MeanKnown), ...
    'Callback','MeanKnown = str2num(get(HC_MN2,''String''));');

HC_KP1 = uicontrol(Hfig_GUI, 'Style','text', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 11, 65, 2], ...
    'String','Concentration param. (Kappa), known ');
HC_KP2 = uicontrol(Hfig_GUI, 'Style','edit', ...
    'Units','characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[80, 11, 8, 2], ...
    'String',num2str(KappaKnown), ...
    'Callback','KappaKnown = str2num(get(HC_KP2,''String''));'); 

HC_2CMP = uicontrol(Hfig_GUI, 'Style','checkbox', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 7, 65, 2], ...
    'Value', Flag4A, ...
    'String','Estimate two components of mixture', ...
    'Callback',['if get(HC_2CMP,''Value'') == 1,' ...
                  'Flag4A = str2num(''1'');,' ...
                'else,' ...
                  'Flag4A = str2num(''0'');,'...
                'end']);

HPushCls = uicontrol(Hfig_GUI, 'Style','push', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[10, 2, 10, 2], ...
    'String','OK', ...
    'Callback','RtrnCode = str2num(''0'');close');

HPushCncl = uicontrol(Hfig_GUI, 'Style','push', ...
    'Units', 'characters', 'Horizontalalignment','left', ...
    'FontSize',str2num('10'),  'FontUnits','points', ...
    'Position',[80, 2, 10, 2], ...
    'String','Cancel', ...
    'Callback','RtrnCode = str2num(''9'');close');

uiwait




