function [H, x_log] = HistISI_ksdensity(TS_uS, varargin)
% A SMOOTHED ISI histogram plot.
% H = HistISI_ksdensity(TS_uS, parameters)
% INPUTS:
%      TS = a single vector of timestamps assumed in uS.
%
% OUTPUTS:
%      H = histogram of ISI
%      x_log = bin centers
%
% PARAMETERS:
%
% If no outputs are given, then plots the figure directly.
%
% Cowen 2023

x_log = nan;
ISI_ms = diff(TS_uS/1000) + eps;
if isempty(ISI_ms) || any(ISI_ms<=eps)
   warning('ISI contains no data or has <=0 values!');
   return
end   

[H, x_log] = ksdensity(log10(ISI_ms));

%-------------------
if nargout == 0 
    plot(x_log,H,'LineWidth',3)
    axis tight
    hold on
    plot([0 0],[get(gca,'YLim')],'r:','LineWidth',3)
    xlabel('log10(msec)'); ylabel('density(p)')
    pubify_figure_axis
end