function [sta, tcorr] = STA(Iapp, spikes, dt, tminus, tplus)
% STA   Spike-triggered-average
%   [sta, tcorr] = STA(Iapp, spikes, dt, tminus, tplus)
%   Returns the spike-triggered average, sta, which is the sum of all the
%   input currents surrounding the time of each spike, divided by the 
%   number of spikes. 
%   The time-window used is described by the set of time-bins in tcorr.
%   
%   STA requires arrays and variables as follows:
% 	
%   Iapp(Ncells,length(t)) noisy applied current to each cell used to
%   produce the response
%   t = 0:dt:tmax is the time array
%   Ncells is the number of cells recorded
%
%   spikes(Ncells,length(t)) is a binart array with the spike times in the 
%   appropriate time bins for each cell indicated as a 1, otherwise entries
%   should be 0.
%
%   dt is the time-step
%
%   tminus is the time to record the stimulus and average it before each
%   spike (default is -75ms)
%
%   tplus is the time to record the stimulus and average it after each
%   spike (default is (+25ms)
%
%   This function is used in Tutorial 2.1 in Chapter 2 of the textbook
%   An Introductory Course in Computational Neuroscience 
%   by Paul Miller, Brandeis University (2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('tminus'))
    tminus = 0.075; % How long before spike time to begin
end
if (~exist('tplus'))
    tplus = 0.025;  % How long after spike time to continue
end

nminus = ceil(tminus/dt); % Number of time points before zero
nplus = ceil(tplus/dt);   % Number of time points after zero
nt = length(Iapp);     % length of original data set
sum_I = zeros(1,nminus+nplus+1);  % STA will accumulate here
tcorr = -nminus*dt:dt:nplus*dt;   % Vector of time points for STA
Iapp = Iapp - mean(Iapp); % Removes mean applied current
spikeposition = find(spikes);  % Time bins for each spike
totalspikes = length(spikeposition)    % Total number of spikes
for spike = 1:totalspikes
    ispike = spikeposition(spike);     % ispike is the bin containing a spike
    imin = max(1,ispike-nminus);       % Bin to start measuring stimulus
    imax = min(nt,ispike+nplus);       % Bin to finish measuring
    % The following lines put the stimulus, Iapp, into bins shifted
    % by the spike time (ispike)
    for i = imin:imax
        sum_I(i-ispike+nminus+1) = sum_I(i-ispike+nminus+1) ...
            + Iapp(i)/totalspikes;
    end
end

% Finally normalize by the number of contributing spikes to obtain the
% average (the mean)
sta = sum_I/totalspikes;
