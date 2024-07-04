function [O] = NPXL_Create_Spike_Seconds_npy_File(spike_times_fname, sFreq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ==0 
    spike_times_fname = 'spike_times.npy';
    sFreq = 30000;
end

O = [];
[p,n] = fileparts(spike_times_fname);

T = readNPY(spike_times_fname); % these are not (*$()*! times - they are records. Geez.
t_sec = double(T)/sFreq;
if any(diff(t_sec)<0)
    error('Timestamps are not in perfect ascending order. This should never happen.')
end
%     save(fullfile(DATA_DIR,'spike_seconds.mat'),'t_sec')
writeNPY(t_sec,fullfile(p,'spike_seconds.npy')); % these are not (*$()*! times - they are records. Geez.



