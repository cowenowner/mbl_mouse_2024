function PHY_fix_synced_spike_times_in_AllSpikes(allspikes_file,spike_seconds_file, spike_clusters_file)
% Fixes my mistake in the Allspikes file when I used the synced_spike_seconds.npy
%
% spike_clusters.npy  spike_seconds.npy
%
if nargin < 1
    allspikes_file = 'AllSpikes.mat';
end
if nargin < 2
    spike_seconds_file = 'spike_seconds.npy';
end
if nargin < 3
    spike_clusters_file = 'spike_clusters.npy';
end

close all

out_allspikes = strrep(allspikes_file,'.mat','_fixed.mat');
load(allspikes_file)
SPold = SP;
C = readNPY(spike_clusters_file);
T_usec = readNPY(spike_seconds_file);
T_usec = T_usec*1e6;


for iC = 1:length(SP)
    sp_ix = C == SP(iC).cluster_id;
    SP(iC).t_uS = T_usec(sp_ix);
    SP(iC).note = 'updated using PHY_fix_synced_spike_times_in_AllSpikes';
end
save(out_allspikes,'SP')

if 1
    % Specific to DANA data.
    %%
    load('Events.mat')
    if contains(pwd,'Rat425')
        stim_times_sec = EVT.stim.time_sec;
        EVT.stimstartstop.time_sec = EVT.fscv.time_sec;
        EVT.fscv.time_sec = EVT.sync3.time_sec;
    end
    all_spikes_t_uS = [];
    all_spikes_old_t_uS = [];
    % Test
    % only for DANA for now.
    for iC = 1:length(SP)
        all_spikes_t_uS = [all_spikes_t_uS; SP(iC).t_uS];
        all_spikes_old_t_uS = [all_spikes_old_t_uS; SPold(iC).t_uS];
    end
    all_spikes_t_uS = unique(all_spikes_t_uS);
    all_spikes_old_t_uS = unique(all_spikes_old_t_uS);
    mean(all_spikes_old_t_uS - all_spikes_t_uS)

    figure;plot(all_spikes_old_t_uS - all_spikes_t_uS)

    figure
    PETH_raster(all_spikes_t_uS/100,EVT.stim.time_sec*10000,4,200,200);
    title('STIM orig spike sync scan')
    saveas(gcf,'STIM orig to sync.png')

    figure
    PETH_raster(all_spikes_old_t_uS/100,EVT.stim.time_sec*10000,4,200,200);
    title('STIM sync spike sync scan')
    saveas(gcf,'STIM sync to sync.png')

    figure
    PETH_raster(all_spikes_t_uS(1:1:end)/100,EVT.fscv.time_sec(1:6:end)*10000,4,200,200);
    title('SCAN orig spike sync scan')
    saveas(gcf,'SCAN orig to sync.png')

    figure
    PETH_raster(all_spikes_old_t_uS(1:1:end)/100,EVT.fscv.time_sec(1:6:end)*10000,4,200,200);
    title('SCAN sync spike sync scan')
    saveas(gcf,'SCAN sync to sync.png')

    if 1


        T_synced_usec = readNPY('synced_spike_seconds.npy');
        T_synced_usec = T_synced_usec*1e6;


        figure
        PETH_raster(T_usec(1:4:end)/100,EVT.fscv.time_sec(1:4:end)*10000,4,200,200);
        title('SCAN ALLSPIKES orig spike sync scan')
        saveas(gcf,'SCAN ALLSPIKES orig to sync.png')

        figure
        PETH_raster(T_synced_usec(1:4:end)/100,EVT.fscv.time_sec(1:4:end)*10000,4,200,200);
        title('SCAN sync spike sync scan')
        saveas(gcf,'SCAN ALLSPIKES sync to sync.png')

    end
end


