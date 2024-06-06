% function [status,cmdout] = NPXL_test_alignment_post_tprime()
%% A script to test alignments.
NI_event_orig_file = '100422_DANA_6500uM_bank0_g0_tcat.nidq.xa_1_0.txt';
NI_event_sync_file = 'synced_100422_DANA_6500uM_bank0_g0_tcat.nidq.xa_1_0.txt';
original_spike_rec_file = '..\spike_times.npy';
original_spike_time_sec_file = '..\spike_seconds.npy';
synced_spike_time_sec_file = '..\synced_spike_seconds.npy';


Orig = load(NI_event_orig_file);
Synced = load(NI_event_sync_file);

figure;plot(Orig-Synced,'o')
ylabel('sec drift')
plot_horiz_line_at_zero
ylabel('sec (orig-sync)')
title('orig - sync square sync pulse diff from imec and NI')
xlabel('sec')

T_sec_orig = readNPY(original_spike_time_sec_file);
T_sec_synced = readNPY(synced_spike_time_sec_file);

figure;plot(T_sec_orig-T_sec_synced,'o')
ylabel('sec drift')
plot_horiz_line_at_zero
ylabel('sec (orig-sync)')
title('orig - sync spike times diff from imec and NI')
xlabel('sec')

% look the output from Phy: original_spike_rec_file

T_rec = readNPY(original_spike_rec_file);
fprintf('Min of recording %d\n',(T_rec(end)/30000)/60);

% For DANA only - see what alignes with the sync pulses the best....
% Unique to DANA
synced_scan_time_sec_file = 'synced_100422_DANA_6500uM_bank0_g0_tcat.nidq.xa_3_0.txt';
orig_scan_time_sec_file = '100422_DANA_6500uM_bank0_g0_tcat.nidq.xa_3_0.txt';

fscv_scan_orig_sec = load(orig_scan_time_sec_file);
fscv_scan_sync_sec = load(synced_scan_time_sec_file);

% 
% figure
% plot(T_sec_orig,ones(size(T_sec_orig)),'k.')
% hold on


figure
PETH_raster(T_sec_orig*10000,fscv_scan_orig_sec*10000,4,200,200);
title('orig spike orig scan')

figure
PETH_raster(T_sec_orig*10000,fscv_scan_sync_sec*10000,4,200,200);
title('orig spike sync scan')

figure
PETH_raster(T_sec_synced*10000,fscv_scan_sync_sec*10000,4,200,200);
title('sync spike sync scan')


