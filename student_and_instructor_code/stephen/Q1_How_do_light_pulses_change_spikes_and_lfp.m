% Do LED blinky lights alter brain activity in the hippocampus?
%% 
clearvars
data_folder = 'D:\Data\';
npxl_top_dir_name = 'PhotoPixelsStrobe_g0';
kilosort_dir_name = 'kilosort_cowen';
figures_dir = './Figures';
SAVE_FIGURES = true;
% Figure out directories and filenames
[D] = NPXL_get_file_names(data_folder,npxl_top_dir_name,kilosort_dir_name );

% Load data
NIDQ = NPXL_Extract_NIDQ(D.nidq_bin_file_path,1:4);

LFP = NPXL_Extract_LFP(D.lfp_bin_file_path,1:20:385);

load(fullfile(D.kilosort_dir,'AllSpikes.mat')); % Loads SP structure. This was generated after spike sorting.
TS_uS = {SP.t_uS}; % A cell array of ts is handy.

EVT = [];
for ii = 1:length(D.event_files)
    EVT.t_sec{ii} = load(D.event_files{ii});
    EVT.name{ii} = 'dontknow';
end
% Depths
figure
plot(LFP.INFO.Ch,LFP.INFO.Depth_uM ,'o-')
ylabel('depth uM');xlabel('ch')
title('Depths per probe')

%  For exploration.

figure
ax(1) = subplot(3,1,1);
plot_raster(TS_uS,[],'time_divisor',1e6)
hold on
plot(EVT.t_sec{1},ones(size(EVT.t_sec{1})),'bo')
plot(EVT.t_sec{2},2*ones(size(EVT.t_sec{2})),'ro')
plot(EVT.t_sec{3},3*ones(size(EVT.t_sec{3})),'go')
xlabel('min')

ax(2) = subplot(3,1,2);
plot_LFP(LFP.data_uV(1:4,1:10:end)',LFP.t_sec(1:10:end))

ax(3) = subplot(3,1,3);
imagesc(LFP.t_sec(1:10:end),[],LFP.data_uV(:,1:10:end))
axis tight
xlabel('sec')
linkaxes(ax,'x')

% Create a neuron x time matrix. Let's call it 'Q'
big_binsize_ms = 500;
[Q, edges_uS] = histcounts_cowen(TS_uS, 'binsize', big_binsize_ms*1000);
Q_uS = edges_uS(1:end-1) + (edges_uS(2) - edges_uS(1))/2;

size(Q)

figure;
imagesc(Q_uS/60e6,[],Q')
ylabel('Neuron');xlabel('time (min)')
colorbar
title(sprintf('Original data binned at %d ms', big_binsize_ms))


% Event analysis: Spikes
E_uS = EVT.t_sec{3}*1e6;
d_ms = diff(E_uS)/1e3;
cnt = 1; end_gix = [];
for ii = 1:length(d_ms)-1
    if d_ms(ii) > d_ms(ii+1) && d_ms(ii) > 100
       end_gix(cnt) = ii;
       cnt = cnt + 1;
    end
end
start_gix = end_gix +1;
% Double check the timestamps to make sure the events are correct.
figure
plot(E_uS,zeros(size(E_uS)),'bo')
hold on
plot(E_uS(start_gix),zeros(size(E_uS(start_gix))),'g+')
plot(E_uS(end_gix),zeros(size(E_uS(end_gix))),'r*')

% Plot some PETHS!

for iC = 1:length(TS_uS)
    figure
    PETH_raster(TS_uS{iC}/100,E_uS(start_gix)/100,100,2000,4000);
    sgtitle(sprintf('Stim start, %1.2f uM',SP(iC).depth_of_probe_tip_uM ))
    if SAVE_FIGURES
        saveas(gca,fullfile(figures_dir,sprintf('PETH_%d.png',iC)));
    end
end

% On the EEG
% function [M, x_axis, fh, OUT, EEG_sec_data] = PETH_EEG(EEG_sec_data, sFreq ,alignments_ts_sec, time_before_sec, time_after_sec, option, option_parameters)

for iCh = 1:Rows(LFP.data_uV)
    figure
    PETH_EEG([LFP.t_sec(:) LFP.data_uV(iCh,:)'],LFP.sFreq,E_uS(start_gix)/1e6, .5,3.2,{'sta_plot'})
    sgtitle(sprintf('Depth %1.2f uM',LFP.INFO.Depth_uM(iCh)))
    if SAVE_FIGURES
        saveas(gca,fullfile(figures_dir,sprintf('LFP_%d.png',iC)));
    end
end
