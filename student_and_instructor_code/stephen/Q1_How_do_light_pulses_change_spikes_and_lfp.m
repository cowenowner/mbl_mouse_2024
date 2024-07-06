% Do LED blinky lights alter brain activity in the hippocampus?
%% 
clearvars
data_folder = 'D:\Data\';
npxl_top_dir_name = 'PhotoPixelsStrobe_g0';
kilosort_dir_name = 'kilosort_cowen';

% Figure out directories and filenames
[D] = WH_load_NPXL_data(data_folder,npxl_top_dir_name,kilosort_dir_name );

% Load data
LFP = NPXL_Extract_LFP(D.lfp_bin_file_path,1:20:385);

load(fullfile(D.kilosort_dir,'AllSpikes.mat')); % Loads SP structure. This was generated after spike sorting.
TS_uS = {SP.t_uS}; % A cell array of ts is handy.

EVT = [];
for ii = 1:length(D.event_files)
    EVT.t_sec{ii} = load(D.event_files{ii});
    EVT.name{ii} = 'dontknow';
end
% 
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
