% Does temperature alter neural activity in the intermediate hippocampus?
% 
% Cowen Anahat Mario 2024
%%
clearvars % Clear out variables in the workspace. Start fresh.
% where is the data?
data_folder = 'D:\M564\'; % D:\M564\M564_2024_07_09_NP2_4_Shanks_Anahat_Mario_Test_01_g0
npxl_top_dir_name = 'M564_2024_07_09_NP2_4_Shanks_Anahat_Mario_Test_01_g0';
probe_depth_adjust_mm = 4; % set this to a non-zero value ONLY if there was a mistake when setting the probe depth in the AllSpikes file.
kilosort_dir_name = 'kilosort'; % a subfolder where the AllSpikes file from spike sorting lives.
figures_dir = './Figures';
SAVE_FIGURES = true;
wheel_hot_cold_neut = [0 5/3; 5/3 10/3; 10/3 5 ];
wheel_hot_cold_neut_labs = {'C' 'H' 'N'};
wheel_hot_cold_neut_colors = [.1 .1 .9; .9 .1 .1; .7 .7 .7];

if ~isfolder(figures_dir)
    mkdir(figures_dir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure out directories and filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D] = NPXL_get_file_names(data_folder,npxl_top_dir_name,kilosort_dir_name );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NIDQ = NPXL_Extract_NIDQ(D.nidq_bin_file_path,1:4);
running_wheel_data = [NIDQ.t_sec(:) NIDQ.data_V(1,:)'];
running_wheel_temp = zeros(size(NIDQ.t_sec(:)));
running_wheel_start_end_sec = [];
IX = false(length(running_wheel_temp),3);
for ii = 1:3
    IX(:,ii) = running_wheel_data(:,2) > wheel_hot_cold_neut(ii,1) & running_wheel_data(:,2) <= wheel_hot_cold_neut(ii,2);
    running_wheel_temp(IX(:,ii)) = ii;
    IX(1,ii) = 0;
    dI = [0; diff(IX(:,ii))];
    dI(1) = 0;    dI(end) = 0;
    stIX = find(dI == 1); edIX = find(dI == -1);
    end_ix = min([length(stIX) length(edIX)]);
    stIX = stIX(1:end_ix); edIX = edIX(1:end_ix); 
    % start and end times when the mouse entered the zone for PETHs and
    % things
    running_wheel_start_end_sec{ii} = [running_wheel_data(stIX,1) running_wheel_data(edIX,1)];
end

time_in_temp_spots_sec = sum(IX)/25000;

figure
bar(time_in_temp_spots_sec)
set(gca,'XTickLabel',wheel_hot_cold_neut_labs)
title('Time spent in each temp')
ylabel('sec')

%% 
CHOOSE_MOVING_TIMES = false;
if CHOOSE_MOVING_TIMES
    figure
    plot(running_wheel_data(:,1), running_wheel_data(:,2))
    axis tight
    hold on

    [x,y] = ginput_cowen;
    start_end_move_sec = [x(1:2:end) x(2:2:end)];
 

else
    start_end_move_sec = [         0.969609783541756           12.604927186043
           19.715398932016          70.1351076761881
           114.73715771911          121.847629465083
          141.886231658279          156.753581672586
           196.83078605898          209.112509983842
          287.974105711906          307.366301382742
          350.675538380941          365.542888395248
           399.80243408039          446.990110212757
          515.509201583042          531.669364642071
          584.674699475688          591.785171221661
          618.934245160831          636.387221264583
          671.939579994448          695.856621321811
          725.591321350426          735.287419185843
          777.950249661681          780.535875751126
          789.585567064183          801.867290989045
           804.45291707849          813.502608391546
          842.590901897799          844.530121464883
          876.204041060581          882.021699761831
           908.52436717864          944.723132430866
          951.833604176839          959.590482445173];

end
[running_wheel_data_restricted ,BADINT] = Restrict(running_wheel_data, start_end_move_sec);
GOODINT = ~BADINT;

IX2 = IX(GOODINT,:);
time_in_temp_spots_good_sec = sum(IX2)/25000;
figure
bar(time_in_temp_spots_good_sec)
set(gca,'XTickLabel',wheel_hot_cold_neut_labs)
title('Time spent in each temp RESTRICTED')
ylabel('sec')

%% 

%% 
% if ~isempty(D.lfp_bin_file_path)
%     LFP = NPXL_Extract_LFP(D.lfp_bin_file_path,1:20:385);
% end

load(fullfile(D.kilosort_dir,'AllSpikes.mat')); % Loads SP structure. This was generated after spike sorting.
% resort the SP data by depth
[~,dpth_six] = sort([SP.neuropixels_depth_uM]);
SP = SP(dpth_six);

TS_uS = {SP.t_uS}; % A cell array of ts is handy.

EVT = [];
for ii = 1:length(D.event_files)
    EVT.t_sec{ii} = load(D.event_files{ii});
    EVT.name{ii} = 'dontknow';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depths along probe for double checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(LFP.INFO.Ch,LFP.INFO.Depth_uM ,'o-')
% ylabel('depth uM');xlabel('ch')
% title('Depths per probe')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raster plot of spiking activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
ax(1) = subplot(2,1,1);
plot_raster(TS_uS,[],'time_divisor',1e6,'make_colorful',false)
axis ij
hold on
if ~isempty(EVT)
    for ii = 1:length(EVT)
        plot(EVT.t_sec{ii},(ii-1)*ones(size(EVT.t_sec{1})),'r>')
        hold on
    end
end
xlabel('sec')
ylabel('Neuron ID')
title('Raster')

% LFP.data_uV_car = LFP.data_uV - mean(LFP.data_uV);

ax(2) = subplot(2,1,2);
t_sec = running_wheel_data(1:10:end,1);
pos = running_wheel_data(1:10:end,2);
heats = running_wheel_temp(1:10:end);
plot(t_sec,pos)
axis tight
hold on
ylabel('Wheel position')
for ii =1:3
    IIX = heats == ii;
    plot(t_sec(IIX),pos(IIX),'.', 'Color',wheel_hot_cold_neut_colors(ii,:))
end
% 
% plot_LFP(LFP.data_uV_car(1:4,1:10:end)',LFP.t_sec(1:10:end))
% title('LFP CAR')
linkaxes(ax,'x')



% Create a neuron x time matrix. Let's call it 'Q'
big_binsize_ms = 500;
[Q, edges_uS] = histcounts_cowen(TS_uS, 'binsize', big_binsize_ms*1000);
Q_uS = edges_uS(1:end-1) + (edges_uS(2) - edges_uS(1))/2;
[~, BADIX] = Restrict(Q_uS, start_end_move_sec*1e6);
Q_good = Q(~BADIX,:);
Q_good_sec = Q_uS(~BADIX)/1e6;


size(Q)

figure;
subplot(2,1,1)
imagesc(Q_uS/60e6,[],Q')
ylabel('Neuron');xlabel('time (min)')
colorbar
title(sprintf('Original data binned at %d ms', big_binsize_ms))
subplot(2,1,2)
imagesc(Q_good_sec/60,[],Q_good')
ylabel('Neuron');xlabel('time (min)')
colorbar
title(sprintf('Restricted data binned at %d ms', big_binsize_ms))

% Pull out the temperature intervals from the Q_good matrix.
Qnew = [];
for ii = 1:3
    [~,IXtemp] = Restrict(Q_good_sec,sortrows(running_wheel_start_end_sec{ii}));
    Qnew{ii} = Q_good(IXtemp,:);
    QMn(ii,:) = mean(Qnew{ii});
end
MN_rates = mean(QMn,2);
figure
bar(MN_rates)

% csvwrite('C:\Temp\junk2.csv',QMn)



figure
for ii = 1:3
    subplot(3,1,ii)
    imagesc(Qnew{ii}')
    title(wheel_hot_cold_neut_labs{ii})
end
figure
bar(QMn')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Occupancy at each position on the running wheel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hist(running_wheel_data(:,2), 100)
xlabel('wheel postion')
ylabel('count')
title('Occupancy at each wheel position')
pubify_figure_axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Event level analysis of the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(EVT)
    % Event analysis: Spikes
    % Choose the event you want to analyze here.
    E_uS = EVT.t_sec{2}*1e6;

    % Plot some PETHS!

    for iC = 1:length(TS_uS)
        figure
        PETH_raster(TS_uS{iC}/100,E_uS/100,100,2000,4000);
        sgtitle(sprintf('Stim start, %1.2f uM',SP(iC).neuropixels_depth_uM ))
        if SAVE_FIGURES
            saveas(gca,fullfile(figures_dir,sprintf('PETH_%d.png',iC)));
        end
    end

    % On the EEG
    % function [M, x_axis, fh, OUT, EEG_sec_data] = PETH_EEG(EEG_sec_data, sFreq ,alignments_ts_sec, time_before_sec, time_after_sec, option, option_parameters)

    for iCh = 1:Rows(LFP.data_uV)
        figure
        PETH_EEG([LFP.t_sec(:) LFP.data_uV_car(iCh,:)'],LFP.sFreq,E_uS/1e6, .5,3.2,{'sta_plot'})
        sgtitle(sprintf('Depth %1.2f uM',LFP.INFO.Depth_uM(iCh)))
        if SAVE_FIGURES
            saveas(gca,fullfile(figures_dir,sprintf('LFP_%d.png',iC)));
        end
    end


    %
    for iCh = 1:Rows(LFP.data_uV)
        PETH_EEG([LFP.t_sec(:) LFP.data_uV_car(iCh,:)'],LFP.sFreq,E_uS/1e6, 2,2,{'power_spectrum'})
        sgtitle(sprintf('Depth %1.2f uM',LFP.INFO.Depth_uM(iCh)))
        if SAVE_FIGURES
            saveas(gca,fullfile(figures_dir,sprintf('LFP_%d.png',iC)));
        end
    end
end

%% figure;
imagesc(Q_uS/60e6,[],Q')
ylabel('Neuron');xlabel('time (min)')
colorbar
title(sprintf('Original data binned at %d ms', big_binsize_ms))

