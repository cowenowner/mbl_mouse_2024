function[]= NP2_LF_DataLoad_v2 (rat,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function NP2_LF_DATA_LOAD(PRM.RAT_NO,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function loads Neuropixel Data for vHC_STIM Project for each rat
% -Loads Stim params from stim file name
% -Loads Chan map -name is kilosortChanMap
% -Loads Meta Data
% -Loads LF data post stim for all channels across each condition - assumes
%   iHC, vHC, Burst Stim Conditions (Time is 150 ms can change parameter)
% -Removes the Noise by subtracting data from the time before stim
% -Sorts LF data on the basis of shank and gets the depth for each site
% -Saves data in .mat file for further processing
%
% NECESSARY INPUT IS RAT NO. - PRM.RAT_NO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS TO CHANGE
%SOFTWARE DIRS
PRM.TPRIME_DIR='C:\DATA\Neuropixel\Software\TPrime-win';

%DATA DIRS
PRM.DATA_DIR='C:\DATA\Neuropixel\NP2';
PRM.STORAGE_DIR='C:\DATA\Neuropixel\NP2\Processed'; %Where data is stored and worked form later
PRM.SUB_FOLDER = {'mPFC_L23_L5_NP2_g0'};

% Data Parameters
PRM.n_stim=25; %no. of stims per burst
PRM.epochs={'iHC','vHC','BurstStim'}; %Different epochs of experiment
PRM.chanFile='mPFC_L23_L5_NP2_g0_tcat.imec0.ap_kilosortChanMap.mat';
PRM.time_extract_postStim=0.15; %Time in sec to extract post stim 120ms
PRM.time_extract_preStim=0.05; %Time in sec to extract post stim 30ms
%SHANK AND AVERAGING DATA SPECS
PRM.shank_no=4; %No. of shanks
PRM.avg_no=5; %Averages across 5 sessions
PRM.max_depth=5125; %Depth of probe In micrometers (Calculated from sterotax coordinates -tip)
%IL PL Separation
PRM.PL_depth=3600; %Depth of PL-IL transition in micrometers

Extract_varargin; % overrides the defaults specified above.
%% STIM PARAMS LOAD
%Fig dirs
PRM.RAT_NO=rat;
PRM.RAT_DATA_DIR=fullfile(PRM.DATA_DIR,PRM.RAT_NO);
PRM.RAT_STORAGE_DIR=fullfile(PRM.STORAGE_DIR,PRM.RAT_NO);
PRM.IMG_DIR=fullfile(PRM.RAT_STORAGE_DIR,'Figures');

if ~exist(PRM.IMG_DIR,"dir")
mkdir(PRM.IMG_DIR)
end

%Load params from Original Stim File
params_path=dir(fullfile(PRM.RAT_DATA_DIR,'params*.mat'));
STIM_PRM=load(fullfile(params_path.folder,params_path.name));
[STIM.sorted_current,STIM.sorted_idx] = sort(STIM_PRM.PRM.output_current_array);
STIM.sorted_current_uA=(STIM.sorted_current-0.0986)./0.0113;                            %%% CHANGE THIS VALUE IF USING A DIFF STIMULATION



%% GET SESS INFO, RUN TPRIME
obj1 = SGLX_Class;
%GET SESSION FILES AND RUN TPRIME
SESS.ROOT=fullfile(PRM.RAT_DATA_DIR,PRM.SUB_FOLDER);                   %Root session folder
SESS.FILE_DIR=fullfile(SESS.ROOT,strrep(PRM.SUB_FOLDER, 'g0','g0_imec0'));          %Sub folder with .bin files
SESS.LF_meta_fname=strrep(PRM.SUB_FOLDER,'g0','g0_tcat.imec0.lf.meta');%LF meta name
SESS.LF_meta=obj1.ReadMeta(SESS.LF_meta_fname{1},SESS.FILE_DIR{1});             %meta data


%Run TPrime
NPXL_Sync_With_TPrime_SS('PRM_TPRIME_DIR',PRM.TPRIME_DIR,'PRM_ROOT_DATA_DIR', SESS.ROOT{1},... 
    'PRM_evt_files_extension',{'*nidq.x*5_0.txt'},'PRM_EVT_generate',false);

%Get Stim Times
STIM.fname=fullfile(SESS.ROOT,sprintf('synced_%s_tcat.nidq.xa_5_0.txt',PRM.SUB_FOLDER{1}));
stimes=readmatrix(STIM.fname{1});
%SPLIT INTO EPOCHS
STIM.epochs=PRM.epochs; %Change if order was switched
STIM.times{1}=stimes(1:PRM.n_stim); %iHC
STIM.times{2}=stimes(PRM.n_stim+1:2*PRM.n_stim); %vHC
STIM.times{3}=stimes(2*PRM.n_stim+5:5:end); %Burst Stim vHC

%Load Channel Map Info
ChanMap=load(fullfile(SESS.FILE_DIR{1},PRM.chanFile));
%% Load LFP DATA
%Initiate arrays
ALL_DATA={};
%LFP sess data
LFP.meta=SESS.LF_meta;
LFP.meta_fname = SESS.LF_meta_fname{1};
LFP.FILE_DIR= SESS.FILE_DIR{1};

%Set Sample Times
time_diff=PRM.time_extract_preStim+PRM.time_extract_postStim; %150 ms window
LFP.nSamp=floor(1.0 * time_diff*str2double(LFP.meta.imSampRate));
LFP.nChan=str2double(LFP.meta.nSavedChans);
ALL_DATA.badChans=[237,299,49];%LFP.meta.imStdby;
%badchans = [237,299,49];%str2num(ALL_DATA(1).badChans);  % Convert the bad channels list to numeric array


for sess=1:length(STIM.times)

    ALL_DATA(sess).Time=nan(length(STIM.sorted_current_uA),LFP.nSamp);
    ALL_DATA(sess).FiltVolts=nan(length(STIM.sorted_current_uA),LFP.nChan,LFP.nSamp); 
    ALL_DATA(sess).FiltData=nan(length(STIM.sorted_current_uA),LFP.nChan,LFP.nSamp);
    
    stim_times=STIM.times{sess};
    %Each Stim pulse
    for stim =1:length(stim_times)
        
        %Calculate sample times
        time_start=stim_times(stim)-PRM.time_extract_preStim; %in s
        LFP.Sample0=floor(time_start*str2double(LFP.meta.imSampRate));
        
        %Load LFP in sample range
        gain=str2double(LFP.meta.imChan0apGain);
        imax=str2double(LFP.meta.imMaxInt);
        vmax=str2double(LFP.meta.imAiRangeMax);

        LFP.binName=strrep(LFP.meta_fname,'.lf.meta','.lf.bin');
        LFP.data=obj1.ReadBin(LFP.Sample0,LFP.nSamp,LFP.meta,LFP.binName,LFP.FILE_DIR);
        LFP.data_Volts=obj1.ReadBinVolts(LFP.Sample0,LFP.nSamp,LFP.meta,LFP.binName,LFP.FILE_DIR);
        %LFP.data_Volts=*vmax/(imax*gain);%in V
        %LFP.data_Volts=obj1.GainCorrectIM(LFP.data,[1:384],LFP.meta);
        
        %Subtract Data from prev 150 ms
        s1=LFP.Sample0-LFP.nSamp;
        base_data=obj1.ReadBin(s1,LFP.nSamp,LFP.meta,LFP.binName,LFP.FILE_DIR);
        %base_data_Volts=obj1.ReadBinVolts(s1,LFP.nSamp,LFP.meta,LFP.binName,LFP.FILE_DIR);
        LFP.data_filt=LFP.data-base_data;
        LFP.data_filt_Volts=LFP.data_filt*vmax/(imax*gain);%in V
        %LFP.data_filt_Volts=LFP.data_Volts-base_data_Volts;
        
%         %Generates fig of data to look at base data
%           figure
%           tiledlayout(1, 3, "TileSpacing", "compact")
%           nexttile
%           imagesc(LFP.data)
%           nexttile
%           imagesc(base_data)
%           nexttile
%           imagesc(LFP.data_filt)
          
        LFP.time=time_start:(time_diff)/(LFP.nSamp-1):time_start+time_diff;%in ms
        LFP.sample_time=(LFP.Sample0+LFP.nSamp)/str2double(LFP.meta.imSampRate);
           
        %StoreLFP Data
        %ALL_DATA(sess).Volts(stim,:,:)= LFP.data_Volts;
        ALL_DATA(sess).Time(stim,:)=LFP.time;
        ALL_DATA(sess).FiltVolts(stim,:,:) = LFP.data_filt_Volts;
        ALL_DATA(sess).FiltData(stim,:,:) = LFP.data_filt;
        
    end
end

%% SORTING DATA BY SHANK
% Extract bad channels from the first session as an example

% Initialize LF_Shank structure
LF_Shank = struct();
badchans=ALL_DATA.badChans;

%% Processing for Each Shank
for ii = 1:PRM.shank_no
    % Initialize or reset the V and Avg fields for the current shank
    LF_Shank(ii).V = [];
    LF_Shank(ii).Avg = [];
    LF_Shank(ii).SortedDepths=[];

    % Extract shank channels and depths, and remove bad channels
    shank_channels = round(ChanMap.chanMap0ind(ChanMap.kcoords == ii));
    depths = ChanMap.ycoords(ChanMap.kcoords == ii);

    % Identify and remove bad channels for current shank
    valid_idx = ~ismember(shank_channels, badchans);
    valid_shank_channels = shank_channels(valid_idx);
    valid_depths = depths(valid_idx);

    % Sort depths 
    [sorted_depths, depth_sort_idx] = sort(valid_depths); 
     sorted_shank_channels = valid_shank_channels(depth_sort_idx) + 1; %Apply the sorted depths to shank channels (maintains order)
    %Calculates the depth of each tetrode (stereotax depth-site depth)
   LF_Shank(ii).SortedDepths =PRM.max_depth- sorted_depths;
    %sorted_shank_channels = min(sorted_shank_channels, size(LF_Volt, 2)); % Ensure indices are within valid range
    

    %% Processing for Each Session
    for sess = 1:2
        reshaped_V=[];
        LF_data = ALL_DATA(sess).FiltData;
        sel_LF = LF_data(STIM.sorted_idx, sorted_shank_channels, :);  % Select LF data based on sorted shank channels

        % Reshape for averaging
        reshaped_V = reshape(sel_LF, [PRM.avg_no, PRM.avg_no, length(sorted_shank_channels), size(sel_LF, 3)]);
        avg_lf = squeeze(mean(reshaped_V, 1)); % Calculate average of every 5 stims

        % Store results and smooth data with mov mean-over time domain
       LF_Shank(ii).V{sess} = movmean(sel_LF,12,3);  
       LF_Shank(ii).Avg{sess} = movmean(avg_lf,12,4);
       
    end
end
%Calculate the IL PL Depth difference chaneel
for shank=1:PRM.shank_no
    differences = abs(LF_Shank(shank).SortedDepths - PRM.PL_depth);
    [~, LF_Shank(shank).IL_end_chan] = min(differences);%index of the minimum difference
end    

%% SAVE IN DATA FILE
save_fname=fullfile(PRM.RAT_STORAGE_DIR,sprintf('%s_DATA.mat',PRM.RAT_NO));
save(save_fname,'STIM','ALL_DATA','ChanMap','LF_Shank','PRM');

fprintf('\n Rat Data for %s saved in %s',PRM.RAT_NO,save_fname)

clear ('ALL_DATA','LFP','STIM','LF_Shank','PRM');

end


























