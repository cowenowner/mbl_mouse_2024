%% Run Code of Single Pulse Stimulation Analysis/Slope Calculation

%% CHANGE THE FOLLOWING VALUES 

PRM_RAT_DATA_DIR='C:\SGL_DATA\vHPC_stim_mPFC_excitability_project\10907\';
PRM_SUB_FOLDER='mPFC_L5_bank0_g0';
PRM_StimFileName='mPFC_L5_bank0_g0_tcat.nidq.xa_5_0.txt';
PRM_HC_FirstLoc='vHC'; %Change this to vHC if vHC was first stimulated then the iHC
PRM_TotalStimCount=26; %CHange this to the total no. of stimulated variables

%load STIM FILE ORDER
load('C:\SGL_DATA\vHPC_stim_mPFC_excitability_project\10899\params_from_stim_exp_08_10_23_Rat10899.mat');

%Setup Data Dirs
PRM_ROOT_DATA_DIR=fullfile(PRM_RAT_DATA_DIR,PRM_SUB_FOLDER);
[~,root_folder] = fileparts(PRM_ROOT_DATA_DIR);
AP_FILE_DIR = fullfile(PRM_ROOT_DATA_DIR, [root_folder '_imec0']);

%Run CatGT for events only
PRM_EVENT_CHANNELS=[1 5];
NPXL_Extract_Events_With_CatGT_SS('PRM_ROOT_DATA_DIR', PRM_ROOT_DATA_DIR,...
     'PRM_EVENT_CHANNELS',PRM_EVENT_CHANNELS,...
        'PRM_INAROW',2,...
        'PRM_XA_THR1',0.8,...
        'PRM_XA_THR2',0.5);

%Load Stim Times
Stim_times_file=fullfile(PRM_ROOT_DATA_DIR,PRM_StimFileName);
Stim_times=readmatrix(Stim_times_file);
if PRM_HC_FirstLoc =='iHC'
    Stim_times_iHC=Stim_times(1:PRM_TotalStimCount); 
    Stim_times_vHC=Stim_times(PRM_TotalStimCount+1:2*PRM_TotalStimCount);
else
    Stim_times_vHC=Stim_times(1:PRM_TotalStimCount); 
    Stim_times_iHC=Stim_times(PRM_TotalStimCount+1:2*PRM_TotalStimCount);
end

%Change this when you want to run vHC vs iHC
Stim_times=Stim_times_vHC;

%Get Sample values from stim timestamps
NIDQ_meta_fname=strrep(PRM_SUB_FOLDER,'g0','g0_tcat.nidq.meta');
obj1 = SGLX_Class;
NIDQ.meta = obj1.ReadMeta(NIDQ_meta_fname,PRM_ROOT_DATA_DIR);

%Get meta files for LFP data
LFP_meta_fname=strrep(PRM_SUB_FOLDER,'g0','g0_t0.imec0.lf.meta');
LFP.meta=obj1.ReadMeta(LFP_meta_fname,AP_FILE_DIR);

%Initialize Slopes
Slope=zeros(length(Stim_times),385);

%% Calculate Slopes -Abhi Code
for stim =1:length(Stim_times)
    %Calculate sample times
    time_start=Stim_times(stim); %in s
    time_diff=0.1; %300 ms window
    
    LFP_Sample0=floor(time_start*str2double(LFP.meta.imSampRate));
    LFP_Sample_startafterstim = LFP_Sample0+8;
    LFP_nSamp=floor(1.0 * time_diff*str2double(LFP.meta.imSampRate));

    %AP_Sample0=str2double(AP.meta.firstSample)+(time_start*str2double(AP.meta.imSampRate));
    %AP_nSamp=time_diff*str2double(AP.meta.imSampRate);

    %Load LFP in sample range
    LFP.gain=250;
    LFP_binName=strrep(LFP_meta_fname,'.lf.meta','.lf.bin');
    LFP.data=obj1.ReadBin(LFP_Sample0,LFP_nSamp,LFP.meta,LFP_binName,AP_FILE_DIR);
    LFP.data_Volts=LFP.data.*(LFP.meta.imAiRangeMax/(LFP.meta.imMaxInt*LFP.gain));%in uV
    
    LFP.time=[time_start:(time_diff)/(LFP_nSamp-1):time_start+time_diff];%in ms
    LFP_sample_time=(LFP_Sample0+LFP_nSamp)/str2double(LFP.meta.imSampRate);

    All_LFP_data_V(stim,:,:)= LFP.data_Volts;
    
    %Plot post stim time
    %LFP.samples=1:LFP_nSamp;
    %plot(LFP.samples,LFP.data(1,:)+0.5*stim,'k');
    %xline(LFP.time(1,1),'r');
    %caxis([-100 100])
    
    %Find slope for each row
    [channel_cnt,~]=size(LFP.data_Volts);
    for ii=1:channel_cnt
        
        [~,min_idx]=min(LFP.data_Volts(ii,25:end));
        [~,max_idx]=max(LFP.data_Volts(ii,25:end));

        if min_idx > max_idx
            [~,min_idx] = min(LFP.data_Volts(ii,max_idx+24:end));
            min_idx = min_idx + (max_idx+23);
            max_idx = max_idx+24;
        elseif max_idx > min_idx
            [~,max_idx] = max(LFP.data_Volts(ii,min_idx+24:end));
            max_idx = max_idx + (min_idx+23);
            min_idx = min_idx+24;
        elseif max_idx == min_idx
            [~,max_idx] = max(LFP.data_Volts(ii,min_idx+24:end));
            max_idx = max_idx + (min_idx+23);
            min_idx = min_idx+24;
        end

        %Slope Calculation
        min_point=[LFP.time(min_idx),LFP.data_Volts(ii,min_idx)];
        max_point=[LFP.time(max_idx),LFP.data_Volts(ii,max_idx)];
        
        LFP_sec = LFP.time-LFP.time(1);

        min_point_all(stim,ii) = min_point(2); 
        min_time_all(stim,ii) = LFP_sec(min_idx);
        max_point_all(stim,ii) = max_point(2); 
        max_time_all(stim,ii) = LFP_sec(max_idx);
        Slope(stim,ii)=(max_point(2)-min_point(2))/(max_point(1)-min_point(1));
    end
end

%% Filter 60hz noise
%Remove 60hz Noise
All_LFP_data_V_filtered=zeros(size(All_LFP_data_V));

%notch filter
Fs=LFP.meta.imSampRate;
notch_filt= designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
%Filter
for stim_no =1: size(All_LFP_data_V,1)
    for chan_no=1:size(All_LFP_data_V,2)
        All_LFP_data_V_filtered(stim_no,chan_no,:)=filter(notch_filt,All_LFP_data_V(stim_no,chan_no,:));
    end
end

%% PLOT STUFF
%create video of raw lfp across different stim currents
[sorted_current,sorted_idx] = sort(PRM.output_current_array);
% x_s = LFP.time-LFP.time(1);
% XIX = x_s < 0.1;
figure
for ii = 1:26
    figure
    caxis([-1.5 1.5 ])
    per_stim(:,:) = All_LFP_data_V_filtered(sorted_idx(ii),:,:);
    imagesc(LFP.time-LFP.time(1),[],per_stim)
    colorbar('Limits',[-2 2])
    xlabel('sec')
    ylabel('channels')
    %     hold on
    plot_vert_line_at_zero(.05)
    %     pause(.4)
    title(sprintf('mPFC L5 vHPC Stimulation Voltage %0.2f',PRM.output_current_array(sorted_idx(ii))))
    axis xy
end

%% time of max response in the dorsal inframlimbic area

infra_idx = [140:1:180];
infra_idx2 = [50:1:150];
why = mean(max_time_all(:,infra_idx2),2);
why=why(sorted_idx);
figure
scatter(sorted_current, why)
ylabel('mean time to LFP response from stim (s)')
xlabel('vHPC stimulation current')
title('Response time and voltage relationship mPFC L23 vHPC stim')
pubify_figure_axis
pf = polyfit(sorted_current, why,3);
x1 = linspace(min(sorted_current), max(sorted_current),16);
pv = polyval(pf,x1);
hold on
plot(x1,pv)
lsline



figure
scatter(PRM.output_current_array, nanmean(max_point_all(:,infra_idx2),2))
ylabel('max LFP response')
xlabel('vHPC stimulation current')
title('Voltage and max LFP response relationship mPFC L23 vHPC stim Rat10899')
lsline
pubify_figure_axis

figure
scatter(PRM.output_current_array, nanmean(Slope(:,infra_idx2),2))
ylabel('Slope')
xlabel('vHPC stimulation current')
title('Voltage and Slope relationship mPFC L23 vHPC stim Rat10899')
lsline
pubify_figure_axis

%% Calculate 70% amplitude
max_x=3.68;
min_x=1.18;
x70=((max_x-min_x)*0.7)+min_x