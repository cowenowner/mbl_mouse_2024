function OUT = Q9_Does_neural_activity_corr_with_treadmill()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% this was tested on rat 320 day 3. C:\Users\Stephen Cowen\Box\Cowen Laboratory\Data\LID_Ketamine_Single_Unit_R56\Rat320\03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_IT = true;
SES = LK_Session_Info();
OUT = [];
[GP,E,EVT,SP,RHD,META,DEPTHS,TS] = LK_Load_Important_Things();
TRD = LK_Treadmill_Speed(EVT, META);
POS = LK_Load_and_Clean_POS;
IMU = LK_Load_and_Process_IMU;
% Let's determine if the IMU data correlates in some way to the treadmill
% data...
imu_in_trd = interp1(IMU.t_uS,IMU.speed,TRD.t_uS);
figure
IX1 = TRD.tread_direction == 1;
IX2 = TRD.tread_direction == -1;
plot(TRD.speed(IX1),imu_in_trd(IX1),'.',TRD.speed(IX2),imu_in_trd(IX2),'.')
lsline
xlabel('tread speed');ylabel('IMU speed');
% Let's determine if the PC of the spike data has some relationship with
% treadmill speed....
TSr = Restrict(TS,TRD.tread_intervals_uS);
bin_ms = 20;
% Mex file for Bin_ts_array does not work on mac. Load data from Abhi
% instead.
[Q,Qt_uS] = Bin_ts_array(TSr,bin_ms*1000);
%load('matlab.mat');
bin_t_uS = mean(Qt_uS,2);
Qf = conv_filter(Q,hanning(5));
[~,sc,lat] = pca(Qf);

scintrd = interp1_matrix(bin_t_uS,sc(:,1:2),TRD.t_uS);

% Sarah and Abhi.. Let's determine if the FR of the spike data has some 
% relationship with treadmill speed....

qfintrd = interp1_matrix(bin_t_uS,Qf(:,:),TRD.t_uS); % interpolate Qf data to match speed. 
%trdinqf = interp1_matrix(TRD.t_uS,TRD.speed,bin_t_uS);  % interpolate speed data to match Qf

% Find average firing rate for speed of interest during speed increase
% There is some noise and the speed might slightly decreasing at
% some time sampled during an interval when the speed is increasing 
% Plot this to see what I am talking about: 
% figure; plot(TRD.speed(393:456)); hold on; plot(TRD.speed(261:324)); 
% Blue curve shows decrease.
% To get around this, smooth the time series to find the intervals 
trdspeed = movmean(TRD.speed,5);
trdspeed_accel = cat(1, NaN,diff(trdspeed)); % accelation at time t, is the speed increases from the previous time point. Add NaN at beginning because speed before start is unknown. 
vect_speed_increasing = (trdspeed_accel > 0); % find time points when speed increases (aka difference in speed between time point and the previous time is positive)

indices_increasing = nan(1,2);
indices_decreasing = nan(1,2);

prevVal = NaN;
for idx = 1:length(vect_speed_increasing)
   thisVal = vect_speed_increasing(idx);
   if idx == 1
       if thisVal == 1 
           indices_increasing(end+1,1) = idx;
       else
           indices_decreasing(end+1,1) = idx;
       end
   elseif thisVal ~= prevVal
        if thisVal == 1 
           indices_decreasing(end,2) = idx-1;
           indices_increasing(end+1,1) = idx;
        else
           indices_increasing(end,2) = idx-1;
           indices_decreasing(end+1,1) = idx;
       end
       
   elseif idx == length(vect_speed_increasing)
      if thisVal == 1 
          indices_increasing(end,2) = idx;
      else
          indices_decreasing(end,2) = idx;
      end
   end
    
   prevVal = thisVal;
end

% Now delete any windows that are 1 sample or nan
indices_increasing(indices_increasing(:,1)-indices_increasing(:,2)==0|all(isnan(indices_increasing),2),:) = [];
indices_decreasing(indices_decreasing(:,1)-indices_decreasing(:,2)==0|all(isnan(indices_decreasing),2),:) = [];

% Had made a matrix call speed change intervals so combine indices into a
% single matrix so you can use the logic you already wrote

speedChangeIntervals = cat(1, cat(2, ones(size(indices_increasing,1),1), indices_increasing), cat(2, zeros(size(indices_decreasing,1),1), indices_decreasing));
% SpeedChangeIntervals is a matrix that has the start and end index from 
% vect_speed_increasing for any intervals where there is a block of consecutive 0s or 1s 
% The first column indices whether the block was 1 or zero (1 = increasing
% speed). The second column contains the start index and the third is the
% end index. 

% There is some noise and the speed might slightly decreasing at
% some time sampled during an interval when the speed is increasing 
% Plot this to see what I am talking about: 
% figure; plot(TRD.speed(393:456)); hold on; plot(TRD.speed(261:324)); 
% Blue curve shows decrease.
% To fix this, we want to loop through the speed change intervals and
% combine any blocks of zeros of ones that are separate by less than 5(chose this aribtrarily for now)
newSpeedChangeIntervals = [];
prevInterval = [];

[~,idx] = sort(speedChangeIntervals(:,2));
speedChangeIntervals = speedChangeIntervals(idx,:);
speedChangeInterval_duration = speedChangeIntervals(:,3) - speedChangeIntervals(:,2); % subtract start index from end index to determine how long the block lasts (there is some 

% Now find the intervals that are shorted than 20 samples and see if you
% need to combine it with the interval before or after
idxSmall = find(speedChangeInterval_duration < 20);
idxToDelete = [];
for idx = 1:length(idxSmall)
    if idx == 1 % if first interval
        curr = speedChangeIntervals(idxSmall(idx));
        prev = NaN;
        next = speedChangeIntervals(idxSmall(idx) + 1);
    elseif idx == length(speedChangeIntervals) % if last interval
        curr = speedChangeIntervals(idxSmall(idx));
        prev = speedChangeIntervals(idxSmall(idx) - 1);
        next = NaN;
    else
        curr = speedChangeIntervals(idxSmall(idx));
        prev = speedChangeIntervals(idxSmall(idx) - 1);
        next = speedChangeIntervals(idxSmall(idx) + 1);
    end
    if curr == prev && curr ~= next
        speedChangeIntervals(idxSmall(idx)-1, 3) = speedChangeIntervals(idxSmall(idx), 3); % set the interval before this one to include this smaller interval
        idxToDelete = cat(1, idxToDelete, idxSmall(idx));
    elseif curr== next && curr ~= prev
        speedChangeIntervals(idxSmall(idx)+1, 2) = speedChangeIntervals(idxSmall(idx), 2); % set the interval after this one to include this smaller interval
        idxToDelete = cat(1, idxToDelete, idxSmall(idx));
    elseif curr == next && curr == prev
        
    end
end

speedChangeIntervals(idxToDelete,:) = []; % Now delete the rows that you combined with the row above or below
speedChangeIntervals(speedChangeIntervals(:,3)-speedChangeIntervals(:,2)<40,:) = []; % now delete any intervals that are less than 40 samples long

% Now add the direction as a column (indices from speeed are the same indices as for direction)
direction = cellfun(@(x,y) mean(TRD.tread_direction(x:y,1)), num2cell(speedChangeIntervals(:,2)),num2cell(speedChangeIntervals(:,3)),'uniformOutput',1);

speedChangeIntervals(:,end+1) = direction;
% Alright, now we can look at FR for increasing and decreasing speed
% intervals for CW and CCW
FR_CWIncreasingSpeed = speedChangeIntervals(speedChangeIntervals(:,1)==1 & speedChangeIntervals(:,4)==1,:);
FR_CCWIncreasingSpeed = speedChangeIntervals(speedChangeIntervals(:,1)==1 & speedChangeIntervals(:,4)==-1,:);
FR_CWDecreasingSpeed = speedChangeIntervals(speedChangeIntervals(:,1)==0 & speedChangeIntervals(:,4)==1,:);
FR_CCWDecreasingSpeed = speedChangeIntervals(speedChangeIntervals(:,1)==0 & speedChangeIntervals(:,4)==-1,:); % first col = increasing decreasing 4th col = avg direction over window


FRCWIncrease = {};
FRCCWIncrease = {};

FRCWDecrease = {};
FRCCWDecrease = {};



for i = 1:size(qfintrd,2)
    
    thisFRCWIncrease = cellfun(@(x,y) qfintrd(x:y,i)', num2cell(FR_CWIncreasingSpeed(:,2)),num2cell(FR_CWIncreasingSpeed(:,3)),'uniformOutput',0);
    thisFRCWDecrease = cellfun(@(x,y) qfintrd(x:y,i)', num2cell(FR_CWDecreasingSpeed(:,2)),num2cell(FR_CWDecreasingSpeed(:,3)),'uniformOutput',0);
    thisFRCCWIncrease = cellfun(@(x,y) qfintrd(x:y,i)', num2cell(FR_CCWIncreasingSpeed(:,2)),num2cell(FR_CCWIncreasingSpeed(:,3)),'uniformOutput',0);
    thisFRCCWDecrease = cellfun(@(x,y) qfintrd(x:y,i)', num2cell(FR_CCWDecreasingSpeed(:,2)),num2cell(FR_CCWDecreasingSpeed(:,3)),'uniformOutput',0);


    FRCWIncrease = cat(2, FRCWIncrease, thisFRCWDecrease);
    FRCWDecrease = cat(2, FRCWDecrease, thisFRCWIncrease);
    FRCCWIncrease = cat(2, FRCCWIncrease, thisFRCCWDecrease);
    FRCCWDecrease = cat(2, FRCCWDecrease, thisFRCCWIncrease);
    
end

% The intervals are not consistent sizes so for now, just append nan to end
% so we can plot easily
maxlength = max(max(cellfun(@numel, cat(1, FRCWIncrease, FRCWDecrease,FRCCWIncrease,FRCCWDecrease)))); % find max interval length
FRCWIncrease = cellfun(@(x) [x, nan(1,maxlength-numel(x))], FRCWIncrease, 'UniformOutput', false);
FRCWDecrease = cellfun(@(x) [x, nan(1,maxlength-numel(x))], FRCWDecrease, 'UniformOutput', false);
FRCCWIncrease = cellfun(@(x) [x, nan(1,maxlength-numel(x))], FRCCWIncrease, 'UniformOutput', false);
FRCCWDecrease = cellfun(@(x) [x, nan(1,maxlength-numel(x))], FRCCWDecrease, 'UniformOutput', false);

% Plot all cells
for idx = 1:size(FRCWIncrease,2)
    Cell_FRCWIncrease = cell2mat(FRCWIncrease(:,idx));
    Cell_FRCWDecrease = cell2mat(FRCWDecrease(:,idx));
    Cell_FRCCWIncrease = cell2mat(FRCCWIncrease(:,idx));
    Cell_FRCCWDecrease = cell2mat(FRCCWDecrease(:,idx));
    figure;
    subplot(2,2,1);
    imagesc(Cell_FRCWIncrease);
    xlabel('time');
    ylabel('bout');
    title(sprintf('Neuron %d\n Clockwise  : Increasing Speed',idx));
    subplot(2,2,2);
    imagesc(Cell_FRCWDecrease);
    xlabel('time');
    ylabel('bout');
    title(sprintf('Neuron %d\n Clockwise  : Decrease Speed',idx));
    
    subplot(2,2,3);
    imagesc(Cell_FRCCWIncrease);
    xlabel('time');
    ylabel('bout');
    title(sprintf('Neuron %d\n Counter-Clockwise  : Increasing Speed',idx));
    subplot(2,2,4);
    imagesc(Cell_FRCCWDecrease);
    xlabel('time');
    ylabel('bout');
    title(sprintf('Neuron %d\n Counter-Clockwise  : Decrease Speed',idx));
    
end

%% LFP beta and treadmill speed and direction analysis
% load LFP
fqs = 1:.5:150;
load('Processed_Data/best_channels_gamma_80.mat','OUT')
if exist(fullfile('LFP',OUT.best_non_reref),'file')
    % load if stored locally.
    LFP = LK_Load_and_Clean_LFP('./LFP',OUT.best_non_reref);
else
    global DIRS
    lfp_dir = fullfile(DIRS.LFP_Dir,SES.rat_str,SES.session_str,'LFP');
    LFP = LK_Load_and_Clean_LFP(lfp_dir, OUT.best_non_reref);
end

filts = SPEC_create_filters({'beta' [3 6.5]}, LFP.sFreq);
treadmill_start_min = E.MinFromStart(E.EventID == 'TreadmillStart');
treadmill_win = [0 30] + treadmill_start_min;
treadmill_base_win = [-25 -5] + treadmill_start_min;
IXeff = LFP.t_uS >treadmill_win(1)*60e6 & LFP.t_uS < treadmill_win(2)*60e6;

% Create an index of beta power.
L = zeros(length(LFP.LFP),length(filts));
for ii = 1:length(filts)
    L(:,ii) = filtfilt(filts{ii},LFP.LFP);
end
%powbeta = abs(hilbert(L(:,1))) - abs(hilbert(L(:,2)))./ abs(hilbert(L(:,1))) + abs(hilbert(L(:,2)));
powbeta = (abs(hilbert(L(:,1))));

pow_low_th = abs(hilbert(L(:,2)));

%%
if PLOT_IT
    figure
    subplot(1,3,1)
    pwelch(LFP.LFP(IXeff),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
    grid off
    pubify_figure_axis
    title('beta during treadmill win')
    subplot(1,3,2)
    pwelch(LFP.LFP(IXbase),LFP.sFreq,LFP.sFreq/2,fqs,LFP.sFreq)
    grid off
    pubify_figure_axis
    title('prior to treadmill')
    subplot(1,3,3)
    cla
    pwelch(powbeta(IXeff),LFP.sFreq,LFP.sFreq/2,1:.2:40,LFP.sFreq)
    hold on
    pwelch(powbeta(IXbase),LFP.sFreq,LFP.sFreq/2,1:.2:40,LFP.sFreq)
    
    grid off
    pubify_figure_axis
    h = gca;
    h.Children(1).Color = [0 0 0 ];
    title('psd of power changes at beta')
    
end

figure
plot(LFP.t_uS(IXeff)/1e6, L(IXeff,1))
hold on
plot(LFP.t_uS(IXeff)/1e6, powbeta(IXeff))
plot(LFP.t_uS(IXeff)/1e6, pow_low_th(IXeff))

plot(LFP.t_uS(IXeff)/1e6, LFP.LFP(IXeff))
legend('fbeta','pbeta','pthet','orig');
%% We probably need to find the start and end time of each
% interval better than what I did here, but this is good enough to at least
% visualize the data. 
% Also looks like the firing rate was not recorded during the middle
% chunk of bouts (?) need to look into this more. 
figure
subplot(3,2,1:2)
imagesc(bin_t_uS,[],Qf')
subplot(3,2,3:4)
plot(bin_t_uS,sc(:,1),bin_t_uS,sc(:,2))
axis tight
subplot(3,2,5)
%plot(TRD.speed(IX1),scintrd(IX1,1),'.',TRD.speed(IX2),scintrd(IX2,2),'.')
plot(TRD.speed(IX1),scintrd(IX1))
lsline
xlabel('tread speed');ylabel('sc spikes');
