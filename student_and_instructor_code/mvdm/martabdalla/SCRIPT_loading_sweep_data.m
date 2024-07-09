%%
restoredefaultpath;
cd('C:\Users\mattm\Documents\GitHub\mbl_mouse_2024\wavesurfer');
installWavesurferForNow; % need this to use wavesurfer data loader

cd('C:\data\nsb2024\MartAbdalla');
fn = 'M6_20240706_baseline_0001.h5';

%%
data = loadDataFile(fn);

firstSweepNumber = data.header.NextSweepIndex;
lastSweepNumber = data.header.NextSweepIndex + data.header.NSweepsPerRun - 1; % not sure if this always works

Fs = data.header.AcquisitionSampleRate;

%% get one sweep
[this_sweep_tvec, this_sweep_data] = getSweepData(data, firstSweepNumber);

%% get all sweeps and concatenate
all_data = []; all_tvec = [];
for iS = firstSweepNumber:lastSweepNumber % will be slow because not pre-allocating variables

    [this_sweep_tvec, this_sweep_data] = getSweepData(data, iS);

    all_data = cat(1,all_data,this_sweep_data);
    all_tvec = cat(1,all_tvec,this_sweep_tvec);

end 


%% TODO: get artifact size from a given sweep


%% TODO: get baseline from a given sweep

%%

function [tvec, data_out] = getSweepData(data_in, sweepNumber)
% retrieve timestamps and data for specified sweep

field_no = cat(2, 'sweep_0', num2str(sweepNumber)); % build field name of requested sweep
data_out = data_in.(field_no).analogScans;

nSamples = size(data_out, 1); % build matching vector of timestamps
t0 = data_in.(field_no).timestamp;

Fs = data_in.header.AcquisitionSampleRate;
tvec = (0:nSamples-1)'*(1/Fs) + t0;

end