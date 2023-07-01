function Downsample_LFP(DATA_DIR)

% determine resampling rate for LFP
LFP_sFreq = 1500;

% Read info.rhd to get original sampling rate
IF = INTAN_Read_RHD_file(fullfile(DATA_DIR,'info.rhd')); %IF will contain meta data on the Intan session
fs_initial = IF.frequency_parameters.amplifier_sample_rate;


% Save LFP data if required...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir(fullfile(DATA_DIR,'LFP_sfreq1500'))
ff = find_files(fullfile(DATA_DIR,'amp-*.dat'));
if isempty(ff)
    % Assume that we are re-running the code and the aux files have been
    % zipped up. So, unzip them.
    disp('Unzipping AMP Files')
    unzip(fullfile(DATA_DIR,'amp-B-013.dat.zip'),DATA_DIR); % CHanGE THIS
    unzip(fullfile(DATA_DIR,'amp-B-054.dat.zip'),DATA_DIR); % CHanGE THIS
    ff = find_files(fullfile(DATA_DIR,'amp-*.dat'));
end



nrecs = [];
LFP = [];
cnt = 1;
for iF = 1:length(ff)
    [~,fname] = fileparts(ff{iF});
    fp = fopen(ff{iF},'rb');
    D = fread(fp,'int16');
    fclose(fp);
    LFP.data = int16(resample(D,LFP_sFreq, fs_initial));
    LFP.to_uV_conversion = 0.195; % to microvolts
    LFP.LFP_sFreqj = LFP_sFreq;
    LFP.original_sFreq = fs_initial;
    LFP.fname = ff{iF};
    
    nrecs(cnt) = length(D);
    save(fullfile(DATA_DIR,'LFP_sfreq1500',[fname '_LFP']),'LFP')
    fprintf('.')
    cnt = cnt + 1;
end
% Save timestamps for each record in the data.
if min(abs(diff(nrecs))) > 2
    error('rec size not equal among files')
end
n_lfp_recs = length(LFP.data);
t_uS = 1e6*((0:(n_lfp_recs-1))/LFP_sFreq + .5/LFP_sFreq);

save(fullfile(DATA_DIR,'LFP_sfreq1500','LFP_times'),'t_uS')

disp('Saved LFP data')

delete('*.dat')


