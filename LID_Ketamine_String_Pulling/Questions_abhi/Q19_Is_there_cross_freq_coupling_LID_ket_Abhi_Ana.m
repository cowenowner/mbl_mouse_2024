%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q19_Is_there_cross_freq_coupling_LID_ket_Abhi_Ana';
ses_to_ana = 'Q19_Is_there_cross_freq_coupling_LID_ket_Abhi';

GP.Analysis_Dir = 'E:\Temp_6.6.23';

PLOT_IT = true;
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
cm = lines(5);

% Initialize
SES = [];
ses_cnt = 1;

for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));

    % Get the structure
    cTBL = Dset.cTBL;
    
    % Creating meta table for phase and power modulation information foe
    % each neuron
    if Dset.aborted == false
        for ii = 1:length(cTBL)
            % general
            
            SES(ses_cnt).Rat = Dset.SES.Rat;
            SES(ses_cnt).Session = Dset.SES.Session;
            if Dset.SES.RatType == '6ODHA_LID'
                SES(ses_cnt).Group = '6OHDA_LID';
            elseif Dset.SES.RatType == 'SHAM'
                SES(ses_cnt).Group = 'SHAM';
            end
            SES(ses_cnt).LFP_filename = cTBL(ii).LFP_filename;
            SES(ses_cnt).Drugs = Dset.SES.Drugs;
            SES(ses_cnt).Interval = cTBL(ii).Interval;
            SES(ses_cnt).Depth_lf = cTBL(ii).Depth;
            if cTBL(ii).Depth < 2800
                SES(ses_cnt).BrainRegion = 'M1';
            else
                SES(ses_cnt).BrainRegion = 'Striatum';
            end
            
            SES(ses_cnt).Tetrode = cTBL(ii).TT;
            if cTBL(ii).TT < 9
                SES(ses_cnt).Hemisphere = 'L';
            elseif cTBL(ii).TT > 8
                SES(ses_cnt).Hemisphere = 'R';
            end
            
            SES(ses_cnt).method = cTBL(ii).method;
            SES(ses_cnt).CM = cTBL(ii).CM;
            SES(ses_cnt).low_fq = cTBL(ii).low_fq_range;
            SES(ses_cnt).high_fq = cTBL(ii).high_fq_range;
            

            ses_cnt = ses_cnt + 1;
            
        end
    else
    end
end
% save to table
TBL = struct2table(SES);

%% Plot stuff
low_fq = TBL.low_fq(1,:);
high_fq = TBL.high_fq(1,:);
method = TBL.method(1);

IX_LK_LS = categorical(TBL.Drugs) == 'LDOPA&Ketamine' | categorical(TBL.Drugs) == 'LDOPA&Saline';
IX_R_80 = IX_LK_LS & categorical(TBL.Hemisphere) == 'R' & TBL.Interval == 2;
temp = TBL(IX_R_80,:);
for ii = 1:sum(IX_R_80)
    figure
    imagesc(low_fq, high_fq, temp.CM{ii,1})
    ylabel('High fq')
    xlabel('Low fq')
    title(sprintf('Rat %d Session %d depth %.2f',temp.Rat(ii), temp.Session(ii), temp.Depth_lf(ii)))
    axis xy
    colorbar
end



