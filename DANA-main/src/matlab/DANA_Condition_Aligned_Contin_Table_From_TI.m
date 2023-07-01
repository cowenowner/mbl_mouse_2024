function [TBL, INFO] = DANA_Condition_Aligned_Contin_Table_From_TI(M,data_labels, TI,rat_num, ALLCV, varargin)
% M is multi-col matrix where col 1 is time (seconds) and the other cols
% are the data that are to be aligned.
% data_labes is a cell array of the names of each column
LV = unique(TI.LV_group);
Hz = unique(round(TI.Hz));
Hz = Hz(Hz < 50);

sec_before = 35;
sec_after = 20;

Extract_varargin;

FSCV_LV = unique(floor(ALLCV.LV*10)/10);
FSCV_HZ = unique(floor(ALLCV.Hz));

sFreq = 1/median(diff(M(:,1)));

INFO.FSCV_x_sec = ALLCV.x_sec;

INFO.x_axis_s = [];
ALL_M = []; ALL_Mnorm = []; ALL_Hz_G = []; ALL_LV_G = []; AFTIX = []; BEFIX = [];V = [];
cnt = 1;
for iG = 1:length(Hz)
    for iLV = 1:length(LV)
        IX = TI.Hz_group == Hz(iG) & TI.LV_group == LV(iLV);
        IXFSCV = ALLCV.Hz == Hz(iG) & floor(ALLCV.LV*10)/10 == LV(iLV);
        FSCV = ALLCV.M(IXFSCV,:);
        for iC = 2:Cols(M)
            %         PETH_raster(T_uS{iN}/100,TI.end_times_sec(ix)*10000,20,time_before_sec*1000,time_after_sec*1000);
            % due to drift and other things, probably would be good to
            % subtract pre stim baseline.
            [M2, ~ , x_axis_s] = PETH_EEG_simple(M(:,[1 iC]),TI.end_times_sec(IX),sec_before*sFreq, sec_after*sFreq, sFreq,false);
            % [FSCV, ~, x_sec_FSCV] = PETH_EEG_simple(FSCV_DATA, NP.Stim_Start_sec(IX), sec_before*sFreq, sec_after*sFreq, sFreq,false);
            start_s = -10;
            NORMIX = x_axis_s > start_s-21 & x_axis_s < start_s - 1
            M2n = M2-mean(M2(:,NORMIX),2);
            % figure; imagesc(x_axis_s,[], M2)
            % figure; imagesc(x_axis_s,[], M2n)
            if isempty(INFO.x_axis_s)
                INFO.x_axis_s = x_axis_s;
            end


            for iTrial = 1:Rows(M2)


                V(cnt).RatNum = rat_num;
                V(cnt).UniqueID = rat_num*10000 + iC-1*100 + iTrial;
                V(cnt).TrialNum = iTrial;
                V(cnt).Hz = Hz(iG);
                V(cnt).LV = LV(iLV);

                V(cnt).data_type = data_labels{iC-1};
                V(cnt).data = M2(iTrial,:);
                V(cnt).FSCV = FSCV(iTrial,:);

                ix = find(IX);
                % the actual stimulation sequence.
                V(cnt).Stim_sequence_sec = TI.Stim_sequence_sec{ix(1)}';


                cnt = cnt + 1;
            end

        end

    end

end

if isempty(V)
    TBL = [];
else
    TBL = struct2table(V);
    TBL.data_type = categorical(TBL.data_type);
end
% Do as above, but with ensemble activity.

