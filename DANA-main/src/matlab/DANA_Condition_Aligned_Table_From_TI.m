function [TBL, INFO] = DANA_Condition_Aligned_Table_From_TI(SP,TI,rat_num, varargin)
Hz = [10 20];
LV = [0 .3 1 1.2];
base_time_sec = -11;
peth_bin_ms = 550;
time_before_sec = 25;
time_after_sec = 50;

Extract_varargin;

T_uS = {SP.t_uS}; WV = [SP.WaveformBest]';


INFO.PETH_x_axis_ms = [];
ALL_M = []; ALL_Mnorm = []; ALL_Hz_G = []; ALL_LV_G = []; AFTIX = []; BEFIX = [];V = [];
cnt = 1;
for iG = 1:length(Hz)
    for iLV = 1:length(LV)
        IX = TI.Hz_group == Hz(iG) & TI.LV_group == LV(iLV);
        for iN = 1:length(T_uS)
            %         PETH_raster(T_uS{iN}/100,TI.end_times_sec(ix)*10000,20,time_before_sec*1000,time_after_sec*1000);
            % due to drift and other things, probably would be good to
            % subtract pre stim baseline.
            [M, x_axis_ms, A_msec, h ] = PETH_raster(T_uS{iN}/100,TI.end_times_sec(IX)*10000,peth_bin_ms,time_before_sec*1000,time_after_sec*1000);
            BEFIX = (x_axis_ms/1000) > base_time_sec-14 & (x_axis_ms/1000) < base_time_sec;
            
            if isempty(INFO.PETH_x_axis_ms)
               INFO.PETH_x_axis_ms = x_axis_ms;
            end

            ALL_M(cnt,:) = mean(M);
            ALL_Mnorm(cnt,:) = mean(M - mean(M(:,BEFIX),2));
            ALL_Hz_G(cnt) = Hz(iG);
            ALL_LV_G(cnt) = LV(iLV);

            V(cnt).RatNum = rat_num;
            V(cnt).CellID = rat_num*100 + iN;
            V(cnt).CellIDWithinRat = iN;
            V(cnt).WV = WV(iN,:);
            V(cnt).Depth_corrected_uM = SP(iN).Depth_corrected_uM;
            ix = find(IX);
            % the actual stimulation sequence.
            V(cnt).Stim_sequence_sec = TI.Stim_sequence_sec{ix(1)};

            V(cnt).mean_PETH = mean(M);
            V(cnt).mean_PETHnorm = mean(M - mean(M(:,BEFIX),2));
            V(cnt).Hz = Hz(iG);
            V(cnt).LV = LV(iLV);

            % Compute LV for each interval.
            IVsec{1} = [TI.end_times_sec(IX)-25 TI.end_times_sec(IX)-11];
            IVsec{2} = [TI.end_times_sec(IX)+1 TI.end_times_sec(IX)+15];
            IVsec{3} = [TI.end_times_sec(IX)+15 TI.end_times_sec(IX)+29];
            IVsec{4} = [TI.end_times_sec(IX)+29 TI.end_times_sec(IX)+33];
            
            % Compute LV, CV for each interval
            tmpLV = nan(1,length(IVsec));
            tmpCV = nan(1,length(IVsec));
            %
            V(cnt).PeakISIs = nan;
            V(cnt).PeakRatesHz = nan;

            for iI = 1:length(IVsec)
                t = Restrict(T_uS{iN},IVsec{iI}*1e6);
                d = diff(t/1000);
                d = d(d<5000);
                if length(d) > 20
                    [~,tmpLV(iI) ] = LocalVariance(d,5000);
                    tmpCV(iI) = std(d)/mean(d);
                    BA = Burst_analysis(cumsum(d));
                    if ~isempty(BA.PeakISIs)
                        V(cnt).PeakISIs = BA.PeakISIs(1);
                        V(cnt).PeakRatesHz  = BA.PeakRatesHz(1);
                    end
                end

                % TBL.BurstHistISI(iR,:) = BA.HistISI;

            end
            V(cnt).LocalVarianceAroundStim = tmpLV;
            V(cnt).CVAroundStim = tmpCV;

            cnt = cnt + 1;
        end
        
    end

end
if isempty(V)
    TBL = [];
else
    TBL = struct2table(V);
end
% Do as above, but with ensemble activity.

