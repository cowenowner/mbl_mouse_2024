function [new_SP] = DANA_process_spikes(SP, quality_threshold)
if nargin == 0
    load('AllSpikes.mat','SP')
end
min_spikes = 20;
if nargin < 2
    quality_threshold = 2;
end
% NOTE - I think people are storing the tip depth using the wrong units.
% Correct with: [SP(1).depth_of_probe_tip_uM]/1000 - 
T_uS = []; All_T = []; WV = [];gix = [];
cnt = 1;
for ii = 1:length(SP)
    if isnan(SP(ii).Quality) || (SP(ii).Quality >= quality_threshold && length(SP(ii).t_uS) > min_spikes)
        gix(cnt) = ii;
        T_uS{cnt} = SP(ii).t_uS;
        All_T = [All_T; SP(ii).t_uS(:)];
        [~,ix] = max(max(SP(ii).WV.mWV,[],1));
        WV(cnt,:) = SP(ii).WV.mWV(:,ix)';
        SP(ii).WaveformBest = SP(ii).WV.mWV(:,ix);
        % Correct incorrect depth.
        tip_dpth = SP(1).depth_of_probe_tip_uM;
        if tip_dpth >1000 && tip_dpth < 9000
            % Probably OK as is
        else
            tip_dpth = tip_dpth/1000;
        end

        SP(ii).Depth_corrected_uM = tip_dpth - SP(ii).depth_on_probe;
        cnt = cnt + 1;
    end
end
% [SP(1).depth_of_probe_tip_uM tip_dpth]
new_SP = SP(gix);