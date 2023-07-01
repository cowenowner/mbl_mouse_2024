function [objective]  = DA_LV_opt(ISI, lv_tgt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize so that local variance is varied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isi_thresh = 0.003; % could have this as an input.
v = (ISI(1:(end-1))-ISI(2:end))./(ISI(1:(end-1))+ISI(2:end));
LV = 3*sum(v.^2)/(length(ISI)-1);
% penalize any negative ISIs or ISIs shorter than some interval.
negpen = double(sum(ISI<isi_thresh))*1.9;
objective = (LV-lv_tgt).^2 + negpen;
