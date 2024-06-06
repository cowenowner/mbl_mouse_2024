function [F, POW, ix_fq] = instfrq_cwt_cowen(CWT,fqs,desired_band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An alternative to the built in instfq function. I have been unhappy with
% the built in function - the frequency resolution can change radomly.
%
% INPUT: A wavelet spectrum (col = time, row = fq). Make sure it is smooth.
% I sometimes find that you get discrete frequency clusters if the CWT is
% say in singles or not smoothed. Smoothing with convn fixes this.
% frequencies for each col. This just finds the fq of max power for
% each time point. This function is only as good as the spec passed in. Do
% your best to get rid of the 1/f component beforehand (e.g., with a
% regression, baseline normalization, or something.)
% 
% fqs: the frequency for each row of CWT
%
% desired_band: upper and lower limit of the 'good' band. Peaks outside of
% this band will be thrown out.
%
% For example, here is my function for the CWT: 
%   [CWT,phase] = SPEC_cwt_simple_cowen(LFP.data, LFP.sFreq, psd_fqs);
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BAND_IX = fqs>=desired_band(1) & fqs<=desired_band(2);
[POW,ix_fq] = max(double(CWT));
F = fqs(ix_fq);

mnPOW = mean(mean(CWT(BAND_IX,:)));
sdPOW = std(mean(CWT(BAND_IX,:)));

% If it catches power outside of the core band, delete
% If the identified frequency was low power, delete.
BIX = (F <= desired_band(1)) | (F >= desired_band(2)) | (POW < mnPOW + sdPOW);

F(BIX) = nan; 
ix_fq(BIX) = nan;

if nargout ==0 
    figure
    imagesc([],fqs(:),CWT)
    hold on
    plot(F,'.w-')
    axis xy
end