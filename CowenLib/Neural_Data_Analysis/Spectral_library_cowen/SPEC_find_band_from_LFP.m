function [bandHz, peakHz, fqs, PW] = SPEC_find_band_from_LFP(L, sFreq, target_rangeHz, wider_rangeHz )
% Looks for a peak in the PSD and finds a reasonable window around that
% peak.
% Cowen 2023.
fqs = wider_rangeHz(1):.1:wider_rangeHz(2);
PW = pwelch(L,round(sFreq/4),[],fqs,sFreq); %  pburg(L,141,fqs,sFreq);
PW = log10(PW);
GIX = fqs<target_rangeHz(1) | fqs > target_rangeHz(2);
m = fitlm(fqs(GIX)',PW(GIX)');
PW = PW(:) - predict(m,fqs(:));
[pv,ix] = findpeaks(PW);
new_fqs = fqs(ix);

gix = new_fqs > target_rangeHz(1) & new_fqs < target_rangeHz(2);
if sum(gix) == 0 || max(pv(gix)) <=0 % If < 0, that means its < than baseline so NO peak in traget so ignore
    peakHz = nan;
    bandHz = [nan nan];
    disp('no peak found in target_range.')
else
    new_fqs = new_fqs(gix); 
    pv = pv(gix);
    [mx,ix] = max(pv);
    peakHz = new_fqs(ix(1));
    GIX = fqs>target_rangeHz(1) & fqs < target_rangeHz(2);
    PWs = PW(GIX);
    fqs_s = fqs(GIX);
    ix = find_intervals(PWs,mx/2,mx/4); % Cowen - added mean - seems reasonable - but needs verifiation.
    if isempty(ix) || ix(end) == length(PWs) % the last condition is when there is no clear peak as the right edge abuts the target_rangeHz  
        peakHz = nan;
        bandHz = [nan nan];
        disp('no peak found in target_range.')
    else
        bandHz = fqs_s(ix);
    end
    % figure;plot(fqs,PW); hold on; plot_markers_simple(bandHz)
end
