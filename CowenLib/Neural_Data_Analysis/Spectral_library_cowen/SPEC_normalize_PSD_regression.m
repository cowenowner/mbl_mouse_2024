function [PW] = SPEC_normalize_PSD_regression(fqs,PSD,control_bands_Hz)
% Looks for a peak in the PSD and finds a reasonable window around that
% peak.
% Cowen 2023.

if Rows(PSD) > 1
    PW = zeros(Rows(PSD),length(fqs));
    for iR = 1:Rows(PSD)
        PW(iR,:) = SPEC_normalize_PSD_regression(fqs,PSD(iR,:), control_bands_Hz);
    end
    return
end
GIX = false(size(fqs));
for iR = 1:Rows(control_bands_Hz)
    GIX = GIX |  fqs >= control_bands_Hz(iR,1) & fqs <= control_bands_Hz(iR,2);
end
m = fitlm(fqs(GIX)',PSD(GIX)');
PW = PSD(:)- predict(m,fqs(:));

% figure;plot(fqs,PSD);hold on;plot(fqs,PW)
