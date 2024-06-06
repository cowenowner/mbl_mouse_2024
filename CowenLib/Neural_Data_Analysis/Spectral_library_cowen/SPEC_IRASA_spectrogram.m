function [Spec] = SPEC_IRASA_spectrogram(L, sFreq, Frange, WinSizeSec, FrangeForPowerLaw)
% FROM: https://purr.purdue.edu/publications/1987/1
% Also see https://pubmed.ncbi.nlm.nih.gov/36639900/
% INPUT: 
% A single column vector signal.
% sFreq= the sampling rate
% Frange - 2 element vector indicating the min and max frequencies analyzed
% WinSizeSec - function chunks over these windows (like pwelch). Choose a
% reasonable value to get at the frequency that you are interested.
% FrangeForPowerLaw - often you will not want some very low or high
% frequencies for power law fitting as strong slow oscillations or plateaus
% will corrupt the measure. This range allows you to focus the range to a
% region where you can be confident of a decent power law estimate.
% Defaults to the Frange
%
% This was used in the Halje paper on psychedelics. I think the results
% look good and it gets rid of the 1/f nicely.
%
% Run with no outputs to get plots
%
% note, your frequency resolution will depend a bit on your window size -
% smaller window, the lower the frequency resolution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen - wrapper for the amri_ functions from link above.
% 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    Frange = [10 sFreq/4]; % this method cannot go up to sFreq/2 so you may need to sample at a higher rate for these frequecnies.
end
if nargin < 4
    WinSizeSec = 8; % used in a Halje paper on psychedelics
end
if nargin < 5
    FrangeForPowerLaw = Frange;
end
WinSamp = round(WinSizeSec*sFreq);
% NOTE: it is best to chop up your continuous data in to chunks.
L = L(:);
n_res = floor(length(L)/WinSamp);
Lr = reshape(L(1:(n_res*WinSamp)),[WinSamp n_res]);
Spec = amri_sig_fractal(Lr,sFreq,'frange',Frange);
Spec = amri_sig_plawfit(Spec,FrangeForPowerLaw);
% save space as this can get pretty big.
Spec.freq = single(Spec.freq);
Spec.osci = single(Spec.osci);
Spec.frac = single(Spec.frac);
Spec.mixd = []; % For now, I don't see us neededing this so remove to save a lot of space.
Spec.Plaw = []; % For now, I don't see us neededing this so remove to save a lot of space.

if nargout == 0
    % Display the spectra in log-log scale
    SPEC_IRASA_spectrogram_plot(Spec);

    figure
    subplot(2,1,1)
    PW = pwelch(L, WinSizeSec*sFreq, WinSizeSec*sFreq/2, Frange(1):.2:Frange(2), sFreq);
    plot(Frange(1):.2:Frange(2),log10(PW))
    title('Pwelch for comparison')
    % legend('Frac.osci','pwelch'); legend boxoff
    subplot(2,1,2)
    spectrogram(L, WinSizeSec*sFreq, WinSizeSec*sFreq/2, Frange(1):.2:Frange(2), sFreq,'yaxis')
    title('spectrogram')
end