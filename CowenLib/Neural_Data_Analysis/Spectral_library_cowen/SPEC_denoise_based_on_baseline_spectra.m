function [denoised_signal,fft_sig_to_denoise] = SPEC_denoise_based_on_baseline_spectra(baseline_sig, sig_to_denoise, sFreq, threshold )
%% The threshold is the hard part - you need to figure this out.
% TODO: probably a good idea to pad the start and end of the sig_to_denoise
% to get rid of edge artifacts.
% Cowen 2024
what_was_added = [];
if length(baseline_sig) ~= length(sig_to_denoise)
    error('Currently this dumb function requires the length of the noise and teh output to be the same. Sorry.')
end
if nargin == 0
    x = 0:.01:2*pi;
    baseline_sig = sin(x*20) + .5*cos(x*20+.2).^2;
    x2 = 2:.01:20*pi;
    x2 = x;
    sig_to_denoise = sin(x2*20) + .5*cos(x2*20+.2).^2 + 2*sin(x2*4+.2);
    sig_to_denoise = sig_to_denoise + .4*randn(size(sig_to_denoise));
    what_was_added = 2*sin(x2*4+.2);
    sFreq = 1/.01;
    threshold = 3;
end
n =  length(baseline_sig); % order of the fft.
fft_baseline = fft(baseline_sig,n);
fft_sig_to_denoise = fft(sig_to_denoise, n); % you have to ensure the order of the fft is the same so that it fits to the signal.
% Identify the noise frequencies (e.g., frequencies below a certain threshold)
noise_freq_IX = abs(fft_baseline) >= threshold;
% noise_freq_indices = find(abs(fft_signal) <= threshold); % sometimes this
% could be appropriate??
fft_filter = fft_baseline;
fft_filter(noise_freq_IX) = 0;
fft_sig_to_denoise_post_filter = fft_sig_to_denoise;
fft_sig_to_denoise_post_filter(noise_freq_IX) = 0;

% Compute the inverse Fourier transform
denoised_signal = ifft(fft_sig_to_denoise_post_filter,n,'symmetric');
% upscale it again so that it matches the origianl data.
% x = linspace(1,length(sig_to_denoise),length(baseline_sig));
% denoised_signal_up = interp1(x, denoised_signal, 1:length(sig_to_denoise));


if nargout ==0 
    figure
    subplot(2,2,1)
    plot(baseline_sig)
    title('orig signal')
    subplot(2,2,2)
    plot(sig_to_denoise)
    hold on
    if ~isempty(what_was_added)
        plot(what_was_added)
    end
    title('sig to denoise')
    subplot(2,2,3)
    hold on 
    plot(abs(fft_baseline))
    plot(abs(fft_sig_to_denoise)+5)
    plot(abs(fft_filter)+15)
    plot(abs(fft_sig_to_denoise_post_filter)+25)
    legend('base','orig sig', 'filt','postfilt' )
    title('fft of the signals')
    subplot(2,2,4)
    hold on
    plot(denoised_signal)
    % plot(denoised_signal_up)
end


