function [wv,phase] = SPEC_cwt_simple_cowen(LFP, sFreq, fqs, varargin)
% Get a wavelet spectrogram, the 'simple' way.
% INPUT: LFP vector, sFreq, fqs: the desired frequencies you want to see.
% OUTPUT: wv - the power in each band over time, phase - phase in each
% band.
% phase- I am not as confident as this as there is TOO much phase alignment
% across bands - this should not be. Different frequences should be
% completely out of phase.
% 
% Cowen 2023
wavelet_type = 'bump';
voices = 32;
interp_meth = 'spline';
% interp_meth = 'linear';
% interp_meth = 'cubic'; % Never works as it expects uniform spacing between frequencies.

Extract_varargin;

[wt,wv_f] = cwt(double(LFP),wavelet_type, sFreq,'VoicesPerOctive',voices); %,'VoicesPerOctive',32

wv_f = wv_f(end:-1:1); 
wt = wt(end:-1:1,:);
wv = zeros(length(fqs),Cols(wt));

for iC = 1:Cols(wt)
    % NOTE: you many NEED to do a spline or logarithmic inerpolation to more
    % accurately estimate freqwuencies here.

    wv(:,iC) = interp1(wv_f, wt(:,iC),fqs(:),interp_meth);
end
if nargout > 1
    phase = angle(wv); % again - this seems like adjacent frequencies are too phase aligned, but maybe because cwt has more poor frequency resolution than filter hilber?
end
wv = abs(wv);

% wv = convn(wv,hanning(3)/sum(hanning(3)),'same'); % this fixes the discrete problem with indices.