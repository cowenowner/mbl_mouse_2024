function [fqs, new_t] = instfrq_cowen(x,t,fr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An alternative to the built in instfq function. It was originally motivated by
% apparent low time resolution of instfrq - but later I found that you
% could overcome this with pspectrum. See below.
%
% NOTE: if you want to increase the time resolution of instfrq(), see this...
% https://www.mathworks.com/matlabcentral/answers/474062-how-to-change-time-interval-value-when-using-instfreq-function-in-matlab
%
%
% INPUT: provide, for example, a wavelet spectrum (row = time, col = fq)
% and the frequencies for each col. This just finds the fq of max power for
% each time point. This function is only as good as the norm_spec passed in
% so be very careful that it reflects data that is valid for computing
% oscillatory power (e.g., is of sufficient amplitude and does not have the
% 1/f component in it.
%
%
% 
% Cowen 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    fr = 40;
end
[p,fd,td] = pspectrum(x,t,'spectrogram','FrequencyResolution',fr); 
% the last argument is to customize the time interval
% you get the range of argument when you execute the command if it exceeds the range
[fqs,new_t]=instfreq(p,fd,td); 