function [AC,lags_sec] = Auto_corr_continuous(C,sFreq,sec_window)
% Cowen 2023
n_lags = round(sec_window*sFreq);
AC = xcorr(C,C,n_lags);
AC = AC((n_lags+2):end);
lags_sec = (1:n_lags)/sFreq;

if nargout == 0
    plot(lags_sec,AC,'LineWidth',4)
    xlabel('sec')
    pubify_figure_axis;

end
