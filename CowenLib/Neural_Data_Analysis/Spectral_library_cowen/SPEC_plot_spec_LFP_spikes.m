function SPEC_plot_spec_LFP_spikes(t, CWT, fqs, pk_fqs, LFP, ST, time_divisor, time_label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% t = time for each col in CWT 
% CWT = wavelet spectrogram - one point for each point in LFP and t. Row =
%          Hz, Col = time.
% fqs = frequencies for each row in CWT.
% pk_fqs = the instantaneous frequency for each col in CWT. Make unknown
%          fqs into nans.
% LFP = LFP signal - 
% ST = cell array of spike times in same units of t.
% time_divisor = divide time units by this (e.g., to convert into seconds)
% time_label = text for the x label
%
% inputs can come from the following functions
% [pk_fqs] = instfrq_cwt_cowen(CWT,psd_fqs,[78 98]);
% [CWT,phase] = SPEC_cwt_simple_cowen(LFP.data_uV, LFP.sFreq, psd_fqs);
%
% OUTPUT: plot that can be scrolled. 
% The last 'neuron' is the aggregate of all neurons.
%
% NOTE: the scrolling often crashes on matlab version 2022 and earlier.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ST = Restrict(ST, t(1), t(end));

ST2 = cell(size(ST));
ALLST = [];
for ii = 1:length(ST)
    ST2{ii} = ST{ii}/time_divisor;
    ALLST = [ALLST;ST2{ii}(:)];
end
ST2{end+1} = unique(ALLST);

figure
a(1) = subplot(6,1,1:2);
plot(t/time_divisor,LFP,'color',[.7 .7 .7]); axis tight
% ylabel('uV')
yyaxis right
plot_raster(ST2,[],'make_colorful',false)
box off
ylabel('Neuron')
% set(gca,'XTickLabel',[])

a(2) = subplot(6,1,3:4);
imagesc(t/time_divisor,fqs(:),CWT);
colormap(jet)
axis xy
ylabel('Hz')
hold on
plot(t/time_divisor,pk_fqs,'.m')
box off
% set(gca,'XTickLabel',[])

a(3) = subplot(6,1,5:6);
plot(t/time_divisor,LFP); axis tight
box off
% ylabel('uV')

linkaxes(a,'x')
xlabel(time_label)

