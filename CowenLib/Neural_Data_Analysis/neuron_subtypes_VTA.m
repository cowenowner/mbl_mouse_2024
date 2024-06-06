function [G,INFO] = neuron_subtypes_VTA(spikes_uS_ca, WV_IN, WV_sFreq, varargin)
% function [G,INFO] = neuron_subtypes_VTA(spikes_uS_ca, WV_IN, WV_sFreq, varargin)
%
% Features (e.g., peak to trough width, half width, amplitude ratios, will
% be computed from the data you provide so provide good data.
%
% For the Hill et al paper, fixing the random seed for clustering to 2 (rng
% (2) - to help with replication.
%
% General procedure: Step 1) Cluster unsupervised based on waveform and
% waveform parameters like half width and peak to trough distance.
% Cluster into k groups (2-3 typically).
% Step 2) Once done, take the group with the
% largest mean waveform and pull out neurons that have firing rates above a
% reasonable amount since we know that DA neurons fire at lower rates. Also
% take out the subset of these neurons that also have waveform widths less
% than the global median. These two extracted groups then get their own
% categories and are typically not analyzed as they are small groups and
% their cell type is rather ambiguouss.
% All narrow waveforms that are identified in stage 1 are analyzed and go
% to the Narrow category which really means narrow waveform cells.
%
%
% INPUT: cell array of action potential times for each neuron in uS.
%        WV_IN = one row for each neuron - the waveform. Presumed to be
%        aligned on the max or min (we found it best to align on the max of
%        the ABSOLUTE value of the waveform).
%        Since the polarity of waveforms can invert if you are at the
%        initial segment vs. dendrites, it MAY be optimal to invert the
%        waveform if the max is > min of the waveform. Also assumes that
%        the waveform is in some units (normalized - but specify if not uV).
%        WV_sFreq = the sampling frequency of the passed in waveform.
%
% OUTPUT: G clusters
%         INFO more data on each cluster like mean rate and stuff.
%
% Things to read:
%
% Ungless, M.A., Grace, A.A., 2012. Are you or aren t you? Challenges associated with physiologically identifying dopamine neurons. Trends Neurosci 35, 422–430. https://doi.org/10.1016/j.tins.2012.02.003
%  Argues for bi or tri-phasic. Notch sometimes. 2-10Hz. >2.2ms duration
%  overall. the tri-phasic is only really seen with high-pass filters at
%  50000Hz. We don't have this. Data is band-pass filtered I believe to
%  upper limit of 8000Hz at most. Probably 6000Hz and samppled at 30000Hz
%  argues that filter settings matter - a high pass of 300Hz (typicall) can
%  mask parts of the DA waveform.
%  NOTE: Indeed, the often-reported triphasic dopamine neuron action potential is actually reflective of a filtering artefact, since with open filters dopamine neuron action potentials are biphasic when recorded extracellularly (see Box 1b) [39].
%  NOTE: A more conservative criterion of ≥1.1 ms from start to trough of the waveform was found to be adequate to identify putative dopamine neurons even with filtering [24]. Indeed, applying this criterion to invivojuxtacellular labelling studies shows that ≈ 90% of neurons fitting this criterion are dopaminergic [24–28, 47] (Table 1).
%  From my initial pass at our data - applying this criterion yields 35
%  neurons out of 572 - this seems unlikely. It does bug me though and
%  makes me want to not call anything a DA neuron explicitly, just WWLF
%  (for wide waveform low firing).
%   -- doing this we could be more agnostic. WWLF, WWHF (omit likely due to
%   small group), NW (no need to split into subgroups as we'll analyze FR
%   later.)
%

% The following labeling intracellular + extra seem to confirm that low
% firing and wide waveform is sufficient to ID da cells. More nuances to
% this are discussed in teh Ungless paper above so this is debated.
% Brown, M.T.C., Henny, P., Bolam, J.P., Magill, P.J., 2009. Activity of neurochemically heterogeneous dopaminergic neurons in the substantia nigra during spontaneous and driven changes in brain state. J Neurosci 29, 2915–2925. https://doi.org/10.1523/JNEUROSCI.4423-08.2009
% Pinault, D., 1996. A novel single-cell staining procedure performed in vivo under electrophysiological control: morpho-functional features of juxtacellularly labeled thalamic cells and other central neurons with biocytin or Neurobiotin. J Neurosci Methods 65, 113–136. https://doi.org/10.1016/0165-0270(95)00144-1
%
% Bakkum, D.J., Obien, M.E.J., Radivojevic, M., Jäckel, D., Frey, U., Takahashi, H., Hierlemann, A., 2019. The Axon Initial Segment is the Dominant Contributor to the Neuron’s Extracellular Electrical Potential Landscape. Adv Biosyst 3, e1800308. https://doi.org/10.1002/adbi.201800308
%  important as it shows that the polarity of the AP is completely
%  dependent on electrode position AND that largest component of all APs
%  (not just DA) is the potential in the axon initial segment.
%
% López-Jury, L., Meza, R.C., Brown, M.T.C., Henny, P., Canavier, C.C., 2018. Morphological and Biophysical Determinants of the Intracellular and Extracellular Waveforms in Nigral Dopaminergic Neurons: A Computational Study. J. Neurosci. 38, 8295–8310. https://doi.org/10.1523/JNEUROSCI.0651-18.2018
%  Talks about the notch in the AP of DA cells. They also admit there is a lot of variability here so a clear algorithm for this is not very obvious.
%  Was a modeling study of SNc DA cells
%  FROM PAPER: There is strong support in the literature for the two-component shape of the intracellular and extracellular AP in DA neurons (Grace and Bunney, 1983a, b; Kuhr et al., 1987; Bean, 2007; Henny et al., 2012; Tucker et al., 2012).
%  Not reliably seen in extracellular waveform:  extracellular AP (Brown et al., 2009; Meza et al., 2018).
%  Does a nice simulation of the extracellular AP. HUGE variability.
%  BUT - you do see that strong biphasic pattern typically. Perhaps a good
%  measure would be both a large amplitude ratio AND a large peak-to-trough
%  and/or half width. The problem with this is that there is a NEGATIVE
%  relationship in our data between the amp ratio and the width which is
%  confusing as both are supposed to be indicators of DA cells.
%
% Roesch MR, Calu DJ, Schoenbaum G. Dopamine neurons encode the better option in rats deciding between differently delayed or sized rewards. Nature Neuroscience. 2007;10:1615–1624.
%  this is the paper we originally used - it's older and I am not a huge
%  fan of the amplitude ratio as it seems to pull out a small subset of
%  cells that I think were possibly inverted. Later studies seemed to do a
%  better job of validating with intracellular - and also verifying that
%  meassures of width were good predictors.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_IT = true;
VALIDATE_IT = false;
% Key parameters for clustering...
USE_FIXED_RNG = true; % For the hill paper we used rng seed of 2 (see rng)
USE_WV_ENERGY = true; % Most recent addition and potentially most optional: Used for the hill paper as it indicates flux in the waveform and biphasicness. This is our standard measure in spike sorting in general so seems like a good feature.
USE_WV_ENERGY_DIFF = true; % measures the degree of fluctuation in the waveform- should highlight multi-phasic neurons
USE_ACORR = false; % For DA neuron identification we did not use ACORR
USE_AMP_RATIO = false; % For DA neuron identification we did not use amp ratio
FLIP_BIG_POSITIVE = true; % True for Hill paper as sometimes the peak was > trough. If not flipped, the estimation of spike width was messed up and this could be due to inversion based on electrode placement relative to the axon initial segment.

JUST_PLOT_DA_NARROW = true;
critical_range_ms = [-.35 .9]; % range for clustering (aside from features like half width etc..). Helpful as it limits to a reasonable range of the waveform shape.
AC_binsize_ms = 4; % While the AC may not be used, we do use it for visualization.
DA_neuron_FR_cutoff_Hz = 11; % 11Hz was the max found in a set of papers so seems reasonable.
DA_neuron_width_percentile_cutoff = 50; % Once you have an initial set of clusters, use this to omit narrow waveforms from the wide waveform putative DA cluster
AC_duration_ms = 200;
WV_units = 'uV'; % whatever the units are that are passed in. We do assume non-normalized data are pased in.
max_k = 8; % max number of clusters.
Cat_groups = {'pDA' 'Narrow' 'WideHROther' 'Ambig'}; % only the first 2 are used in the paper.
% Thinking about renaming to WWLF NW WWHF to be more agnostic and not make
% strong claims about DA neurons. I think now that we have an indication of
% increased bursting in these cells, we can make stronger claims though.
Cat_colors = [ 0.4940    0.1840    0.5560; 0.4660    0.6740    0.1880; 0 1 1 ;0     0     0];
Cat_mkr = {'o' 'o' 's' 's'};

Extract_varargin

if USE_FIXED_RNG
    rng(2) % ensure we get consistent results. Just do this once things are working iff you want to reproduce a plot.
end

% x axis of the waveform.
WV_x_ms = (0:(Cols(WV_IN)-1)) * 1000/WV_sFreq;
% Align x axis to the peak - so 0 = peak.
[~,ix] = max(mean(abs(WV_IN),'omitmissing'));
WV_x_ms = WV_x_ms - WV_x_ms(ix);
% Indices of region we will use for clustering.
WV_IX = WV_x_ms >= critical_range_ms(1) & WV_x_ms <= critical_range_ms(2);
wv_st_ix = find(WV_IX); 
wv_st_ix = wv_st_ix(1);
% Flip the sign if necessary. A positive peak > the trough indicates an
% inverted waveorm. If so, invert it.
if FLIP_BIG_POSITIVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This seems quite reasonable, but it is not typical that I am aware of
    % so perhaps should be omitted. If we omit though - then we REALLY need
    % to combine more than one group into DA neurons as the flipped ones
    % will cluster right away.
    % I tried without this and wow - BIG change in clustering.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mx = max(WV_IN,[],2,'omitmissing');
    minv = min(WV_IN,[],2,'omitmissing');
    IX = mx > minv*-1;
    WV_IN(IX,:) = WV_IN(IX,:)*-1;
    disp(['Flipped ' num2str(sum(IX)) ' waveforms.'])
end

% WV_IN(isnan(WV_IN)) = 0; % fill nan with zero as assume this is before the start or after the end of the waveform.
mn = mean(WV_IN(:,1:wv_st_ix),2,'omitmissing');
mn(isnan(mn)) = 0;
WV_IN = WV_IN - mn;
WV2 = WV_IN;
% % For each waveform, find the last point and smooth it down to zero.
end_ix = Cols(WV_IN);
for iR = 1:Rows(WV_IN)
    start_nan_ix = find(isnan(WV_IN(iR,:)) & WV_x_ms > .1, 1,'first');
    if ~isempty(start_nan_ix)
        fill_ix = start_nan_ix:min([end_ix start_nan_ix+10]);
        WV2(iR,fill_ix) = interp1([start_nan_ix-1 fill_ix(end)],[WV_IN(iR,start_nan_ix-1) 0],fill_ix,'linear');
    end

end
WV_IN = WV2;
%
WV_IN(isnan(WV_IN)) = 0;


% rescale so that the distance from the first point to the minumum is zero.
WV_norm_to_min = WV_IN./abs(min(WV_IN,[],2,'omitmissing'));

% Determine the width of each waveform (peak to trough etc...)figuffasdf
% Mult by -1 since Spike_width assumes upward spike.
[peak_half_width_ms, peak_to_trough_ms, WVup] = Spike_width(WV_norm_to_min*-1,10);
WVup = WVup*-1; % unused for now. upsampled.
% Another measure of half width could be the half width of the absolute
% value of the waveform as this would handle the bi-phasic waves and give
% them a wide waveform (which they should have).
% IX  = (peak_to_trough_ms > 27 & peak_to_trough_ms < 28)
peak_half_width_ms = 1000*( peak_half_width_ms/WV_sFreq);
peak_to_trough_ms = 1000*( peak_to_trough_ms/WV_sFreq);

if USE_AMP_RATIO

    % The ungless paper suggested start to trough of 1.1 ms. A more conservative criterion of ≥1.1 ms from start to trough of the waveform was found to be adequate to identify putative dopamine neurons even with filtering [24]. Indeed, applying this criterion to invivojuxtacellular labelling studies shows that ≈ 90% of neurons fitting this criterion are dopaminergic [24–28, 47] (Table 1).
    % what is the 'start'. We have peak to trough so just add start to peak and
    % we're done.
    sd = max(std(abs(WV_IN(:,1:5)),'omitmissing'));
    mn = max(mean(abs(WV_IN(:,1:5)),'omitmissing'));
    V = abs(WV_IN) > (mn+2*sd);
    for iV = 1:Rows(V)
        first_ix(iV) = find(V(iV,:),1,'first');
    end
    start_ms = WV_x_ms(first_ix)';
    start_to_trough_ms = abs(start_ms);
    % figure;imagesc(V)
    % histogram(start_to_trough_ms,100)

    % two measures of amplitude ratio: the first is the ratio of the max BEFORE
    % the dip. The second is the max AFTER the dip.
    % I am not really sold on this measure as it taps into the multi-phasic
    % features of DA neurons whih is somewhat unreliable. (see comments at top)

    min_v = abs(min(WV_IN,[],2));

    BEF_IX = WV_x_ms<0 &  WV_x_ms > critical_range_ms(1);
    % BEF_IX = WV_x_ms<0;
    max_v = max(WV_IN(:,BEF_IX),[],2);
    % amp_ratio_bef = (min_v-max_v)./(min_v+max_v);
    amp_ratio_bef = max_v./min_v;

    AFT_IX = WV_x_ms > 0;
    max_v = max(WV_IN(:,AFT_IX),[],2);
    % amp_ratio_aft = (min_v-max_v)./(min_v+max_v);
    amp_ratio_aft = max_v./min_v;
end
% NOTE: because of the flipping above, this will always be > 0. The
% positive values will not be larger than the dip values. This argues
% against this measure perhaps or perhaps not. I am still for it as we make
% the reasonable assumption that an AP needs to dip negative - if not, then
% the electrode was recording somewhere where the waveform was reversed.
% This may not be true all of the time, but it is safer give that we know
% the waveform does flip if you move from the axon initial segment to the
% dendrites.
%
Firing_Rate_Hz = nan(length(spikes_uS_ca),1);
LV = nan(length(spikes_uS_ca),1);
histisi_x_log_ms = linspace(log10(1),log10(3000),500);
ac = []; 
histisi_ks_logms = nan(length(spikes_uS_ca),length(histisi_x_log_ms));
histisi_ks_logms = nan(length(spikes_uS_ca),length(histisi_x_log_ms));
histisi_logms = nan(length(spikes_uS_ca),length(histisi_x_log_ms));
for iN = 1:length(spikes_uS_ca)
    [tmp,ac_lag_x_ms] = AutoCorr(spikes_uS_ca{iN}/1e3, AC_binsize_ms, AC_duration_ms/AC_binsize_ms); % 100x faster than Cross_cor. Same result.
    ac(iN,:) = movmean(tmp,3);
    % Same thing but SLOOOWER [tmp2,ac_lag_x_ms2] = Auto_corr(spikes_uS_ca{iN}/1e3, AC_binsize_ms, AC_duration_ms/AC_binsize_ms); % 100x faster than Cross_cor. Same result.

    % tmp = histcounts(log10(diff(spikes_uS_ca{iN}/1e3)),x_log);
    d = diff(spikes_uS_ca{iN}/1000);
    d(d>2000) = []; d(d<1) = [];
    [histisi_ks_logms(iN,:)] = ksdensity(log10(d),histisi_x_log_ms);
    [histisi_logms(iN,:)] = histcounts(log10(d),[histisi_x_log_ms histisi_x_log_ms(end)+log10(1)]);

    duration_sec = (spikes_uS_ca{iN}(end)- spikes_uS_ca{iN}(1))/1e6;

    Firing_Rate_Hz(iN) = length(spikes_uS_ca{iN})/duration_sec;
    [~,LV(iN)] = LocalVariance(diff(spikes_uS_ca{iN})/1e3);
end

ac_norm = ac./sum(ac,2);
ac_norm(isnan(ac_norm)) = 0; % sometimes no spikes so divide by zero error.

if VALIDATE_IT
    figure;plot_neurons(spikes_uS_ca, WV_IN, WV_x_ms);
    figure;plot_neurons(spikes_uS_ca, WV_norm_to_min, WV_x_ms);
    figure;plot(peak_half_width_ms,peak_to_trough_ms,'.'); xlabel('HW ms'); ylabel('PT ms')
    figure;plot(amp_ratio_bef,peak_to_trough_ms,'.'); xlabel('amp_ratio_bef'); ylabel('PT ms')
    figure;plot(amp_ratio_bef,amp_ratio_aft,'.'); xlabel('amp_ratio_bef'); ylabel('amp_ratio_aft')
    figure;
    subplot(2,2,1)
    tsne_plot([WV_norm_to_min(:,WV_IX)]);
    subplot(2,2,2)
    tsne_plot([WV_norm_to_min(:,WV_IX) Z_scores([peak_half_width_ms peak_to_trough_ms])]);
    subplot(2,2,3)
    tsne_plot([WV_norm_to_min(:,WV_IX) peak_half_width_ms peak_to_trough_ms ac_norm]);
    subplot(2,2,4)
    % tsne_plot([ALL_WV_norm(:,WV_IX)]);
    tsne_plot([WV_IN(:,WV_IX)]);
    tsne_plot([peak_half_width_ms peak_to_trough_ms ac_norm]);
    tsne_plot([amp_ratio_bef amp_ratio_aft peak_half_width_ms peak_to_trough_ms ]);
    tsne_plot([amp_ratio_bef amp_ratio_aft peak_to_trough_ms ]);
end

% FINAL DECISION: NOTE: adding amp_ratio and peak2trough2 seemed to
% mess thigns up.
M = [peak_half_width_ms peak_to_trough_ms WV_norm_to_min(:,WV_IX) diff(WV_norm_to_min(:,WV_IX),[],2) ];
% Clustering on the log - results in essentially only one cluster.
% M = [log10(peak_half_width_ms) log10(peak_to_trough_ms WV_norm_to_min(:,WV_IX) diff(WV_norm_to_min(:,WV_IX),[],2) ];

if USE_ACORR
    M = [M ac_norm];
end
if USE_AMP_RATIO
    M = [M amp_ratio_bef amp_ratio_aft];
end
energy = sum(WV_norm_to_min(:,WV_IX).^2,2);
energy_diff = sum(diff(WV_norm_to_min(:,WV_IX),[],2).^2,2);
if USE_WV_ENERGY
    M = [M energy];
end
if USE_WV_ENERGY_DIFF
    M = [M energy_diff];
end


M = Z_scores(M);
% Hierarchical clustering optimal K
klist=2:max_k;%the number of clusters you want to try
eva = evalclusters(M,'linkage','CalinskiHarabasz','klist',klist); % there are other measures like silhouette
% eva = evalclusters(M,'linkage','gap','klist',klist); % I get a lot of clusters with GAP
C = clusterdata(M,'Linkage','ward','SaveMemory','off','Maxclust',eva.OptimalK);
% C = clusterdata(M,'Linkage','ward','SaveMemory','off','Maxclust',3);
% If there is a cluster with < thresh members, then merge it to the next
% closest in terms of peak width.

CNTS = histcounts(C);

[PCA_Mpc,PCA_Msc,PCA_Mlat] = pca(M);

mn_hw = grpstats(peak_half_width_ms,C,{'mean'});
mn_pt = grpstats(peak_to_trough_ms,C,{'mean'});

[~,tmpDA_ix] = max(mn_hw);
[~,tmpDA_ix2] = max(mn_pt);
DA_ix = unique([tmpDA_ix tmpDA_ix2]);

if length(DA_ix) > 1
    error('more than one group is likely dopamine: TODO- accomodate this in the code. Hope this is rare')
end

if DA_ix ~= 1
    % Make sure the DA neurons are in group 1
    tmp_C = C;
    % Flip the 2 groups
    C(C==DA_ix) = 1;
    C(tmp_C==1) = DA_ix;
end
% Add one more group: Dopamine neurons that EXCEED the threshold.
BIX = C==1 & Firing_Rate_Hz > DA_neuron_FR_cutoff_Hz;
C(BIX) = 0; % a new cluster for now that has these highFRwidewave cells.

th_hw = prctile(peak_half_width_ms,DA_neuron_width_percentile_cutoff);
th_pt = prctile(peak_to_trough_ms,DA_neuron_width_percentile_cutoff);

BIX = C==1 & peak_half_width_ms < th_hw & peak_to_trough_ms < th_pt;
C(BIX) = -1; % a new cluster for now that has these highFRwidewave cells.

% figure;hist(C)

% redo group stats.
mn_hw = grpstats(peak_half_width_ms,C,{'mean'});
mn_pt = grpstats(peak_to_trough_ms,C,{'mean'});
mn_fr = grpstats(Firing_Rate_Hz,C,{'mean'});
mn_lv = grpstats(LV,C,{'mean'});
% Assign proper group names.
G = categorical(repmat({'Unknown'},length(spikes_uS_ca),1));

G(C==1) = Cat_groups{1};
%let's go through the other groups now.
u = unique(C);
% For now, we're going to combine all Non-DA into one group.
% This could change if we find evidence that we should further sub-divide
% the OTHER group.
other_cats = u(u > 1);
for iG = 1:length(other_cats)
    IX = C == other_cats(iG);
    G(IX) = Cat_groups{2};
    % G(IX) = sprintf('Narrow%d', iG);
end

G(C==0) = Cat_groups{3};
G(C==-1) = Cat_groups{4};
% Do a final cleanup of few narrow cells that have a WIDE half width since that is strange and qualifies as 'ambiguous'.
BIX = G =='Narrow' & peak_half_width_ms > .3;
G(BIX) = Cat_groups{4}; % to catch the few narrows that had a VERY wide half width.

G = removecats(G); % cleans unused categories.
% G = reordercats(G,Cat_groups); % put DA first.

% Some basic stats...
G1IX = G =='pDA';
G2IX = G =='Narrow';

[h,FR_p] = ttest2(Firing_Rate_Hz(G1IX),Firing_Rate_Hz(G2IX));
FR_Cohens_d = Cohens_d(Firing_Rate_Hz(G1IX),Firing_Rate_Hz(G2IX));

[h,peak_to_trough_p] = ttest2(peak_to_trough_ms(G1IX),peak_to_trough_ms(G2IX));
peak_to_trough_Cohens_d = Cohens_d(peak_to_trough_ms(G1IX),peak_to_trough_ms(G2IX));

[h,LV_p] = ttest2(LV(G1IX),LV(G2IX));
LV_Cohens_d = Cohens_d(LV(G1IX),LV(G2IX));



% Plot.
if PLOT_IT
    % just DA other
    % NEED TO FINALIZE THE PLOT IN THE FIGURE!!!

    if JUST_PLOT_DA_NARROW
        cats_to_plot = Cat_groups(1:2);
    else
        cats_to_plot = Cat_groups;
    end
    GIX = false(size(G));
    for iG = 1:length(cats_to_plot)
        GIX = GIX | G == cats_to_plot{iG};
    end

    figure;
    subplot(3,2,1);
    GG = G(GIX);
    GG = removecats(GG);

    histogram(GG,'FaceColor','k', 'FaceAlpha',1); pubify_figure_axis; ylabel('Count')

    subplot(3,2,2)
    mn = [];
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        if any(IX)
            R = Firing_Rate_Hz(IX);
            mn(ii) = mean(R);
            [ks,x] = ksdensity(R,[0:.2:30], 'Support'  ,'positive');
            plot(x,ks,'Color',Cat_colors(ii,:),'LineWidth',4)
            % histogram(R,0:2:26,'Normalization','probability','EdgeColor',clrs(ii,:),'FaceAlpha',.1, 'LineWidth',4)
            hold on
        end
    end
    axis tight
    for ii = (length(cats_to_plot)):-1:1
        plot_vert_line_at_zero('where_to_plot', mn(ii),'color',Cat_colors(ii,:))
    end
    % legend(Cat_labs(1:end-1));legend boxoff
    xlabel('Hz');ylabel('p')
    axis tight
    pubify_figure_axis
    % xlim([0 30])

    subplot(3,2,3)
    % Boxplot_points(Firing_Rate_Hz(GIX),GG)
    Boxplot_points(log10(Firing_Rate_Hz(GIX)),GG)
    ylabel('log10 Hz')
    pubify_figure_axis


    subplot(3,2,4)
    Boxplot_points(log10(peak_to_trough_ms(GIX)),GG)
    ylabel('log10 peak_to_trough_ms')
    pubify_figure_axis

    subplot(3,2,5)
    Boxplot_points(log10(peak_half_width_ms(GIX)),GG)
    ylabel('log10 peak_half_width_ms')
    pubify_figure_axis

    subplot(3,2,6)
    Boxplot_points(LV(GIX),GG)
    ylabel('LV')
    pubify_figure_axis

    % subplot(3,2,6)
    % Boxplot_points(amp_ratio_aft,G)
    % ylabel('amp_ratio_aft')

    % subplot(3,2,6)
    % Boxplot_points(amp_ratio_bef,G)
    % ylabel('amp_ratio_bef')


    set(gcf,'Position',[377    55   410   916])

    figure
    tsne_plot(M(GIX,:),'Groups',GG);
    xlabel('tsne dim 1'); ylabel('tsne dim 2')
    pubify_figure_axis

    sz = 8;

    figure
    for ii = (length(cats_to_plot)):-1:1
        IX = G == cats_to_plot{ii};
        % hh(ii) = plot(peak_half_width_ms(IX),peak_to_trough_ms(IX),mkr,'MarkerSize',sz,'Color',clr,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',clr);
        % hh(ii) = plot(peak_half_width_ms(IX),peak_to_trough_ms(IX),Cat_mkr{ii},'MarkerSize',sz,'Color',Cat_colors(ii,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Cat_colors(ii,:));
        hh(ii) = plot(log10(peak_half_width_ms(IX)),log10(peak_to_trough_ms(IX)),Cat_mkr{ii},'MarkerSize',sz,'Color',Cat_colors(ii,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Cat_colors(ii,:));
        hold on
    end
    % xlabel(' peak_half_width_ms ms'); ylabel(' peak_to_trough_ms ms')
    xlabel('log10 peak_half_width_ms ms'); ylabel('log10 peak_to_trough_ms ms')
    pubify_figure_axis



    figure
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        hh(ii) = plot3(peak_half_width_ms(IX),peak_to_trough_ms(IX),PCA_Msc(IX,1), Cat_mkr{ii},'MarkerSize',sz,'Color',Cat_colors(ii,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Cat_colors(ii,:));
        hold on
    end
    xlabel(' peak_half_width_ms ms'); ylabel(' peak_to_trough_ms ms'); zlabel('PC1')
    pubify_figure_axis

    if USE_AMP_RATIO

        figure
        for ii = (length(cats_to_plot)):-1:1
            IX = G==cats_to_plot{ii};
            % hh(ii) = plot(peak_half_width_ms(IX),peak_to_trough_ms(IX),mkr,'MarkerSize',sz,'Color',clr,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',clr);
            hh(ii) = plot(amp_ratio_bef(IX)+amp_ratio_aft(IX),peak_to_trough_ms(IX),Cat_mkr{ii},'MarkerSize',sz,'Color',Cat_colors(ii,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Cat_colors(ii,:));
            hold on
        end
        xlabel('amp ratio');ylabel('peak to trough')
        pubify_figure_axis
    end
    figure
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        % hh(ii) = plot(peak_half_width_ms(IX),peak_to_trough_ms(IX),mkr,'MarkerSize',sz,'Color',clr,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',clr);
        hh(ii) = plot(PCA_Msc(IX,1),peak_to_trough_ms(IX) + peak_half_width_ms(IX),Cat_mkr{ii},'MarkerSize',sz,'Color',Cat_colors(ii,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Cat_colors(ii,:));
        hold on
    end
    xlabel('PC1');ylabel('peak to trough + half width')
    pubify_figure_axis

    figure
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        % hh(ii) = plot(peak_half_width_ms(IX),peak_to_trough_ms(IX),mkr,'MarkerSize',sz,'Color',clr,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',clr);
        hh(ii) = plot(PCA_Msc(IX,1),PCA_Msc(IX,2),Cat_mkr{ii},'MarkerSize',sz,'Color',Cat_colors(ii,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',Cat_colors(ii,:));
        hold on
    end
    xlabel('PC1');ylabel('PC2')



    figure
    subplot(1,2,1)
    XIX = WV_x_ms > -.6 & WV_x_ms < .8;
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        wv = WV_norm_to_min(IX,XIX);
        plot(WV_x_ms(XIX),wv,':','Color',Cat_colors(ii,:))
        hold on
    end

    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        wv = WV_norm_to_min(IX,XIX);
        plot_confidence_intervals(WV_x_ms(XIX),wv,[],Cat_colors(ii,:))
        hold on
    end
    ylabel('norm to trough')
    xlabel('ms')
    pubify_figure_axis

    subplot(1,2,2)
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        wv = WV_norm_to_min(IX,XIX);
        plot_confidence_intervals(WV_x_ms(XIX),wv,[],Cat_colors(ii,:))
        hold on
    end
    ylabel('norm to trough')
    xlabel('ms')
    pubify_figure_axis
    set(gcf,'Position',[ 554                     363.5                     746.5                     371.5])

    figure
    subplot(1,2,1)
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        % sum(IX)
        plot_confidence_intervals(ac_lag_x_ms,ac(IX,:),[],Cat_colors(ii,:))
        hold on
    end
    ylabel('ac')
    xlabel('ms lag')
    pubify_figure_axis

    subplot(1,2,2)
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        % sum(IX)
        plot_confidence_intervals(ac_lag_x_ms,ac_norm(IX,:),[],Cat_colors(ii,:))
        hold on
    end
    ylabel('norm units')
    xlabel('ms lag')
    pubify_figure_axis
    sgtitle('Autocorrelogram')


    figure
    subplot(1,2,1)
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        % sum(IX)
        plot_confidence_intervals(histisi_x_log_ms,histisi_ks_logms(IX,:),[],Cat_colors(ii,:))
        hold on
    end
    ylabel('count')
    xlabel('log10 ms')
    pubify_figure_axis
    subplot(1,2,2)
    histisi_logms_norm = histisi_ks_logms./sum(histisi_ks_logms,2);
    for ii = (length(cats_to_plot)):-1:1
        IX = G==cats_to_plot{ii};
        % sum(IX)
        plot_confidence_intervals(histisi_x_log_ms,histisi_logms_norm(IX,:),[],Cat_colors(ii,:))
        hold on
    end
    ylabel('norm count')
    xlabel('log10 ms')
    pubify_figure_axis
    sgtitle('Hist ISI')

end

% Assumptions: cluster withe largest half width and/or peak-trough are big dopamine cells.
% Group 1 is da.

% Provide useful informaiton...
if nargout >1
    INFO.WV_norm_to_min = WV_norm_to_min;
    INFO.WV_x_ms = WV_x_ms;
    INFO.WVupsampled = WVup;
    INFO.firing_rate_Hz = Firing_Rate_Hz;
    INFO.LocalVariance = LV;
    INFO.peak_half_width_ms = peak_half_width_ms;
    INFO.peak_to_trough_ms = peak_to_trough_ms;
    INFO.energy = energy;
    INFO.energy_diff = energy_diff;
    
    INFO.ac = ac;
    INFO.ac_norm = ac_norm;
    INFO.ac_lag_x_ms = ac_lag_x_ms;
    INFO.FR_p = FR_p;
    INFO.FR_Cohens_d = FR_Cohens_d;
    INFO.LV_p = LV_p;
    INFO.LV_Cohens_d = LV_Cohens_d;
    INFO.peak_to_trough_p = peak_to_trough_p;
    INFO.peak_to_trough_Cohens_d = peak_to_trough_Cohens_d;
    INFO.histisi_ks_logms = histisi_ks_logms;
    INFO.histisi_logms = histisi_logms;
    INFO.histisi_x_log_ms = histisi_x_log_ms ;
    INFO.OptimalK = eva.OptimalK;

    INFO.USE_FIXED_RNG = USE_FIXED_RNG; % For the hill paper we used rng seed of 2 (see rng)
    INFO.USE_WV_ENERGY = USE_WV_ENERGY; % Most recent addition and potentially most optional: Used for the hill paper as it indicates flux in the waveform and biphasicness. This is our standard measure in spike sorting in general so seems like a good feature.
    INFO.USE_ACORR = USE_ACORR; % For DA neuron identification we did not use ACORR
    INFO.USE_AMP_RATIO = USE_AMP_RATIO; % For DA neuron identification we did not use amp ratio

    [INFO.PCA_WVpc,INFO.PCA_WVsc,INFO.PCA_WVlat] = pca(WV_norm_to_min(:,WV_IX));
end


if VALIDATE_IT
    figure;plot_neurons(spikes_uS_ca, WV_IN, WV_x_ms,'Groups',G,'SORT_ROWS',true);

    figure
    subplot(4,1,1)
    boxplot(log10(Firing_Rate_Hz),C)
    ylabel('LOG Firing rate Hz')
    subplot(4,1,2)
    boxplot(peak_half_width_ms,C)
    ylabel('Half WIdth')
    subplot(4,1,3)
    boxplot(LV,C)
    ylabel('Local Variance')
    subplot(4,1,4)
    hist(C)

    figure
    subplot(2,2,1)
    imagesc(ac_lag_x_ms,[],sort_matrix( ac)); colorbar
    subplot(2,2,3)
    plot_confidence_intervals(ac_lag_x_ms,ac)
    xlabel('lag ms')
    subplot(2,2,2)
    imagesc(ac_lag_x_ms,[],sort_matrix( ac_norm)); colorbar
    subplot(2,2,4)
    plot_confidence_intervals(ac_lag_x_ms,ac_norm)

end

% Check the waveforms
% Strange waveforms... 86 garbage 283? 315 is misaligned 353? 396 397 garbage
% Problem is that few return to baseline - always shifted up or something.
if VALIDATE_IT

    figure
    for iW = 1:Rows(WV_IN)
        subplot(2,2,1)
        plot(WV_x_ms, WV_IN(iW,:))
        hold on
        ylabel('dff')
        subplot(2,2,2)
        plot(ac_lag_x_ms, ac(iW,:))
        hold on
        subplot(2,2,3)
        plot(ac_lag_x_ms, ac_norm(iW,:))
        hold on
        subplot(2,2,4)
        plot(WV_x_ms, [0 diff(WV_IN(iW,:))],'-')
          hold on
        if mod(iW,10) == 0
            clf
        end
        title(num2str(iW))
        pause
    end
    sub_ix = [86 283 315 353 396 397];
    figure
    plot(WV_x_ms, WV_IN(sub_ix,:))
    legend
end
%
