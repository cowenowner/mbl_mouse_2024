%% population activity and firing rate analysis.

% Limit to M1 and neurons that fired at least a little.
GIX = sum(AllQ_all)>500 & TBL.Depth_uM' < 2000;

binsize_ms = 50; 
han_filt = hanning(20)/sum(hanning(20));

% All_TS_al
edges = -15*60e6:binsize_ms*1000:60*60e6;
[Q,bins_uS] = Bin_ts_array(All_TS_al, edges);
Q = Q(:,GIX); % filter out crappy neurons and neurons not in M1.
bins_uS = mean(bins_uS,2);

BASEIX = bins_uS < -15e6;
Qn = standardize_range(Q);
% Qn = standardize_range(Qn);
Qn = Qn - nanmean(Qn(BASEIX,:));
Qn = conv_filter(Qn,han_filt);
% Qn = standardize_range(Qn);
[~,six] = sort(mean(Qn,1));

M = Qn(:,six);


figure
x = bins_uS/60e6;
subplot(4,1,1:2)
imagesc(x,[],M')
ylabel('Neuron ID')
pubify_figure_axis
colorbar_label
set(gca,'XTickLabel',[])

plot_vert_line_at_zero()
subplot(4,1,3:4)
plot_confidence_intervals(x,M')
plot_vert_line_at_zero()
xlabel('Minutes')
pubify_figure_axis

%%
[SD,PCA]= Sliding_dynamics(Q, (1*60)/(binsize_ms/1000));
lbl = {'R_matrix_mn_r' 'S_matrix_skew_r' 'nEffDim' 'Rolls_Treves' 'CV' 'prop_active' 'Net_state_PC1' };

for ii = 1:length(lbl)
    figure
    plot(linspace(x(1),x(end),length(SD.(lbl{ii}))), nanmean( SD.(lbl{ii}),2))
    title(lbl{ii})
    pubify_figure_axis
    ylabel(lbl{ii})
    axis tight
    plot_vert_line_at_zero
end

%% Mean firing rate analyses...



binsize_sec = diff(Dset.Qbase_bins_uS(1:2))/1e6;
mn_b = mean(AllQ{1},1)/binsize_sec;
mn_p1 = mean(AllQ{2},1)/binsize_sec;
mn_p2 = mean(AllQ{3},1)/binsize_sec;

figure
subplot(2,2,1)
plot(mn_p1,mn_b,'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2, mn_b,'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch')

subplot(2,2,3:4)
p1 = signrank(mn_p1-mn_b);
mn1 = mean(mn_p1-mn_b);
p2 = signrank(mn_p2-mn_b);
mn2 = mean(mn_p2-mn_b);

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1-mn_b mn_p2-mn_b},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))
%% Popluation distance comparison (as in that Mohgaddahm paper we went over...

figure
subplot(1,2,1)
p1 = pdist([mn_b;mn_p1],'euclidean');
p2 = pdist([mn_b;mn_p2],'euclidean');
bar([p1 p2])
set(gca,'YLim',[12 19])
ylabel('Euclidean Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B'})
pubify_figure_axis

subplot(1,2,2)
p1 = pdist([mn_b;mn_p1],'spearman');
p2 = pdist([mn_b;mn_p2],'spearman');
bar([p1 p2])
ylabel('Spearman r Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B'})
pubify_figure_axis
