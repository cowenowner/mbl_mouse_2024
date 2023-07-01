%% population activity and firing rate analysis.
figure
% Limit to M1 and neurons that fired at least a little.
fsize = 40;
binsize_sec = bin_sec;
% GIX = sum(AllQ_all)' > 400 & TBL.Depth_uM < 2000 & categorical(TBL.Group) == '6OHDA_LID';
GIX = sum(AllQ_all)' >200 & TBL.Depth_uM > 3000 & categorical(TBL.Group) == '6OHDA_LID';
% GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
% GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
sum(GIX), sum(GIXleft), sum(GIXright)

GIXleftPK = GIXleft & TBL.Neuron_type > 1;
GIXleftIN = GIXleft & TBL.Neuron_type == 1;

GIXrightPK = GIXright & TBL.Neuron_type > 1;
GIXrightIN = GIXright & TBL.Neuron_type == 1;

x = Qall_x_uS(:,1)-Qall_x_uS(1,1);
x = x/60e6 - 25; %alingment to first inj; why is this hard coded??
%  M = log(double(AllQ_all'));
M = double(AllQ_all);
M = standardize_range(M);

M = M - nanmean(M(1:1200,:));

M = conv_filter(M,hanning(fsize)/sum(hanning(fsize)));

M = M/binsize_sec;

Mall = M(:,GIX);
Mleft = M(:,GIXleft);
Mright = M(:,GIXright);

[~,six_all] = sort(mean(Mall,1));
[~,six_left] = sort(mean(Mleft,1));
[~,six_right] = sort(mean(Mright,1));

Mall = Mall(:,six_all);
Mleft = Mleft(:,six_left);
Mright = Mright(:,six_right);

MAll = [Mleft Mright];
figure
subplot(4,1,1:2)
imagesc(x,[],MAll')
plot_horiz_line_at_zero(Cols(Mleft),2,'w','-');
plot_horiz_line_at_zero(Cols(Mleft),1,'k');
ylabel('Neuron')
pubify_figure_axis
colorbar_label
set(gca,'XTickLabel','')
plot_vert_line_at_zero()
plot_vert_line_at_zero(60)


subplot(4,1,3:4)
fs = 500;
Mleft2 = conv_filter(Mleft,hanning(fs)/sum(hanning(fs)));
Mright2 = conv_filter(Mright,hanning(fs)/sum(hanning(fs)));
plot_confidence_intervals(x,Mleft2',[],GP.Colors.LeftPaw)
plot_confidence_intervals(x,Mright2',[],GP.Colors.RightPaw)
plot_vert_line_at_zero()
plot_vert_line_at_zero(60)
ylabel('Mean firing rate')
xlabel('Time (min)')
pubify_figure_axis
% legend_color_text({'Un-Lesioned' 'Lesioned'},{GP.Colors.LeftPaw GP.Colors.RightPaw})
legend_color_text({'Intact' 'Sham-Lesioned'},{GP.Colors.LeftPaw GP.Colors.RightPaw})




%% Just count numbers of cells that increase or decrease.
d = mean(AllQ{3}) - mean(AllQ{2});
dleft = d(GIXleft);
dright = d(GIXright);
updownneith = [sum(dleft>0) sum(dleft<0); sum(dright>0) sum(dright<0)];
figure
bar(updownneith)
set(gca,'XTickLabel',{'Left Hem', 'Right Hem'})
ylabel('# neurons')
pubify_figure_axis
legend('Increase','Decrease');legend boxoff

%%

figure
[~,six] = sort(mean(AllQ{1},1));
lbls = {'Baseline' 'LID' 'Post Ketamine'};
han_size = 50;
rr = [];
for ii = 1:3
    AllQ{ii} = double(AllQ{ii});
    MAll = AllQ{ii}(:,six)/(Dset.binsize_Q_ms/1000);
    MAll = conv_filter(MAll,hanning(han_size)/sum(hanning(han_size)));
    x_min = (1:Rows(AllQ{ii}))*Dset.binsize_Q_ms/60e3;
    
    figure(111)
    subplot(1,3,ii)
    imagesc(x_min,[],MAll')
    xlabel('Min')
    if ii == 1
        ylabel('Neuron')
    end
    title(lbls{ii})

    pubify_figure_axis
    
    figure(112)
    subplot(1,3,ii)
    plot_confidence_intervals(x_min,MAll')
    pubify_figure_axis
    xlabel('min')
    ylabel('mean rate')
    
    
    
end
equalize_axes
figure(111)
equalize_color_axes
colorbar_label

%%
figure
for ii = 1:2
    MAll = (double(AllQ{ii+1})-mean(double(AllQ{1})))'/bin_sec;
    MAll = conv_filter(MAll',hanning(han_size)/sum(hanning(han_size)))';
    x_min = (1:Rows(AllQ{+1}))*Dset.binsize_Q_ms/60e3;
    
    
    [~,six] = sort(mean(MAll,2));
    subplot(1,2,ii)
    imagesc(x_min,[], MAll(six,:))
    xlabel('Min')
    ylabel('Neuron')
    title(sprintf('Diff from base Post %d',ii))
    pubify_figure_axis
    caxis([-2 2])
    colormap(hotcold_colormap)
end
equalize_color_axes
colorbar_label
%% Mean firing rate analyses...
binsize_sec = diff(Qbase_bins_uS(1:2))/1e6;
mn_b = mean(AllQ_b,1)/binsize_sec;
mn_p1 = mean(AllQ_p1,1)/binsize_sec;
mn_p2 = mean(AllQ_p2,1)/binsize_sec;
mn_p3 = mean(AllQ_p3,1)/binsize_sec;

figure
subplot(2,2,1)
plot(mn_p1(GIXright),mn_b(GIXright),'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2(GIXright), mn_b(GIXright),'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch right hemisphere')

subplot(2,2,3:4)
p1 = signrank(mn_p1(GIXright)-mn_b(GIXright));
mn1 = mean(mn_p1(GIXright)-mn_b(GIXright));
p2 = signrank(mn_p2(GIXright)-mn_b(GIXright));
mn2 = mean(mn_p2(GIXright)-mn_b(GIXright));
p3 = signrank(mn_p3(GIXright)-mn_b(GIXright));
mn3 = mean(mn_p3(GIXright)-mn_b(GIXright));
p2mp1 = signrank(mn_p2(GIXright)-mn_p1(GIXright));
mn2m1 = mean(mn_p2(GIXright)-mn_p1(GIXright));

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1(GIXright)-mn_b(GIXright) mn_p2(GIXright)-mn_b(GIXright)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))


figure
subplot(2,2,1)
plot(mn_p1(GIXleft),mn_b(GIXleft),'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2(GIXleft), mn_b(GIXleft),'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch left hemisphere')

subplot(2,2,3:4)
p1 = signrank(mn_p1(GIXleft)-mn_b(GIXleft));
mn1 = mean(mn_p1(GIXleft)-mn_b(GIXleft));
p2 = signrank(mn_p2(GIXleft)-mn_b(GIXleft));
mn2 = mean(mn_p2(GIXleft)-mn_b(GIXleft));
p3 = signrank(mn_p3(GIXleft)-mn_b(GIXleft));
mn3 = mean(mn_p3(GIXleft)-mn_b(GIXleft));
p2mp1 = signrank(mn_p2(GIXleft)-mn_p1(GIXleft));
mn2m1 = mean(mn_p2(GIXleft)-mn_p1(GIXleft));

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1(GIXleft)-mn_b(GIXleft) mn_p2(GIXleft)-mn_b(GIXleft)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))


%% comparing left sham hem and un-lesioned hem during LDOPA
GIXShamL = sum(AllQ_all)' > 400 & TBL.Depth_uM < 2000 & categorical(TBL.Group) == 'SHAM' & categorical(TBL.Hemisphere) == 'L';
GIXLIDL = sum(AllQ_all)' > 400 & TBL.Depth_uM < 2000 & categorical(TBL.Group) == '6OHDA_LID' & categorical(TBL.Hemisphere) == 'L';

p1 = ranksum(mn_p1(GIXLIDL),mn_p1(GIXShamL));
p1_d = Cohens_d(mn_p1(GIXLIDL),mn_p1(GIXShamL));

%% Popluation distance comparison (as in that Mohgaddahm paper we went over...

figure
subplot(2,2,1)
p1b = pdist([mn_b(GIXright);mn_p1(GIXright)],'euclidean');
p2b = pdist([mn_b(GIXright);mn_p2(GIXright)],'euclidean');
p3b = pdist([mn_b(GIXright);mn_p3(GIXright)],'euclidean');
p1p2 = pdist([mn_p1(GIXright);mn_p2(GIXright)],'euclidean');
bar([p1b p2b p3b p1p2])
% set(gca,'YLim',[12 19])
ylabel('Euclidean Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B', 'P3-B', 'P2-P1'})
title('Right hemisphere')
pubify_figure_axis

subplot(2,2,2)
p1b = pdist([mn_b(GIXleft);mn_p1(GIXleft)],'euclidean');
p2b = pdist([mn_b(GIXleft);mn_p2(GIXleft)],'euclidean');
p3b = pdist([mn_b(GIXleft);mn_p3(GIXleft)],'euclidean');
p1p2 = pdist([mn_p1(GIXleft);mn_p2(GIXleft)],'euclidean');
bar([p1b p2b p3b p1p2])
% set(gca,'YLim',[12 19])
ylabel('Euclidean Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B', 'P3-B', 'P2-P1'})
title('Left hemisphere')
pubify_figure_axis

subplot(2,2,3)
p1b = pdist([mn_b(GIXright);mn_p1(GIXright)],'spearman');
p2b = pdist([mn_b(GIXright);mn_p2(GIXright)],'spearman');
p3b = pdist([mn_b(GIXright);mn_p3(GIXright)],'spearman');
p1p2 = pdist([mn_p1(GIXright);mn_p2(GIXright)],'spearman');
bar([p1b p2b p3b p1p2])
ylabel('Spearman r Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B', 'P3-B', 'P2-P1'})
pubify_figure_axis


subplot(2,2,4)
p1b = pdist([mn_b(GIXleft);mn_p1(GIXleft)],'spearman');
p2b = pdist([mn_b(GIXleft);mn_p2(GIXleft)],'spearman');
p3b = pdist([mn_b(GIXleft);mn_p3(GIXleft)],'spearman');
p1p2 = pdist([mn_p1(GIXleft);mn_p2(GIXleft)],'spearman');
bar([p1b p2b p3b p1p2])
ylabel('Spearman r Distance from Baseline')
set(gca,'XTickLabel',{'P1-B','P2-B', 'P3-B', 'P2-P1'})
pubify_figure_axis

