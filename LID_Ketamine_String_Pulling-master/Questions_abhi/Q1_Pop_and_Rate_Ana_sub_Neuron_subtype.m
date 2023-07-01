%% population activity and firing rate analysis.
figure
% Limit to M1 and neurons that fired at least a little.
fsize = 40;
GIX = sum(AllQ_all)' > 400 & TBL.Depth_uM < 2000 & categorical(TBL.Group) == '6OHDA_LID';
% GIX = sum(AllQ_all)' >200 & TBL.Depth_uM > 3000 & categorical(TBL.Group) == 'SHAM';
% GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
% GIXrightPK = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';

GIXleftPK = GIX & categorical(TBL.Hemisphere) == 'L' & TBL.Neuron_type > 1;
GIXleftIN = GIX & categorical(TBL.Hemisphere) == 'L' & TBL.Neuron_type == 1;

GIXrightPK = GIX & categorical(TBL.Hemisphere) == 'R' & TBL.Neuron_type > 1;
GIXrightIN = GIX & categorical(TBL.Hemisphere) == 'R' & TBL.Neuron_type == 1;

sum(GIX), sum(GIXleftIN), sum(GIXrightIN)

x = Qall_x_uS(:,1)-Qall_x_uS(1,1);
x = x/60e6 - 25; %alingment to first inj; why is this hard coded??
%  M = log(double(AllQ_all'));
M = double(AllQ_all);
M = standardize_range(M);

M = M - nanmean(M(1:1200,:));

M = conv_filter(M,hanning(fsize)/sum(hanning(fsize)));
binsize_sec = 0.0200;
M = M/binsize_sec;
%% Neuron type 
Mleft_PK = M(:,GIXleftPK);
Mleft_IN = M(:,GIXleftIN);
[~,six_left_PK] = sort(mean(Mleft_PK,1));
[~,six_left_IN] = sort(mean(Mleft_IN,1));
Mleft_PK = Mleft_PK(:,six_left_PK);
Mleft_IN = Mleft_IN(:,six_left_IN);

Mright_PK = M(:,GIXrightPK);
Mright_IN = M(:,GIXrightIN);
[~,six_right_PK] = sort(mean(Mright_PK,1));
[~,six_right_IN] = sort(mean(Mright_IN,1));
Mright_PK = Mright_PK(:,six_right_PK);
Mright_IN = Mright_IN(:,six_right_IN);

Mleft_NT = [Mleft_IN Mleft_PK];

figure
subplot(4,1,1:2)
imagesc(x,[],Mleft_NT')
plot_horiz_line_at_zero(Cols(Mleft_IN),2,'w','-');
plot_horiz_line_at_zero(Cols(Mleft_IN),1,'k');
ylabel('Neuron')
pubify_figure_axis
colorbar_label
set(gca,'XTickLabel','')
plot_vert_line_at_zero()
plot_vert_line_at_zero(60)


subplot(4,1,3:4)
fs = 500;
Mleft_IN2 = conv_filter(Mleft_IN,hanning(fs)/sum(hanning(fs)));
Mleft_PK2 = conv_filter(Mleft_PK,hanning(fs)/sum(hanning(fs)));
plot_confidence_intervals(x,Mleft_IN2',[],GP.Colors.LeftPaw)
plot_confidence_intervals(x,Mleft_PK2',[],GP.Colors.RightPaw)
plot_vert_line_at_zero()
plot_vert_line_at_zero(60)
ylabel('Mean z-score firing rate')
xlabel('Time (min)')
pubify_figure_axis
legend_color_text({'Interneuron' 'Principal cells'},{GP.Colors.LeftPaw GP.Colors.RightPaw})

Mright_NT = [Mright_IN Mright_PK];

figure
subplot(4,1,1:2)
imagesc(x,[],Mright_NT')
plot_horiz_line_at_zero(Cols(Mright_IN),2,'w','-');
plot_horiz_line_at_zero(Cols(Mright_IN),1,'k');
ylabel('Neuron')
pubify_figure_axis
colorbar_label
set(gca,'XTickLabel','')
plot_vert_line_at_zero()
plot_vert_line_at_zero(60)


subplot(4,1,3:4)
fs = 500;
Mright_IN2 = conv_filter(Mright_IN,hanning(fs)/sum(hanning(fs)));
Mright_PK2 = conv_filter(Mright_PK,hanning(fs)/sum(hanning(fs)));
plot_confidence_intervals(x,Mright_IN2',[],GP.Colors.LeftPaw)
plot_confidence_intervals(x,Mright_PK2',[],GP.Colors.RightPaw)
plot_vert_line_at_zero()
plot_vert_line_at_zero(60)
ylabel('Mean z-score firing rate')
xlabel('Time (min)')
pubify_figure_axis
legend_color_text({'Interneurons' 'PK'},{GP.Colors.LeftPaw GP.Colors.RightPaw})

%% Mean firing rate analyses...
binsize_sec = diff(Qbase_bins_uS(1:2))/1e6;
mn_b = mean(AllQ_b,1)/binsize_sec;
mn_p1 = mean(AllQ_p1,1)/binsize_sec;
mn_p2 = mean(AllQ_p2,1)/binsize_sec;
mn_p3 = mean(AllQ_p3,1)/binsize_sec;
% Right PK
figure
subplot(2,2,1)
plot(mn_p1(GIXrightPK),mn_b(GIXrightPK),'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2(GIXrightPK), mn_b(GIXrightPK),'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch right hemisphere in PK')

subplot(2,2,3:4)
p1 = signrank(mn_p1(GIXrightPK)-mn_b(GIXrightPK));
mn1 = mean(mn_p1(GIXrightPK)-mn_b(GIXrightPK));
p2 = signrank(mn_p2(GIXrightPK)-mn_b(GIXrightPK));
mn2 = mean(mn_p2(GIXrightPK)-mn_b(GIXrightPK));
p3 = signrank(mn_p3(GIXrightPK)-mn_b(GIXrightPK));
mn3 = mean(mn_p3(GIXrightPK)-mn_b(GIXrightPK));
p2mp1 = signrank(mn_p2(GIXrightPK)-mn_p1(GIXrightPK));
mn2m1 = mean(mn_p2(GIXrightPK)-mn_p1(GIXrightPK));

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1(GIXrightPK)-mn_b(GIXrightPK) mn_p2(GIXrightPK)-mn_b(GIXrightPK)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))
%IN
figure
subplot(2,2,1)
plot(mn_p1(GIXrightIN),mn_b(GIXrightIN),'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2(GIXrightIN), mn_b(GIXrightIN),'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch right hemisphere in Interneurons')

subplot(2,2,3:4)
p1 = signrank(mn_p1(GIXrightIN)-mn_b(GIXrightIN));
mn1 = mean(mn_p1(GIXrightIN)-mn_b(GIXrightIN));
p2 = signrank(mn_p2(GIXrightIN)-mn_b(GIXrightIN));
mn2 = mean(mn_p2(GIXrightIN)-mn_b(GIXrightIN));
p3 = signrank(mn_p3(GIXrightIN)-mn_b(GIXrightIN));
mn3 = mean(mn_p3(GIXrightIN)-mn_b(GIXrightIN));
p2mp1 = signrank(mn_p2(GIXrightIN)-mn_p1(GIXrightIN));
mn2m1 = mean(mn_p2(GIXrightIN)-mn_p1(GIXrightIN));

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1(GIXrightIN)-mn_b(GIXrightIN) mn_p2(GIXrightIN)-mn_b(GIXrightIN)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))

% LEft PK
figure
subplot(2,2,1)
plot(mn_p1(GIXleftPK),mn_b(GIXleftPK),'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2(GIXleftPK), mn_b(GIXleftPK),'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch left hemisphere in MSNs')

subplot(2,2,3:4)
p1 = signrank(mn_p1(GIXleftPK)-mn_b(GIXleftPK));
mn1 = mean(mn_p1(GIXleftPK)-mn_b(GIXleftPK));
p2 = signrank(mn_p2(GIXleftPK)-mn_b(GIXleftPK));
mn2 = mean(mn_p2(GIXleftPK)-mn_b(GIXleftPK));
p3 = signrank(mn_p3(GIXleftPK)-mn_b(GIXleftPK));
mn3 = mean(mn_p3(GIXleftPK)-mn_b(GIXleftPK));
p2mp1 = signrank(mn_p2(GIXleftPK)-mn_p1(GIXleftPK));
mn2m1 = mean(mn_p2(GIXleftPK)-mn_p1(GIXleftPK));

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1(GIXleftPK)-mn_b(GIXleftPK) mn_p2(GIXleftPK)-mn_b(GIXleftPK)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))

% Left In
figure
subplot(2,2,1)
plot(mn_p1(GIXleftIN),mn_b(GIXleftIN),'bo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(2,2,2)
plot(mn_p2(GIXleftIN), mn_b(GIXleftIN),'bo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch left hemisphere in Interneurons')

subplot(2,2,3:4)
p1 = signrank(mn_p1(GIXleftIN)-mn_b(GIXleftIN));
mn1 = mean(mn_p1(GIXleftIN)-mn_b(GIXleftIN));
p2 = signrank(mn_p2(GIXleftIN)-mn_b(GIXleftIN));
mn2 = mean(mn_p2(GIXleftIN)-mn_b(GIXleftIN));
p3 = signrank(mn_p3(GIXleftIN)-mn_b(GIXleftIN));
mn3 = mean(mn_p3(GIXleftIN)-mn_b(GIXleftIN));
p2mp1 = signrank(mn_p2(GIXleftIN)-mn_p1(GIXleftIN));
mn2m1 = mean(mn_p2(GIXleftIN)-mn_p1(GIXleftIN));

edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1(GIXleftIN)-mn_b(GIXleftIN) mn_p2(GIXleftIN)-mn_b(GIXleftIN)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz)')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))


%% Difference between Pk and In
p1 = ranksum(mn_p1(GIXrightPK), mn_p1(GIXrightIN));
p2 = ranksum(mn_p2(GIXrightPK), mn_p2(GIXrightIN));
p2_d = Cohens_d(mn_p2(GIXrightPK), mn_p2(GIXrightIN));
p11 = ranksum(mn_p1(GIXleftPK), mn_p1(GIXleftIN));
p22 = ranksum(mn_p2(GIXleftPK), mn_p2(GIXleftIN));
p22_d = Cohens_d(mn_p2(GIXleftPK), mn_p2(GIXleftIN));




