%% Firing rate analysis averaged by rat
% Limit to M1 and neurons that fired at least a little.
fsize = 40;
Rat_LID = [320 342 378 380 396];
Rat_SHAM = [350 352 376 392 393];
Rat_num = [320 342 378 380 396 350 352 376 392 393];

GIX = sum(AllQ_all)' >400 & TBL.Depth_uM < 2000;
% % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
% % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
 GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
 GIXright = GIX & categorical(TBL.Hemisphere) == 'R';

sum(GIX), sum(GIXleft), sum(GIXright)

for ii = 1:length(Rat_num)
   GIX_Rat_left(ii,:) = GIXleft & TBL.Rat == Rat_num(ii); 
   GIX_Rat_right(ii,:) = GIXright & TBL.Rat == Rat_num(ii);     
end

x = Dset.Qall_x_uS(:,1)-Dset.Qall_x_uS(1,1);
x = x/60e6 - 25; %alingment to first inj; why is this hard coded??
%  M = log(double(AllQ_all'));
M = double(AllQ_all);
M = standardize_range(M);

M = M - nanmean(M(1:1200,:)); % baseline subtraction 1200 = 20 minutes at 500 ms binned data

M = conv_filter(M,hanning(fsize)/sum(hanning(fsize)));

for ii = 1:length(Rat_LID)
    GIX = sum(AllQ_all)' >500 & TBL.Depth_uM < 2000 & TBL.Rat == Rat_LID(4);
    XIX = mean(M,1) > -0.1;
    GIX = GIX & XIX';
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
    Mall = M(:,GIX);
    Mleft = M(:,GIXleft);
    Mright = M(:,GIXright);
    
    [~,six_all] = sort(mn_p1(GIX));
    [~,six_left] = sort(mn_p1(GIXleft));
    [~,six_right] = sort(mn_p2(GIXright)-mn_p1(GIXright));
%     [~,six_all] = sort(mean(Mall,1));
%     [~,six_left] = sort(mean(Mleft,1));
%     [~,six_right] = sort(mean(Mright,1));
    
    Mall = Mall(:,six_all);
    Mleft = Mleft(:,six_left);
    Mright = Mright(:,six_right);
    
    M_LR = [Mleft Mright];
    
    
    figure
    subplot(2,1,1:2)
    imagesc(x,[],Mright(IX)')
    plot_horiz_line_at_zero(Cols(Mleft),2,'w','-');
    plot_horiz_line_at_zero(Cols(Mleft),1,'k');
    ylabel('Neuron')
    pubify_figure_axis
    colorbar_label
    set(gca,'XTickLabel','')
    plot_vert_line_at_zero()
    title(sprintf('LID Rat number %1.0f',Rat_LID(4)))
    
    
    subplot(4,1,3:4)
    fs = 500;
    Mleft2 = conv_filter(Mleft,hanning(fs)/sum(hanning(fs)));
    Mright2 = conv_filter(Mright,hanning(fs)/sum(hanning(fs)));
    plot_confidence_intervals(x,Mleft2',[],GP.Colors.LeftPaw)
    plot_confidence_intervals(x,Mright2',[],GP.Colors.RightPaw)
    plot_vert_line_at_zero()
    xlabel('Time (min)')
    pubify_figure_axis
    legend_color_text({'L' 'R'},{GP.Colors.LeftPaw GP.Colors.RightPaw})
    %
    
    All_M_left{ii} = Mleft2;
    All_M_right{ii} = Mright2;

end

clrs = lines(5);
figure
for iu = 1:length(Rat_LID)
    plot_confidence_intervals(x,All_M_left{iu}',[],clrs(iu,:))
    
end
xlabel('Minutes')
pubify_figure_axis


%% Statistics on mean Firing rate difference
binsize_sec = diff(Dset.Qbase_bins_uS(1:2))/1e6;
mn_b = mean(AllQ_b,1)/binsize_sec;
mn_p1 = mean(AllQ_p1,1)/binsize_sec;
mn_p2 = mean(AllQ_p2,1)/binsize_sec;
mn_p3 = mean(AllQ_p3,1)/binsize_sec;

for ii = 1:length(Rat_num)
    mn_b_rat(ii,2) = mean(mn_b(GIX_Rat_left(ii,:)));
    mn_p1_rat(ii,2) = mean(mn_p1(GIX_Rat_left(ii,:)));
    mn_p2_rat(ii,2) = mean(mn_p2(GIX_Rat_left(ii,:)));
    mn_p3_rat(ii,2) = mean(mn_p3(GIX_Rat_left(ii,:)));
    
    mn_b_rat(ii,1) = mean(mn_b(GIX_Rat_right(ii,:)));
    mn_p1_rat(ii,1) = mean(mn_p1(GIX_Rat_right(ii,:)));
    mn_p2_rat(ii,1) = mean(mn_p2(GIX_Rat_right(ii,:)));
    mn_p3_rat(ii,1) = mean(mn_p3(GIX_Rat_right(ii,:)));
    
    if ii < 6
        mn_b_rat(ii,3) = 1;
        mn_p1_rat(ii,3) = 1;
        mn_p2_rat(ii,3) = 1;
        mn_p3_rat(ii,3) = 1;
    elseif ii > 5
        mn_b_rat(ii,3) = 2;
        mn_p1_rat(ii,3) = 2;
        mn_p2_rat(ii,3) = 2;
        mn_p3_rat(ii,3) = 2;
    end
end

figure
subplot(3,2,1)
gscatter(mn_p1_rat(:,1),mn_b_rat(:,1),mn_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(3,2,2)
gscatter(mn_p2_rat(:,1),mn_b_rat(:,1),mn_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch right hemisphere by rat')

subplot(3,2,3:4)
numRows = size(mn_b_rat, 1);
selectedRows = 1:floor(numRows/2); 

p1 = signrank(mn_p1_rat(selectedRows,1)-mn_b_rat(selectedRows,1));
mn1 = mean(mn_p1_rat(selectedRows,1)-mn_b_rat(selectedRows,1));
p2 = signrank(mn_p2_rat(selectedRows,1)-mn_b_rat(selectedRows,1));
mn2 = mean(mn_p2_rat(selectedRows,1)-mn_b_rat(selectedRows,1));


edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1_rat(selectedRows,1)-mn_b_rat(selectedRows,1) ...
    mn_p2_rat(selectedRows,1)-mn_b_rat(selectedRows,1)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz) LID')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))

subplot(3,2,5:6)
selectedRows2 = ceil(numRows/2)+1:numRows;

p1 = signrank(mn_p1_rat(selectedRows2,1)-mn_b_rat(selectedRows2,1));
mn1 = mean(mn_p1_rat(selectedRows2,1)-mn_b_rat(selectedRows2,1));
p2 = signrank(mn_p2_rat(selectedRows2,1)-mn_b_rat(selectedRows2,1));
mn2 = mean(mn_p2_rat(selectedRows2,1)-mn_b_rat(selectedRows2,1));


edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1_rat(selectedRows2,1)-mn_b_rat(selectedRows2,1) ...
    mn_p2_rat(selectedRows2,1)-mn_b_rat(selectedRows2,1)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz) Sham')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))

% Left hem
figure
subplot(3,2,1)
gscatter(mn_p1_rat(:,2),mn_b_rat(:,2),mn_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line

subplot(3,2,2)
gscatter(mn_p2_rat(:,2),mn_b_rat(:,2),mn_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
sgtitle('Change in mean rates by epoch left hemisphere by rat')

subplot(3,2,3:4)
numRows = size(mn_b_rat, 1);
selectedRows = 1:floor(numRows/2); 

p1 = signrank(mn_p1_rat(selectedRows,2)-mn_b_rat(selectedRows,2));
mn1 = mean(mn_p1_rat(selectedRows,2)-mn_b_rat(selectedRows,2));
p2 = signrank(mn_p2_rat(selectedRows,1)-mn_b_rat(selectedRows,2));
mn2 = mean(mn_p2_rat(selectedRows,2)-mn_b_rat(selectedRows,2));


edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1_rat(selectedRows,2)-mn_b_rat(selectedRows,2) ...
    mn_p2_rat(selectedRows,2)-mn_b_rat(selectedRows,2)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz) LID')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))

subplot(3,2,5:6)
selectedRows2 = ceil(numRows/2)+1:numRows;

p1 = signrank(mn_p1_rat(selectedRows2,2)-mn_b_rat(selectedRows2,2));
mn1 = mean(mn_p1_rat(selectedRows2,2)-mn_b_rat(selectedRows2,2));
p2 = signrank(mn_p2_rat(selectedRows2,2)-mn_b_rat(selectedRows2,2));
mn2 = mean(mn_p2_rat(selectedRows2,2)-mn_b_rat(selectedRows2,2));


edges = (-3.1:.2:3.1); % need to add the .1 to make sure the middle bin straddles zero.
histogram_cowen({mn_p1_rat(selectedRows2,2)-mn_b_rat(selectedRows2,2) ...
    mn_p2_rat(selectedRows2,2)-mn_b_rat(selectedRows2,2)},edges)
legend('P1-B','P2-B')
legend boxoff
xlabel('Change in firing rate (Hz) Sham')
plot_vert_line_at_zero
set(gca,'XLim',[-3 3])
title(sprintf('m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1,mn2,p2))



%% look at cells that fire at least once every 100 ms
XIX = min([mn_b; mn_p1; mn_p2])>0.1;

for ii = 1:length(Rat_LID)
    GIX = XIX' & TBL.Depth_uM < 2000 & TBL.Rat == Rat_LID(ii);
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
    mn_p2p1_right(ii) = mean(mn_p2(GIXright)-mn_p1(GIXright));
    mn_p2p1_left(ii) = mean(mn_p2(GIXleft)-mn_p1(GIXleft));
    
    mn_p1_right(ii) = mean(mn_p1(GIXright)-mn_b(GIXright));
    mn_p1_left(ii) = mean(mn_p1(GIXleft)-mn_b(GIXleft));
    
    mn_p2_right(ii) = mean(mn_p2(GIXright)-mn_b(GIXright));
    mn_p2_left(ii) = mean(mn_p2(GIXleft)-mn_b(GIXleft));
end

p2p1_r = signrank(mn_p2p1_right);
p2p1_l = signrank(mn_p2p1_left);

p1_r = signrank(mn_p1_right);
p1_l = signrank(mn_p1_left);

p2_r = signrank(mn_p2_right);
p2_l = signrank(mn_p2_left);

%% LID rats looking at the pre and post diff histograms
% Checking to see if we get a bimodal distribution and seeing these changes
% are significant across rats?

for ii = 1:length(Rat_LID)
    GIX = TBL.Depth_uM < 2000 & TBL.Rat == Rat_LID(ii);
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
    p2p1_right{ii} = mn_p2(GIXright)-mn_p1(GIXright);
    p2p1_left{ii} = mn_p2(GIXleft)-mn_p1(GIXleft);
    
    p1_right{ii} = mn_p1(GIXright)-mn_b(GIXright);
    p1_left{ii} = mn_p1(GIXleft)-mn_b(GIXleft);
    
    p2_right{ii} = mn_p2(GIXright)-mn_b(GIXright);
    p2_left{ii} = mn_p2(GIXleft)-mn_b(GIXleft);
end

figure
histogram_cowen({p2p1_right{1} p2p1_right{2} p2p1_right{4} p2p1_right{4}}, .5) 

    figure
for ii = 1:length(Rat_LID)
    GIX = XIX' & TBL.Depth_uM < 2000 & TBL.Rat == Rat_LID(ii);
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
    subplot(4,1,ii)
    histogram_cowen({mn_b(GIXright) mn_p1(GIXright) mn_p2(GIXright)},.5)
    title(sprintf('Rat %1.0f', Rat_LID(ii)))
    
end

    figure
for ii = 1:length(Rat_LID)
    GIX = XIX' & TBL.Depth_uM < 2000 & TBL.Rat == Rat_LID(ii);
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
%     [~,six_b_r] = sort(mn_b(GIXright));
    
end

%% Sham rats
for ii = 1:length(Rat_SHAM)
    GIX = sum(AllQ_all)' >500 & TBL.Depth_uM < 2000 & TBL.Rat == Rat_SHAM(ii);
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
    Mall = M(:,GIX);
    Mleft = M(:,GIXleft);
    Mright = M(:,GIXright);
    
    [~,six_all] = sort(mean(Mall,1));
    [~,six_left] = sort(mean(Mleft,1));
    [~,six_right] = sort(mean(Mright,1));
    
    Mall = Mall(:,six_all);
    Mleft = Mleft(:,six_left);
    Mright = Mright(:,six_right);
    
    M_LR = [Mleft Mright];
    
    
    figure
    subplot(4,1,1:2)
    imagesc(x,[],M_LR')
    plot_horiz_line_at_zero(Cols(Mleft),2,'w','-');
    plot_horiz_line_at_zero(Cols(Mleft),1,'k');
    ylabel('Neuron ID')
    pubify_figure_axis
    colorbar_label
    set(gca,'XTickLabel','')
    plot_vert_line_at_zero()
    title(sprintf('Sham Rat number %1.0f',Rat_SHAM(ii)))
    
    
    subplot(4,1,3:4)
    fs = 500;
    Mleft2 = conv_filter(Mleft,hanning(fs)/sum(hanning(fs)));
    Mright2 = conv_filter(Mright,hanning(fs)/sum(hanning(fs)));
    plot_confidence_intervals(x,Mleft2',[],GP.Colors.LeftPaw)
    plot_confidence_intervals(x,Mright2',[],GP.Colors.RightPaw)
    plot_vert_line_at_zero()
    xlabel('Minutes')
    pubify_figure_axis
    legend_color_text({'L' 'R'},{GP.Colors.LeftPaw GP.Colors.RightPaw})
    %
    
    All_M_left{ii} = Mleft2;
    All_M_right{ii} = Mright2;

end

clrs = lines(5);
figure
for iu = 1:length(Rat_SHAM)
    plot_confidence_intervals(x,All_M_left{iu}',[],clrs(iu,:))
    
end
xlabel('Minutes')
pubify_figure_axis


%% Statistics on mean Firing rate difference
binsize_sec = diff(Dset.Qbase_bins_uS(1:2))/1e6;
mn_b = mean(AllQ{1},1)/binsize_sec;
mn_p1 = mean(AllQ{2},1)/binsize_sec;
mn_p2 = mean(AllQ{3},1)/binsize_sec;
mn_p3 = mean(AllQ{4},1)/binsize_sec;
% look at cells that fire at least once every 100 ms
XIX = min([mn_b; mn_p1; mn_p2])>0.1;

mn_p2p1_right = [];
mn_p2p1_left = [];

mn_p1_right = [];
mn_p1_left =[];

mn_p2_right = [];
mn_p2_left = [];

for ii = 1:length(Rat_SHAM)
    
    GIX = sum(AllQ_all)' >500 & TBL.Depth_uM < 2000 & TBL.Rat == Rat_SHAM(ii);
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
    mn_p2p1_right(ii) = mean(mn_p2(GIXright)-mn_p1(GIXright));
    mn_p2p1_left(ii) = mean(mn_p2(GIXleft)-mn_p1(GIXleft));
    
    mn_p1_right(ii) = mean(mn_p1(GIXright)-mn_b(GIXright));
    mn_p1_left(ii) = mean(mn_p1(GIXleft)-mn_b(GIXleft));
    
    mn_p2_right(ii) = mean(mn_p2(GIXright)-mn_b(GIXright));
    mn_p2_left(ii) = mean(mn_p2(GIXleft)-mn_b(GIXleft));
end

p2p1_r = signrank(mn_p2p1_right);
p2p1_l = signrank(mn_p2p1_left);

p1_r = signrank(mn_p1_right);
p1_l = signrank(mn_p1_left);

p2_r = signrank(mn_p2_right);
p2_l = signrank(mn_p2_left);


