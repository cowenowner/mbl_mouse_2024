%% Local variance by rat
Rat_LID = [320 342 378 380 396];
Rat_SHAM = [350 352 376 392 393];
Rat_num = [320 342 378 380 396 350 352 376 392 393];
% creating a matrix with first colum with right hem data for all rats (first
% 5 LID and second 5 Sham) and second column is left hem data for all rats,
% third column has grouping for rat type
for ii = 1:length(Rat_num)
    LV_b_rat(ii,2) = nanmean(TBL.LocVar_base(GIX_Rat_right(ii,:)));
    LV_p1_rat(ii,2) = nanmean(TBL.LocVar_post1(GIX_Rat_right(ii,:)));
    LV_p2_rat(ii,2) = nanmean(TBL.LocVar_post2(GIX_Rat_right(ii,:)));
    LV_p3_rat(ii,2) = nanmean(TBL.LocVar_post3(GIX_Rat_right(ii,:)));
    
    LV_b_rat(ii,1) = nanmean(TBL.LocVar_base(GIX_Rat_left(ii,:)));
    LV_p1_rat(ii,1) = nanmean(TBL.LocVar_post1(GIX_Rat_left(ii,:)));
    LV_p2_rat(ii,1) = nanmean(TBL.LocVar_post2(GIX_Rat_left(ii,:)));
    LV_p3_rat(ii,1) = nanmean(TBL.LocVar_post3(GIX_Rat_left(ii,:)));
    
    if ii < 6
        LV_b_rat(ii,3) = 1;
        LV_p1_rat(ii,3) = 1;
        LV_p2_rat(ii,3) = 1;
        LV_p3_rat(ii,3) = 1;
    elseif ii > 5
        LV_b_rat(ii,3) = 2;
        LV_p1_rat(ii,3) = 2;
        LV_p2_rat(ii,3) = 2;
        LV_p3_rat(ii,3) = 2;
    end
end

numRows = size(LV_b_rat, 1);
selectedRows = 1:floor(numRows/2); 
selectedRows2 = ceil(numRows/2)+1:numRows;

figure
p1 = signrank(LV_p1_rat(selectedRows,1) -LV_b_rat(selectedRows,1));
mn1 = nanmean(LV_p1_rat(selectedRows,1) -LV_b_rat(selectedRows,1));
p2 = signrank(LV_p2_rat(selectedRows,1) -LV_b_rat(selectedRows,1));
mn2 = nanmean(LV_p2_rat(selectedRows,1) -LV_b_rat(selectedRows,1));
p2p1 = signrank(LV_p2_rat(selectedRows,1) -LV_p1_rat(selectedRows,1));
mn2m1 = nanmean(LV_p2_rat(selectedRows,1) -LV_p1_rat(selectedRows,1));

p11 = signrank(LV_p1_rat(selectedRows2,1) -LV_b_rat(selectedRows2,1));
mn11 = nanmean(LV_p1_rat(selectedRows2,1) -LV_b_rat(selectedRows2,1));
p22 = signrank(LV_p2_rat(selectedRows2,1) -LV_b_rat(selectedRows2,1));
mn22 = nanmean(LV_p2_rat(selectedRows2,1) -LV_b_rat(selectedRows2,1));
p2p11 = signrank(LV_p2_rat(selectedRows2,1) -LV_p1_rat(selectedRows2,1));
mn2m11 = nanmean(LV_p2_rat(selectedRows2,1) -LV_p1_rat(selectedRows2,1));

subplot(2,3,1)
gscatter(LV_p1_rat(:,1),LV_b_rat(:,1),LV_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
title(sprintf('Right m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1, mn11, p11))

subplot(2,3,2)
gscatter(LV_p2_rat(:,1),LV_b_rat(:,1),LV_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
title(sprintf('Right m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn2,p2, mn22, p22))

subplot(2,3,3)
gscatter(LV_p2_rat(:,1),LV_p1_rat(:,1),LV_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Post 1 (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
title(sprintf('Right m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn2m1,p2p1, mn2m11, p2p11))

p1 = signrank(LV_p1_rat(selectedRows,2) -LV_b_rat(selectedRows,2));
mn1 = nanmean(LV_p1_rat(selectedRows,2) -LV_b_rat(selectedRows,2));
p2 = signrank(LV_p2_rat(selectedRows,2) -LV_b_rat(selectedRows,2));
mn2 = nanmean(LV_p2_rat(selectedRows,2) -LV_b_rat(selectedRows,2));
p2p1 = signrank(LV_p2_rat(selectedRows,2) -LV_p1_rat(selectedRows,2));
mn2m1 = nanmean(LV_p2_rat(selectedRows,2) -LV_p1_rat(selectedRows,2));

p11 = signrank(LV_p1_rat(selectedRows2,2) -LV_b_rat(selectedRows2,2));
mn11 = nanmean(LV_p1_rat(selectedRows2,2) -LV_b_rat(selectedRows2,2));
p22 = signrank(LV_p2_rat(selectedRows2,2) -LV_b_rat(selectedRows2,2));
mn22 = nanmean(LV_p2_rat(selectedRows2,2) -LV_b_rat(selectedRows2,2));
p2p11 = signrank(LV_p2_rat(selectedRows2,2) -LV_p1_rat(selectedRows2,2));
mn2m11 = nanmean(LV_p2_rat(selectedRows2,2) -LV_p1_rat(selectedRows2,2));

subplot(2,3,4)
gscatter(LV_p1_rat(:,2),LV_b_rat(:,2),LV_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline(Hz)')
xlabel('Post 1 (Hz)')
pubify_figure_axis
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
title(sprintf('Left m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn1,p1, mn11, p11))

subplot(2,3,5)
gscatter(LV_p2_rat(:,2),LV_b_rat(:,2),LV_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
title(sprintf('Left m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn2,p2, mn22, p22))

subplot(2,3,6)
gscatter(LV_p2_rat(:,2),LV_p1_rat(:,2),LV_b_rat(:,3),'br', 'xo')
axis equal
axis square
ylabel('Post 1 (Hz)')
xlabel('Post 2 (Hz)')
set(gca,'XScale','log')
set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
title(sprintf('Left m1=%1.2f,p1=%1.3f, m2=%1.2f, p2=%1.3f',mn2m1,p2p1, mn2m11, p2p11))
sgtitle('Local variance by hemisphere by rat')






%%
for ii = 1:length(Rat_LID)
    GIX = sum(AllQ_all)' >500 & TBL.Depth_uM < 2000 & TBL.Rat == Rat_LID(ii);
    % GIXleft = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'L';
    % GIXright = sum(AllQ_all)'>500 & TBL.Depth_uM < 2000 & categorical(TBL.Hemisphere) == 'R';
    GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
    GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
    
%     figure
%     subplot(2,1,1)
%     histogram_cowen({TBL.LocVar_base(GIXright) TBL.LocVar_post1(GIXright) TBL.LocVar_post2(GIXright)},.1 )
%     title(sprintf('Rat %1.0f Right hemisphere',Rat_LID(ii)))
%     legend('base','post1','post2');legend boxoff
%     subplot(2,1,2)
%     histogram_cowen({TBL.LocVar_base(GIXleft) TBL.LocVar_post1(GIXleft) TBL.LocVar_post2(GIXleft)},.1 )
%     title('Left hemisphere')
%     legend('base','post1','post2');legend boxoff
%     
%     figure
%     subplot(2,1,1)
%     histogram_cowen({(TBL.LocVar_post1(GIXright) -TBL.LocVar_base(GIXright)) (TBL.LocVar_post2(GIXright) -TBL.LocVar_base(GIXright)) ...
%         (TBL.LocVar_post2(GIXright) -TBL.LocVar_post1(GIXright)) },.05 )
%     plot_vert_line_at_zero
%     legend('post1-base','post2-base', 'post2-post1');legend boxoff
%     plot_vert_line_at_zero
%     xlabel('Local Variance')
%     title(sprintf('Rat %1.0f Right hemisphere',Rat_LID(ii)))
    mp2_right(ii) = nanmean(TBL.LocVar_post2(GIXright) -TBL.LocVar_base(GIXright));
    mp1_right(ii) = nanmean(TBL.LocVar_post1(GIXright) -TBL.LocVar_base(GIXright));
    mp2p1_right(ii) = nanmean(TBL.LocVar_post2(GIXright) -TBL.LocVar_post1(GIXright));
    mp3_right(ii) = nanmean(TBL.LocVar_post3(GIXright) -TBL.LocVar_base(GIXright));
    
%     subplot(2,1,2)
%     histogram_cowen({(TBL.LocVar_post1(GIXleft) -TBL.LocVar_base(GIXleft)) (TBL.LocVar_post2(GIXleft) -TBL.LocVar_base(GIXleft)) ...
%         (TBL.LocVar_post2(GIXleft) -TBL.LocVar_post1(GIXleft))},.05 )
%     plot_vert_line_at_zero
%     legend('post1-base','post2-base', 'post2-post1');legend boxoff
%     plot_vert_line_at_zero
%     xlabel('Local Variance')
%     title('Left hemisphere')
    mp2_left(ii) = nanmean(TBL.LocVar_post2(GIXleft) -TBL.LocVar_base(GIXleft));
    mp1_left(ii) = nanmean(TBL.LocVar_post1(GIXleft) -TBL.LocVar_base(GIXleft));
    mp2p1_left(ii) = nanmean(TBL.LocVar_post2(GIXleft) -TBL.LocVar_post1(GIXleft));
    mp3_left(ii) = nanmean(TBL.LocVar_post3(GIXleft) -TBL.LocVar_base(GIXleft));

end

signrank(mp2_right)
signrank(mp1_right)
signrank(mp2p1_right)


signrank(mp2_left)
signrank(mp1_left)
signrank(mp2p1_left)
