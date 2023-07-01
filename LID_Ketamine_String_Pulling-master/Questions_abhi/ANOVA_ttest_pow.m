%% Run the Q10 to load the conditions
% across conditions post ketamine [5 25]
ket_eff_win = [2 30];
IX = O.t_sec > ket_eff_win(1)*60 & O.t_sec < ket_eff_win(2)*60;
peak_pow = [];
% peak_pow = zeros(length(IX(IX==1)));
for ii = 1:length(conditions_to_compare)
%     M = squeeze(ALL.BYFRQ{ii}(6,:,:) - (ALL.BYFRQ{ii}(5,:,:)+ALL.BYFRQ{ii}(7,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(3,:,:) - (ALL.BYFRQ{ii}(2,:,:)+ALL.BYFRQ{ii}(5,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(7,:,:) - (ALL.BYFRQ{ii}(6,:,:)+ALL.BYFRQ{ii}(8,:,:))/2)';
%     M = squeeze(ALL.BYFRQ{ii}(4,:,:)- (ALL.BYFRQ{ii}(1,:,:)+ALL.BYFRQ{ii}(8,:,:))/2)';
    M = squeeze(ALL.BYFRQ{ii}(4,:,:))';
    Ms = movmedian(M,10,2);
    Ms_eff_win = Ms(:,IX)';
    peak_pow{ii} = mean(Ms_eff_win)'; 
end


p2(4) = ranksum(peak_pow{1},peak_pow{2});
d2(4) = Cohens_d(peak_pow{2},peak_pow{1});

[deci,pval_les_unles] = ttest(peak_pow{1}',peak_pow{2}');

peak_pow_mat = cell2mat(peak_pow)';
g1 = ones(1,111);
g2 = 2*ones(1,111);
g3 = 3*ones(1,111);
group = [g1,g2,g3];
[p,tbl,stats] = anovan(peak_pow_mat',group');
[c,m,h,nms] = multcompare(stats);
figure; boxplot(peak_pow_mat',group');

[h,p_CM1_Les] = ttest(peak_pow_mat(:,3),peak_pow_mat(:,2));
[h,p_CM1_UnLes] = ttest(peak_pow_mat(:,3),peak_pow_mat(:,1));
[h,p_Unles_Les] = ttest(peak_pow_mat(:,1),peak_pow_mat(:,2));

[p_C_UnL,h,stats3_1] = ranksum(peak_pow{1},peak_pow{3});
[p_C_L,h,stats3_2] = ranksum(peak_pow{2},peak_pow{3});
[p_L_UnL,h,stats1_2] = ranksum(peak_pow{1},peak_pow{2});

[h,p_C,ci,stats3] = ttest(peak_pow{3});
[h,p_UL,ci,stats2] = ttest(peak_pow{1});
[h,p_L,ci,stats] = ttest(peak_pow{2});

% [p_L_S_K(2),h,stats1_2] = ranksum(peak_pow{1},peak_pow{2});
[p_UL_S_K(3),h,stats1_2] = ranksum(peak_pow{1},peak_pow{2});
d_UL_S_k(2) = Cohens_d(peak_pow{2},peak_pow{1});
% M1 and striatum in Saline and ketamine
% [h,pt(3)] = ttest2(peak_pow{2},peak_pow{3},'Vartype','unequal');
p(1) = ranksum(peak_pow{1},peak_pow{2});
p(2) = ranksum(peak_pow{1},peak_pow{3});
p(3) = ranksum(peak_pow{3},peak_pow{2});
d(1) = Cohens_d(peak_pow{1},peak_pow{2});
d(2) = Cohens_d(peak_pow{1},peak_pow{3});
d(3) = Cohens_d(peak_pow{2},peak_pow{3});


test = ranksum(peak_pow{1},peak_pow{2});
%%%%%%% Striatum - 140 Hz
%  Unlesione and Lesion [2.73374026369400e-05]
%  Unlesion and Control [5.12227276420196e-05]
%  Lesion and Control [2.82307022203554e-09]
%%%%%%% M1 - 40-75 hz
%  Unlesione and Lesion [1.09488540456226e-09]
%  Unlesion and Control [1.25215526836087e-09]
%  Lesion and Control [9.72086189321431e-10]