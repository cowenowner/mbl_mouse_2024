
GIX = sum(AllQ_all)' > 400 & TBL.Depth_uM < 2000 & categorical(TBL.Group) == '6OHDA_LID';
GIXleft = GIX & categorical(TBL.Hemisphere) == 'L';
GIXright = GIX & categorical(TBL.Hemisphere) == 'R';
AllQ_all = double(AllQ_all);
AllQ_b = double(AllQ_b);
AllQ_p1 = double(AllQ_p1);
AllQ_p2 = double(AllQ_p2);


[coeff,score,latent,~,explained] = pca(zscore(AllQ_all(:,GIXright)));

figure
scatter(score(:,1),score(:,2))

figure
scatter((Qall_x_uS/1e6), score(:,1))
xlabel('sec')
ylabel('PC1')
title('PCA over time SHAM-Lesioned hemisphere LID M1 neurons')
pubify_figure_axis
plot_vert_line_at_zero(1410)
plot_vert_line_at_zero(5070)

[~,scoreb, latentb] = pca(zscore(AllQ_b(:,GIXright)));
[~,scorep1, latentp1] = pca(zscore(AllQ_p1(:,GIXright)));
[~,scorep2, latentp2] = pca(zscore(AllQ_p2(:,GIXright)));

figure
plot(latentb)
hold on
plot(latentp1)
hold on 
plot(latentp2)
legend({'b' 'p1' 'p2'})
title('latent values LID right hempishere')
pubify_figure_axis

[~,scorebl, latentbl] = pca(zscore(AllQ_b(:,GIXleft)));
[~,scorep1l, latentp1l] = pca(zscore(AllQ_p1(:,GIXleft)));
[~,scorep2l, latentp2l] = pca(zscore(AllQ_p2(:,GIXleft)));

figure
plot(latentbl)
hold on
plot(latentp1l)
hold on 
plot(latentp2l)
legend({'b' 'p1' 'p2'})
title('latent values LID left hempishere')
pubify_figure_axis


% Stephen's effective dimension code
[Neffb, n1b] = n_effective_dimensions(zscore(AllQ_b(:,GIXright)));
[Neffp1, n1p1]  = n_effective_dimensions(zscore(AllQ_p1(:,GIXright)));
[Neffp2, n1p2]  = n_effective_dimensions(zscore(AllQ_p2(:,GIXright)));

figure
plot(n1b)
hold on
plot(n1p1)
hold on 
plot(n1p2)
legend({'b' 'p1' 'p2'})
title('latent values LID right hempishere')
pubify_figure_axis

[Neffb, n1b] = n_effective_dimensions(zscore(AllQ_b(:,GIXleft)));
[Neffp1, n1p1]  = n_effective_dimensions(zscore(AllQ_p1(:,GIXleft)));
[Neffp2, n1p2]  = n_effective_dimensions(zscore(AllQ_p2(:,GIXleft)));

figure
plot(n1b)
hold on
plot(n1p1)
hold on 
plot(n1p2)
legend({'b' 'p1' 'p2'})
title('latent values LID left hempishere')
pubify_figure_axis


for ii = 1:length(Rat_num)
    nEffb(ii,2) = n_effective_dimensions(zscore(AllQ_b(:,GIX_Rat_left(ii,:))));
    nEffp1(ii,2) = n_effective_dimensions(zscore(AllQ_p1(:,GIX_Rat_left(ii,:))));
    nEffp2(ii,2) = n_effective_dimensions(zscore(AllQ_p2(:,GIX_Rat_left(ii,:))));
    
    nEffb(ii,1) = n_effective_dimensions(zscore(AllQ_b(:,GIX_Rat_right(ii,:))));
    nEffp1(ii,1) = n_effective_dimensions(zscore(AllQ_p1(:,GIX_Rat_right(ii,:))));
    nEffp2(ii,1) = n_effective_dimensions(zscore(AllQ_p2(:,GIX_Rat_right(ii,:))));
    
    if ii < 6
        nEffb(ii,3) = 1;
        nEffp1(ii,3) = 1;
        nEffp2(ii,3) = 1;
    elseif ii > 5
        nEffb(ii,3) = 2;
        nEffp1(ii,3) = 2;
        nEffp2(ii,3) = 2;
    end
end
nEffb = 1;
nEffp2 = 1;
figure
numRows = size(nEffb, 1);
selectedRows = 1:floor(numRows/2); 
selectedRows2 = ceil(numRows/2)+1:numRows;
p1 = signrank(nEffp1(selectedRows,1)-nEffb(selectedRows,1));
mn1 = mean(nEffp1(selectedRows,1)-nEffb(selectedRows,1));

p11 = signrank(nEffp1(selectedRows2,1)-nEffb(selectedRows2,1));
mn11 = mean(nEffp1(selectedRows2,1)-nEffb(selectedRows2,1));

subplot(2,2,1)
gscatter(nEffp1(:,1),nEffb(:,1),nEffb(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline nEff Right Hemipshere')
xlabel('Post 1 nEff')
pubify_figure_axis
% set(gca,'XScale','log')
% set(gca,'YScale','log')
plot_diagonal_line
title(sprintf('LID m1=%1.2f,p1=%1.3f,SHAM m1=%1.2f, p1=%1.3f',mn1,p1,mn11,p11))

p2 = signrank(nEffp2(selectedRows,1)-nEffb(selectedRows,1));
mn2 = mean(nEffp2(selectedRows,1)-nEffb(selectedRows,1));
p22 = signrank(nEffp2(selectedRows2,1)-nEffb(selectedRows2,1));
mn22 = mean(nEffp2(selectedRows2,1)-nEffb(selectedRows2,1));

subplot(2,2,2)
gscatter(nEffp2(:,1),nEffb(:,1),nEffb(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline nEff Right Hemisphere')
xlabel('Post 2 nEff')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
title(sprintf('LID m2=%1.2f,p2=%1.3f, SHAM m2=%1.2f, p2=%1.3f',mn2,p2,mn22,p22))
% Left 

p1 = signrank(nEffp1(selectedRows,2)-nEffb(selectedRows,2));
mn1 = mean(nEffp1(selectedRows,2)-nEffb(selectedRows,2));

p11 = signrank(nEffp1(selectedRows2,2)-nEffb(selectedRows2,2));
mn11 = mean(nEffp1(selectedRows2,2)-nEffb(selectedRows2,2));

subplot(2,2,3)
gscatter(nEffp1(:,2),nEffb(:,2),nEffb(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline nEff Left Hemipshere')
xlabel('Post 1 nEff')
pubify_figure_axis
% set(gca,'XScale','log')
% set(gca,'YScale','log')
plot_diagonal_line
title(sprintf('LID m1=%1.2f,p1=%1.3f,SHAM m1=%1.2f, p1=%1.3f',mn1,p1,mn11,p11))

p2 = signrank(nEffp2(selectedRows,2)-nEffb(selectedRows,2));
mn2 = mean(nEffp2(selectedRows,2)-nEffb(selectedRows,2));
p22 = signrank(nEffp2(selectedRows2,2)-nEffb(selectedRows2,2));
mn22 = mean(nEffp2(selectedRows2,2)-nEffb(selectedRows2,2));

subplot(2,2,4)
gscatter(nEffp2(:,2),nEffb(:,2),nEffb(:,3),'br', 'xo')
axis equal
axis square
ylabel('Baseline nEff Left Hemisphere')
xlabel('Post 2 nEff')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
plot_diagonal_line
pubify_figure_axis
title(sprintf('LID m2=%1.2f,p2=%1.3f, SHAM m2=%1.2f, p2=%1.3f',mn2,p2,mn22,p22))


%% sparness measure

SparsebR = Sparseness_measures(zscore(AllQ_b(:,GIXright)));
Sparsep1R = Sparseness_measures(zscore(AllQ_p1(:,GIXright)));
Sparsep2R = Sparseness_measures(zscore(AllQ_p2(:,GIXright)));

figure;
subplot(3,3,1)
plot(SparsebR.kurtosis)
title('kurtosis')
subplot(3,3,2)
plot(SparsebR.RT)
title('rt')
subplot(3,3,3)
plot(SparsebR.CV)
title('cv')
subplot(3,3,4)
plot(Sparsep1R.kurtosis)
subplot(3,3,5)
plot(Sparsep1R.RT)
subplot(3,3,6)
plot(Sparsep1R.CV)
subplot(3,3,7)
plot(Sparsep2R.kurtosis)
subplot(3,3,8)
plot(Sparsep2R.RT)
subplot(3,3,9)
plot(Sparsep2R.CV)
