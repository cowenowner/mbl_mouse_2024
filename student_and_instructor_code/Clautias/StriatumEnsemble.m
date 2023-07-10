%% Matias and Claudia

load('E:\NSB_Mouse\23241\01\realtest101_23241_g0\realtest101_23241_g0_imec0\AllSpikes.mat');
wheel_times_sec = load('E:\NSB_Mouse\23241\01\realtest101_23241_g0\synced_synced_synced_realtest101_23241_g0_tcat.nidq.xa_2_0.txt');


T_uS = {SP.t_uS}; % This converts the timestamp list for each cell into a more easy to use cell array.
T_uS_sorted = {};
[sorted, sort_ix] = sort([SP.neuropixels_depth_uM]);

for ii = 1:length(T_uS)
    T_uS_sorted{ii} = T_uS{sort_ix(ii)};
end

figure; 
subplot(1,4,1:3)
plot_raster(T_uS_sorted,[],'time_divisor',60e6);xlabel('minutes')
axis ij
ylabel('Neuron')
subplot(1,4,4)
barh(sorted)
axis ij
xlabel('depth uM')



%% Ensemble analyses
binsz_ms = [20 50 150 25 750 1000] ;

for bin = 1:length(binsz_ms)
[Q, edges_uS] = histcounts_cowen(T_uS_sorted, 'binsize', binsz_ms(bin)*1000);
[corr_TuS_sorted,pval] = corr(Q(1:(10*60*1000/binsz_ms(bin)),:));
meanpval(bin,:) = mean(pval,'omitnan');
signpval(bin,:) = sum(pval < 0.05);
figure; imagesc(corr_TuS_sorted,[0 0.1])
end

legend('20ms','50ms','150ms','250ms','500ms','750ms','1000ms')
xlim([0 0.1])

figure; imagesc(pval,[0 0.1]);


[Q, edges_uS] = histcounts_cowen(T_uS, 'binsize', binsz_ms*1000);
corr_TuS = corr(Q(1:(10*60*1000/binsz_ms),:));
corr_TuS_time = corr(Q(1:(10*60*1000/binsz_ms),:)');
figure; imagesc(corr_TuS)
figure; imagesc(corr_TuS_time)

%% z-score
Z = zscore (Q, 1, 2); %(X,Flag,dim) X= matrix, Flag = Sd, dim = (1= column/ 2= row)
figure; imagesc(Z',[-2 2])
Zcorr = corr(Z'); figure; imagesc(Zcorr)