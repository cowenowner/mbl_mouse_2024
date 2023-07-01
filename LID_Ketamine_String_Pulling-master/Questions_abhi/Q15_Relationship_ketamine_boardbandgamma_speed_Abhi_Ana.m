%% Analyze data across sessions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
mfile = 'Q15_Relationship_ketamine_boardbandgamma_speed_Abhi_Ana';
ses_to_ana = 'Q15_Relationship_ketamine_boardbandgamma_IMUspeed_Abhi';

GP.Analysis_Dir = 'E:\Temp\TempAnaResults';

PLOT_IT = true;
adir = fullfile(GP.Analysis_Dir,ses_to_ana);
d = dir(fullfile(adir,'Dset*.mat'));
if isempty(d)
    error('wtf')
end
cm = lines(5);

% Initialize
ANA = [];


for iF = 1:length(d)
    Dset = load(fullfile(adir,d(iF).name));
    
    ANA(iF).Rat = Dset.SES.Rat;
    ANA(iF).Rat_type = categorical(Dset.SES.RatType);
    ANA(iF).Session = Dset.SES.Session;
    ANA(iF).Drugs = categorical(Dset.SES.Drugs);
    ANA(iF).ket_gamma_freq = Dset.ket_gamma_filt;
    ANA(iF).LFP_sfreq = Dset.LFP_sfreq;
    ANA(iF).intervals_evt_min = Dset.intervals_around_evt_min;
    ANA(iF).speed_bin_cm_sec = Dset.speed_bin_cm_sec;
    ANA(iF).freq_by_speed_ket1 = Dset.Bin_freq_by_speed_ket1;
    ANA(iF).speed_bins_center_ket1 = Dset.Speed_bin_center_ket1;
    ANA(iF).freq_by_speed_ket2 = Dset.Bin_freq_by_speed_ket2;
    ANA(iF).speed_bins_center_ket2 = Dset.Speed_bin_center_ket2;
    ANA(iF).power_by_speed_ket1 = Dset.Bin_power_by_speed_ket1;
    ANA(iF).power_by_speed_ket2 = Dset.Bin_power_by_speed_ket2;

    
end

TBL = struct2table(ANA);

for ii = 1:size(TBL)
    if TBL.Drugs(ii) == 'Saline&Ketamine'
        figure
        plot(TBL.speed_bins_center_ket1{ii}, TBL.freq_by_speed_ket1{ii})
        title(sprintf('Rat %g Session# %d Sal&Ket', TBL.Rat(ii), TBL.Session(ii)), 'Fontsize', 10)
    end    
end
cnt352 = 1;
cnt350 = 1;
cnt342 = 1;
cnt320 = 1;

for ii = 1:size(TBL)
    if TBL.Drugs(ii) == 'Saline&Ketamine' && TBL.Rat(ii) == 352
% %         figure
%         plot(TBL.speed_bins_center_ket1{ii}, TBL.freq_by_speed_ket1{ii})
%         hold on
%         title('Rat 352 Sal&Ket', 'Fontsize', 10)
        [beta352(cnt352),~,r352{cnt352}] = regress(TBL.freq_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt352 = cnt352 + 1;
        
    elseif TBL.Drugs(ii) == 'Saline&Ketamine' && TBL.Rat(ii) == 350
        [beta350(cnt350),~,r350{cnt350}] = regress(TBL.freq_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt350 = cnt350 + 1;
    elseif TBL.Drugs(ii) == 'Saline&Ketamine' && TBL.Rat(ii) == 342
        [beta342(cnt342),~,r342{cnt342}] = regress(TBL.freq_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt342 = cnt342 + 1;
    elseif TBL.Drugs(ii) == 'Saline&Ketamine' && TBL.Rat(ii) == 320
        [beta320(cnt320),~,r320{cnt320}] = regress(TBL.freq_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt320 = cnt320 + 1;
    end    
end
[h_SK, p_SK] = ttest([mean(beta320) mean(beta342) mean(beta350) mean(beta352)]);

cnt352 = 1;
cnt350 = 1;
cnt342 = 1;
cnt320 = 1;
for ii = 1:size(TBL)
    if TBL.Drugs(ii) == 'LDOPA&Ketamine' && TBL.Rat(ii) == 352
% %         figure
%         plot(TBL.speed_bins_center_ket1{ii}, TBL.freq_by_speed_ket1{ii})
%         hold on
%         title('Rat 352 Sal&Ket', 'Fontsize', 10)
        [beta352_LK{cnt352},~,r352_LK{cnt352},~,stats352{cnt352}] = regress(TBL.freq_by_speed_ket1{ii}',[ones(length(TBL.speed_bins_center_ket1{ii}),1) TBL.speed_bins_center_ket1{ii}']);
        cnt352 = cnt352 + 1;
        
    elseif TBL.Drugs(ii) == 'LDOPA&Ketamine' && TBL.Rat(ii) == 350
        [beta350_LK{cnt350},~,r350_LK{cnt350},~,stats350{cnt350}] = regress(TBL.freq_by_speed_ket1{ii}',[ones(length(TBL.speed_bins_center_ket1{ii}),1) TBL.speed_bins_center_ket1{ii}']);
        cnt350 = cnt350 + 1;
    elseif TBL.Drugs(ii) == 'LDOPA&Ketamine' && TBL.Rat(ii) == 342
        [beta342_LK{cnt342},~,r342_LK{cnt342},~,stats342{cnt342}] = regress(TBL.freq_by_speed_ket1{ii}',[ones(length(TBL.speed_bins_center_ket1{ii}),1) TBL.speed_bins_center_ket1{ii}']);
        cnt342 = cnt342 + 1;
    elseif TBL.Drugs(ii) == 'LDOPA&Ketamine' && TBL.Rat(ii) == 320
        [beta320_LK{cnt320},~,r320_LK{cnt320},~,stats320{cnt320}] = regress(TBL.freq_by_speed_ket1{ii}',[ones(length(TBL.speed_bins_center_ket1{ii}),1) TBL.speed_bins_center_ket1{ii}']);
        cnt320 = cnt320 + 1;
    end    
end
beta320_LK = cell2mat(beta320_LK);
beta342_LK = cell2mat(beta342_LK);
beta350_LK = cell2mat(beta350_LK);
beta352_LK = cell2mat(beta352_LK);

[h_LK, p_LK] = ttest([mean(beta320_LK(2,:)) mean(beta342_LK(2,:)) mean(beta350_LK(2,:)) mean(beta352_LK(2,:))]);

beta320_AK1 = cat(2,beta320, beta320_LK);
beta342_AK1 = cat(2,beta342, beta342_LK);
beta350_AK1 = cat(2,beta350, beta350_LK);
beta352_AK1 = cat(2,beta352, beta352_LK);

[h_AK1, p_AK1] = ttest([mean(beta320_AK1) mean(beta342_AK1) mean(beta350_AK1) mean(beta352_AK1)]);
d_AK1 = Cohens_d([mean(beta320_AK1) mean(beta342_AK1) mean(beta350_AK1) mean(beta352_AK1)]);

for ii = 1:size(TBL)
    if TBL.Rat(ii) == 320
% %         figure
        plot(TBL.speed_bins_center_ket2{ii}, TBL.freq_by_speed_ket2{ii})
        hold on
        title('Rat 320 All Ket second half', 'Fontsize', 10)
    end
end
cnt352 = 1;
cnt350 = 1;
cnt342 = 1;
cnt320 = 1;
for ii = 1:size(TBL)
    if TBL.Rat(ii) == 352
        [beta352_AK2(cnt352),~,r352_AK2{cnt352}] = regress(TBL.freq_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt352 = cnt352 + 1;
        
    elseif TBL.Rat(ii) == 350
        [beta350_AK2(cnt350),~,r350_AK2{cnt350}] = regress(TBL.freq_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt350 = cnt350 + 1;
    elseif TBL.Rat(ii) == 342
        [beta342_AK2(cnt342),~,r342_AK2{cnt342}] = regress(TBL.freq_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt342 = cnt342 + 1;
    elseif TBL.Rat(ii) == 320
        [beta320_AK2(cnt320),~,r320_AK2{cnt320}] = regress(TBL.freq_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt320 = cnt320 + 1;
    end    
end

[h_AK2, p_AK2] = ttest([mean(beta320_AK2) mean(beta342_AK2) mean(beta350_AK2) mean(beta352_AK2)]);
d_AK2 = Cohens_d([mean(beta320_AK2) mean(beta342_AK2) mean(beta350_AK2) mean(beta352_AK2)]);

[h_AK2mAK1, p_AK2mAK1] = ttest([mean(beta320_AK2)-mean(beta320_AK1) ...
    mean(beta342_AK2)-mean(beta342_AK1) ...
    mean(beta350_AK2)-mean(beta350_AK1) ... 
    mean(beta352_AK2)-mean(beta352_AK1)]);
d_AK2mAK1 = Cohens_d([mean(beta320_AK2)-mean(beta320_AK1) ...
    mean(beta342_AK2)-mean(beta342_AK1) ...
    mean(beta350_AK2)-mean(beta350_AK1) ... 
    mean(beta352_AK2)-mean(beta352_AK1)]);

for ii = 1:size(TBL)
    if TBL.Drugs(ii) == 'LDOPA&Ketamine' && TBL.Rat(ii) == 320
%         figure
        plot(TBL.speed_bins_center_ket1{ii}, TBL.freq_by_speed_ket1{ii})
        hold on
        title('Rat 320 LDO&Ket', 'Fontsize', 10)
    end    
end

for ii = 1:size(TBL)
    if TBL.Rat(ii) == 350
%         figure
        plot(TBL.speed_bins_center_ket1{ii}, TBL.power_by_speed_ket1{ii})
        hold on
        title('Rat 350 All Ket', 'Fontsize', 10)
    end    
end

cnt352 = 1;
cnt350 = 1;
cnt342 = 1;
cnt320 = 1;
for ii = 1:size(TBL)
    if TBL.Rat(ii) == 352
        [pow_beta352_AK1(cnt352),~,pow_r352_AK1{cnt352}] = regress(TBL.power_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt352 = cnt352 + 1;
        
    elseif TBL.Rat(ii) == 350
        [pow_beta350_AK1(cnt350),~,pow_r350_AK1{cnt350}] = regress(TBL.power_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt350 = cnt350 + 1;
    elseif TBL.Rat(ii) == 342
        [pow_beta342_AK1(cnt342),~,pow_r342_AK1{cnt342}] = regress(TBL.power_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt342 = cnt342 + 1;
    elseif TBL.Rat(ii) == 320
        [pow_beta320_AK1(cnt320),~,pow_r320_AK1{cnt320}] = regress(TBL.power_by_speed_ket1{ii}',TBL.speed_bins_center_ket1{ii}');
        cnt320 = cnt320 + 1;
    end    
end

[pow_h_AK1, pow_p_AK1] = ttest([mean(pow_beta320_AK1) mean(pow_beta342_AK1) mean(pow_beta350_AK1) mean(pow_beta352_AK1)]);
pow_d_AK1 = Cohens_d([mean(pow_beta320_AK1) mean(pow_beta342_AK1) mean(pow_beta350_AK1) mean(pow_beta352_AK1)]);

cnt352 = 1;
cnt350 = 1;
cnt342 = 1;
cnt320 = 1;
for ii = 1:size(TBL)
    if TBL.Rat(ii) == 352
        [pow_beta352_AK2(cnt352),~,pow_r352_AK2{cnt352}] = regress(TBL.power_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt352 = cnt352 + 1;
        
    elseif TBL.Rat(ii) == 350
        [pow_beta350_AK2(cnt350),~,pow_r350_AK2{cnt350}] = regress(TBL.power_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt350 = cnt350 + 1;
    elseif TBL.Rat(ii) == 342
        [pow_beta342_AK2(cnt342),~,pow_r342_AK2{cnt342}] = regress(TBL.power_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt342 = cnt342 + 1;
    elseif TBL.Rat(ii) == 320
        [pow_beta320_AK2(cnt320),~,pow_r320_AK2{cnt320}] = regress(TBL.power_by_speed_ket2{ii}',TBL.speed_bins_center_ket2{ii}');
        cnt320 = cnt320 + 1;
    end    
end

[pow_h_AK2, pow_p_AK2] = ttest([mean(pow_beta320_AK2) mean(pow_beta342_AK2) mean(pow_beta350_AK2) mean(pow_beta352_AK2)]);
pow_d_AK2 = Cohens_d([mean(pow_beta320_AK2) mean(pow_beta342_AK2) mean(pow_beta350_AK2) mean(pow_beta352_AK2)]);