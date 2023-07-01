%%% Modulation of neurons to ketmaine %%%
% Looking at only Saline and ketamine sessions
IX_M1 = categorical(TBL_Q_uS.BrainRegion) == 'M1';
IX_STR = categorical(TBL_Q_uS.BrainRegion) == 'Striatum';
IX_SalKet = categorical(TBL_M1.Drugs) == 'Saline&Ketamine';
% Getting logicals for 6-OHDA lesion unlesion and Control for M1 neuron
% types
TBL = 0;
IX_C = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == 'Control';
IX_Lesion = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R';
IX_UnLesion = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L';
IX_Lesion_LDO = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R';
IX_UnLesion_LDO = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L';
IX_Lesion_LK = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R';
IX_UnLesion_LK = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L';

IX_C_PK = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == 'Control' & TBL_M1.NeuronType == 2;
IX_C_PK100 = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == 'Control' & TBL_M1.NeuronType == 3;
IX_C_IN = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == 'Control' & TBL_M1.NeuronType == 1;
IX_Lesion_PK = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 2;
IX_Lesion_PK100 = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 3;
IX_Lesion_IN = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 1;
IX_UnLesion_PK = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 2;
IX_UnLesion_PK100 = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 3;
IX_UnLesion_IN = categorical(TBL_M1.Drugs) == 'Saline&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 1;

IX_Lesion_LDO_PK = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 2;
IX_Lesion_LDO_PK100 = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 3;
IX_Lesion_LDO_IN = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 1;
IX_UnLesion_LDO_PK = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 2;
IX_UnLesion_LDO_PK100 = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 3;
IX_UnLesion_LDO_IN = categorical(TBL_M1.Drugs) == 'LDOPA&Saline' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 1;

IX_Lesion_LK_PK = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 2;
IX_Lesion_LK_PK100 = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 3;
IX_Lesion_LK_IN = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'R' & TBL_M1.NeuronType == 1;
IX_UnLesion_LK_PK = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 2;
IX_UnLesion_LK_PK100 = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 3;
IX_UnLesion_LK_IN = categorical(TBL_M1.Drugs) == 'LDOPA&Ketamine' & categorical(TBL_M1.Group) == '6ODHA_LID' & categorical(TBL_M1.Hemisphere) == 'L' & TBL_M1.NeuronType == 1;

IX_6OHDA = categorical(TBL_M1.Group) == '6ODHA_LID';

category = [IX_C_PK IX_C_PK100 IX_C_IN IX_Lesion_PK IX_Lesion_PK100 IX_Lesion_IN IX_UnLesion_PK IX_UnLesion_PK100 IX_UnLesion_IN];
% zscore of post 1, 2 & 3 and the baseline for each neuron
TBL_M1 = 0;
for iz = 1:length(IX_C_IN)
    % Frate restricted between -1 and 1 by subtracting post from base and
    % divinding by post + base
%     TBL_M1.FR_post1(iz) = TBL_M1.Frate_post1mbase(iz)/(TBL_M1.Frate_post1(iz)+TBL_M1.Frate_base(iz));
%     TBL_M1.FR_post2(iz) = TBL_M1.Frate_post2mbase(iz)/(TBL_M1.Frate_post2(iz)+TBL_M1.Frate_base(iz));
%     TBL_M1.FR_post3(iz) = TBL_M1.Frate_post3mbase(iz)/(TBL_M1.Frate_post3(iz)+TBL_M1.Frate_base(iz));
    % z score of FR by baseline normalizing
    Q_base = TBL_M1.Q_Base{iz,1};
    Q_base = double(Q_base);
    std_base = std(Q_base);
    mn_base = mean(Q_base);
    Q_all = TBL_M1.Q_All{iz,1};
    Q_all = double(Q_all);
    TBL_M1.z_FR{iz,:} = (Q_all-mn_base)/std_base;
end



SalKet_M1_C_FR = [TBL_M1.FR_post1(IX_C,:) TBL_M1.NeuronType(IX_C,:);TBL_M1.FR_post2(IX_C,:) TBL_M1.NeuronType(IX_C,:);...
    TBL_M1.FR_post3(IX_C,:) TBL_M1.NeuronType(IX_C,:)];
Sal = ones(sum(IX_C),1);
Ket = repmat(2, sum(IX_C),1);
Ket2 = repmat(3, sum(IX_C),1);
Condition = vertcat(Sal, Ket, Ket2);
SalKet_M1_C_FR(:,3) = Condition;
M1_SalKet_C_FR = array2table(SalKet_M1_C_FR,"VariableNames",["Frate","NeuronType","Drug"]);
writetable(M1_SalKet_C_FR,'M1_SalKet_C_FR.csv')

SalKet_M1_C_FR(:,4) = ones(length(SalKet_M1_C_FR),1);
SalKet_M1_C_FR(:,5) = ones(length(SalKet_M1_C_FR),1);
M1_SalKet_C_FR = array2table(SalKet_M1_C_FR,"VariableNames",["Frate","NeuronType","Time","Condition","Group"]);


SalKet_M1_Lesion_FR = [TBL_M1.FR_post1(IX_Lesion,:) TBL_M1.NeuronType(IX_Lesion,:);TBL_M1.FR_post2(IX_Lesion,:) TBL_M1.NeuronType(IX_Lesion,:);...
    TBL_M1.FR_post3(IX_Lesion,:) TBL_M1.NeuronType(IX_Lesion,:)];
Sal = ones(sum(IX_Lesion),1);
Ket = repmat(2, sum(IX_Lesion),1);
Ket2 = repmat(3, sum(IX_Lesion),1);
Condition = vertcat(Sal, Ket, Ket2);
SalKet_M1_Lesion_FR(:,3) = Condition;
M1_SalKet_Lesion_FR = array2table(SalKet_M1_Lesion_FR,"VariableNames",["Frate","NeuronType","Drug"]);
writetable(M1_SalKet_Lesion_FR,'M1_SalKet_Lesion_FR.csv')

SalKet_M1_UnLesion_FR = [TBL_M1.FR_post1(IX_UnLesion,:) TBL_M1.NeuronType(IX_UnLesion,:);TBL_M1.FR_post2(IX_UnLesion,:) TBL_M1.NeuronType(IX_UnLesion,:);...
    TBL_M1.FR_post3(IX_UnLesion,:) TBL_M1.NeuronType(IX_UnLesion,:)];
Sal = ones(sum(IX_UnLesion),1);
Ket = repmat(2, sum(IX_UnLesion),1);
Ket2 = repmat(3, sum(IX_UnLesion),1);
Condition = vertcat(Sal, Ket, Ket2);
SalKet_M1_UnLesion_FR(:,3) = Condition;
M1_SalKet_UnLesion_FR = array2table(SalKet_M1_UnLesion_FR,"VariableNames",["Frate","NeuronType","Drug"]);
writetable(M1_SalKet_UnLesion_FR,'M1_SalKet_UnLesion_FR.csv')
    
LDO_M1_UnLesion_FR = [TBL_M1.FR_post1(IX_UnLesion_LDO,:) TBL_M1.NeuronType(IX_UnLesion_LDO,:);TBL_M1.FR_post2(IX_UnLesion_LDO,:) TBL_M1.NeuronType(IX_UnLesion_LDO,:);...
    TBL_M1.FR_post3(IX_UnLesion_LDO,:) TBL_M1.NeuronType(IX_UnLesion_LDO,:)];
LDO = ones(sum(IX_UnLesion_LDO),1);
Sal = repmat(2, sum(IX_UnLesion_LDO),1);
LDO2 = repmat(3, sum(IX_UnLesion_LDO),1);
Condition = vertcat(LDO, Sal, LDO2);
LDO_M1_UnLesion_FR(:,3) = Condition;
M1_LDO_UnLesion_FR = array2table(LDO_M1_UnLesion_FR,"VariableNames",["Frate","NeuronType","Drug"]);
writetable(M1_LDO_UnLesion_FR,'M1_LDO_UnLesion_FR.csv')

LK_M1_UnLesion_FR = [TBL_M1.FR_post1(IX_UnLesion_LK,:) TBL_M1.NeuronType(IX_UnLesion_LK,:);TBL_M1.FR_post2(IX_UnLesion_LK,:) TBL_M1.NeuronType(IX_UnLesion_LK,:);...
    TBL_M1.FR_post3(IX_UnLesion_LK,:) TBL_M1.NeuronType(IX_UnLesion_LK,:)];
LDO = ones(sum(IX_UnLesion_LK),1);
Ket = repmat(2, sum(IX_UnLesion_LK),1);
Ket2 = repmat(3, sum(IX_UnLesion_LK),1);
Condition = vertcat(LDO, Ket, Ket2);
LK_M1_UnLesion_FR(:,3) = Condition;
M1_LK_UnLesion_FR = array2table(LK_M1_UnLesion_FR,"VariableNames",["Frate","NeuronType","Drug"]);
writetable(M1_LK_UnLesion_FR,'M1_LK_UnLesion_FR.csv')

LK_M1_Lesion_FR = [TBL_M1.FR_post1(IX_Lesion_LK,:) TBL_M1.NeuronType(IX_Lesion_LK,:);TBL_M1.FR_post2(IX_Lesion_LK,:) TBL_M1.NeuronType(IX_Lesion_LK,:);...
    TBL_M1.FR_post3(IX_Lesion_LK,:) TBL_M1.NeuronType(IX_Lesion_LK,:)];
LDO = ones(sum(IX_Lesion_LK),1);
Sal = repmat(2, sum(IX_Lesion_LK),1);
LDO2 = repmat(3, sum(IX_Lesion_LK),1);
Condition = vertcat(LDO, Sal, LDO2);
LK_M1_Lesion_FR(:,3) = Condition;
M1_LK_Lesion_FR = array2table(LK_M1_Lesion_FR,"VariableNames",["Frate","NeuronType","Drug"]);
writetable(M1_LK_Lesion_FR,'M1_LK_Lesion_FR.csv')

LDO_M1_Lesion_FR = [TBL_M1.FR_post1(IX_Lesion_LDO,:) TBL_M1.NeuronType(IX_Lesion_LDO,:);TBL_M1.FR_post2(IX_Lesion_LDO,:) TBL_M1.NeuronType(IX_Lesion_LDO,:);...
    TBL_M1.FR_post3(IX_Lesion_LDO,:) TBL_M1.NeuronType(IX_Lesion_LDO,:)];
Sal = ones(sum(IX_Lesion_LDO),1);
Ket = repmat(2, sum(IX_Lesion_LDO),1);
Ket2 = repmat(3, sum(IX_Lesion_LDO),1);
Condition = vertcat(Sal, Ket, Ket2);
LDO_M1_Lesion_FR(:,3) = Condition;
M1_LDO_Lesion_FR = array2table(LDO_M1_Lesion_FR,"VariableNames",["Frate","NeuronType","Drug"]);
writetable(M1_LDO_Lesion_FR,'M1_LDO_Lesion_FR.csv')

%% Trying to plot time series firing of M1 COntrol neurons during Saline and Ketamine
tiger = TBL_M1.z_FR(IX_C,:);
lion = cell2mat(tiger')';
lioncub = standardize_range(lion',[0 1])';
figure; 
imagesc(length(lioncub(1,:)),[],lioncub)
colorbar

%% Local variance burstiness analysis

M1_Burst_1 = table(categorical(TBL_M1.Group(:,:)), categorical(TBL_M1.Drugs(:,:)), ...
    categorical(TBL_M1.Hemisphere(:,:)), TBL_M1.LocVar_post1mbase(:,:), TBL_M1.NeuronType(:,:), ...
    ones(length(IX_C),1), 'VariableNames', {'Group','Drug', 'Hemisphere', 'LocVar' ,'NeuronType', 'Time'});

M1_Burst_2 = table(categorical(TBL_M1.Group(:,:)), categorical(TBL_M1.Drugs(:,:)), ...
    categorical(TBL_M1.Hemisphere(:,:)), TBL_M1.LocVar_post2mbase(:,:), TBL_M1.NeuronType(:,:), ...
    repmat(2,length(IX_C),1), 'VariableNames', {'Group','Drug', 'Hemisphere', 'LocVar' ,'NeuronType', 'Time'});

M1_Burst_3 = table(categorical(TBL_M1.Group(:,:)), categorical(TBL_M1.Drugs(:,:)), ...
    categorical(TBL_M1.Hemisphere(:,:)), TBL_M1.LocVar_post3mbase(:,:), TBL_M1.NeuronType(:,:), ...
    repmat(3,length(IX_C),1), 'VariableNames', {'Group','Drug', 'Hemisphere', 'LocVar' ,'NeuronType', 'Time'});

M1_Burst = vertcat(M1_Burst_1,M1_Burst_2,M1_Burst_3);
writetable(M1_Burst, 'M1_Burst.csv')


%% Looking at pre and post ketamine frate and LocVar
for iz = 1:length(IX_C_IN)
    % Frate restricted between -1 and 1 by subtracting post from base and
    % divinding by post + base
    TBL_M1.FR_postket(iz) = (TBL_M1.Frate_post2(iz)-TBL_M1.Frate_post1(iz))/(TBL_M1.Frate_post2(iz)+TBL_M1.Frate_post1(iz));
    TBL_M1.LV_postket(iz) = TBL_M1.LocVar_post2(iz)-TBL_M1.LocVar_post1(iz);
end

M1_FR_LV = table(categorical(TBL_M1.Group(:,:)), categorical(TBL_M1.Drugs(:,:)), ...
    categorical(TBL_M1.Hemisphere(:,:)), TBL_M1.LV_postket(:,:), TBL_M1.FR_postket(:,:), TBL_M1.NeuronType(:,:), ...
    'VariableNames', {'Group','Drug', 'Hemisphere', 'LocVar' , 'FRate','NeuronType'});

writetable(M1_FR_LV, 'M1_FR_LV.csv')

%% Look at the post LDO or Post saline period relative to baseline
for iz = 1:length(IX_C_IN)
    % Frate restricted between -1 and 1 by subtracting post from base and
    % divinding by post + base
    TBL_M1.FR_post1mbase_restrict(iz) = (TBL_M1.Frate_post1(iz)-TBL_M1.Frate_base(iz))/(TBL_M1.Frate_post1(iz)+TBL_M1.Frate_base(iz));
end
M1_LID_FR_LV = table(categorical(TBL_M1.Group(IX_6OHDA,:)), categorical(TBL_M1.Drugs(IX_6OHDA,:)), ...
    categorical(TBL_M1.Hemisphere(IX_6OHDA,:)), TBL_M1.LocVar_post1mbase(IX_6OHDA,:), TBL_M1.FR_post1mbase_restrict(IX_6OHDA,:), TBL_M1.NeuronType(IX_6OHDA,:), ...
    'VariableNames', {'Group','Drug', 'Hemisphere', 'LocVar' , 'FRate','NeuronType'});

writetable(M1_LID_FR_LV, 'M1_LID_FR_LV.csv')
%% Plotting ketamine aligned frate

for ii = 1:length(IX_C) 
    
    [close_ket_uS,close_ket_idx]=min(abs(TBL_M1.Q_All_uS{ii,1}-TBL_M1.Inj_uS(ii)));
    t_uS = TBL_M1.Q_All_uS{ii,1};
    Ket_aligned = t_uS - t_uS(close_ket_idx);
    TBL_M1.Ket_align_uS{ii,:} = Ket_aligned((close_ket_idx-120):(close_ket_idx+149));
    Q_all = TBL_M1.Q_All{ii,1};
    TBL_M1.Q_All_ket_align{ii,:} = double(Q_all((close_ket_idx-120):(close_ket_idx+149)));

end


M1_C_ketaligned = cell2mat(TBL_M1.Q_All_ket_align(IX_C,1)')';
figure; imagesc(sort_matrix(M1_C_ketaligned,'kmeans',6))
figure; imagesc(TBL_M1.Ket_align_uS{1,1}/60e6,[],sort_matrix(M1_C_ketaligned))
colorbar
caxis([0 200])
ax = gca;

conv_C_ketalign = conv_filter(M1_C_ketaligned',hanning(20));
diff_conv_C_ketalign = diff(conv_C_ketalign);
figure; imagesc(TBL_M1.Ket_align_uS{1,1}/60e6,[],sort_matrix(diff_conv_C_ketalign'))

M1_C_IN_ketaligned = cell2mat(TBL_M1.Q_All_ket_align(IX_C_IN,1)')';
figure; imagesc(sort_matrix(M1_C_IN_ketaligned))
colorbar
caxis([0 200])

M1_C_PK_ketaligned = cell2mat(TBL_M1.Q_All_ket_align(IX_C_PK,1)')';
figure; imagesc(sort_matrix(M1_C_PK_ketaligned))
colorbar
caxis([0 200])

M1_C_PK100_ketaligned = cell2mat(TBL_M1.Q_All_ket_align(IX_C_PK100,1)')';
figure; imagesc(sort_matrix(M1_C_PK100_ketaligned))
colorbar
caxis([0 200])

M1_Les_ketaligned = cell2mat(TBL_M1.Q_All_ket_align(IX_Lesion,1)')';
figure; imagesc(sort_matrix(M1_Les_ketaligned))
colorbar
caxis([0 200])

M1_UnLes_ketaligned = cell2mat(TBL_M1.Q_All_ket_align(IX_UnLesion,1)')';
figure; imagesc(sort_matrix(M1_UnLes_ketaligned))
colorbar
caxis([0 200])