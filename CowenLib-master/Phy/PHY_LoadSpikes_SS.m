sFreq=30000;
groups_to_load = {'good','mua'};

data_dir='F:\Data\vHPC_stim_mPFC_excitability_project\10824\mPFC_L23_2_tetrodeconfig_g0\mPFC_L23_2_tetrodeconfig_g0_imec0';
[SP,~] = PHY_read_spikes(data_dir,sFreq, groups_to_load);
filename=sprintf('%s\\AllSpikes.mat',data_dir);

save(filename,'SP')











F:\Data\vHPC_stim_mPFC_excitability_project\Rat10805\mPFC1_L5_vHPCstim_tetrodeconfig_g0\mPFC1_L5_vHPCstim_tetrodeconfig_g0_imec0