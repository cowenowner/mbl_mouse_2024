function PHY_Remove_Redundant_Files(data_dir)
if nargin < 1
    data_dir = pwd;
end

delete(fullfile(data_dir,'temp_wh.dat'))
delete(fullfile(data_dir,'time.dat'))
delete(fullfile(data_dir,'_phy_spikes*.npy'))
delete(fullfile(data_dir,'CAR.dat'))
delete(fullfile(data_dir,'event_int16.dat'))

rmdir(fullfile(data_dir,'./.phy'),'s')

phy_dir = fullfile(data_dir,'from_phy');
mkdir(phy_dir)
movefile(fullfile(data_dir,'phy.log'),phy_dir)
movefile(fullfile(data_dir,'*.tsv'),phy_dir)
movefile(fullfile(data_dir,'*.npy'),phy_dir)
movefile(fullfile(data_dir,'params.py'),phy_dir)
movefile(fullfile(data_dir,'rez.mat'),phy_dir)
try
    movefile(fullfile(data_dir,'*ksSettings.mat'),phy_dir)
end

