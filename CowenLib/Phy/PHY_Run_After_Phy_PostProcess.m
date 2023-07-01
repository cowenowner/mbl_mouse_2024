function SP = PHY_Run_After_Phy_PostProcess(phy_dir)
if nargin == 0
    phy_dir = pwd;
end
PHY_Process_Sorted_Neurons_to_AllSpikes(phy_dir)
PHY_Remove_Redundant_Files(phy_dir)
