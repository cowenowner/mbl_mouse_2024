# mbl_mouse_2024
Code for the Woods Hole MBL 2024 mouse course. In inherits a lot from the previous 2023 course.

## INSTALL:
(students should not have to worry about this as it should only need to be done once per computer)

Create a C:\Temp\Analysis_results folder

Analysis computers should have the following installed...

GitHub Desktop
Anaconda
phy using the environment.yml file on the cortex/phy github site (takes about 1hr) - run as administrator
Matlab
Kilosort2.5 along with Visual Studio C++ 2019 professional - just install the Desktop C++ development version.

Github for : mbl_mouse_2024

Put these in the C:\ folder (not a subfolder)
SpikeGLX
TCat
Tprime

MATLAB: Create a desktop shortcut for Matlab and change it so that it start in the GitHub directory (e.g., C:\Users\Administrator\Documents\GitHub\mbl_mouse_2024). The specific directory may differ by computer.


# FUN EXCITING DEMOS!

## Demo 1: Ensembles
How to run Stephen's demo: 
1. Start up matlab. Make sure CowenLib and subfolders are in your path.
2. Go to the DEMOS\Ensemble_Demo_1 in Matlab. 
3. Open the DEMO1_does_MFB_stim_affect_ensemble_act.m file and follow the instructions in the file.

## Demo 2: Oscillations
How to run Abhi's demo on how to analyze brain oscillations associated with Parkinson's disease and levodopa-induced dyskinesia.

1. Open the MATLAB icon named Abhi_Demo: This will load the paths you will need to run the demo.
2. Open the script called Demo_NSB: You can do this by typing "edit Demo_NSB" on the command line or double click on the matlab file that can be found here C:\Users\Administrator\Documents\GitHub\mbl_mouse_2024\LID_Ketamine_String_Pulling-master
3. Change the current folder (the top left window in MATLAB most of the time) by clicking on the icon that looks like a folder with a little green arrow (right above the current folder). Change the current folder to the data directory which is: E:\Data_Demo_Abhi\09
4. Now just run the Demo_NSB script(Green arrow or command line))! It will generate some figures showing dyskinesia associated and ketamien associated oscillations in the spectrogram and power spectral density figures.

