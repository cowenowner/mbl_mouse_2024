
:: You can call CatGT three ways:
::
:: 1) > CatGT cmd-line-parameters
:: 2) > runit.bat cmd-line-parameters
:: 3a) Edit parameters in runit.bat, then call it ...
:: 3b) > runit.bat
::
:: This script effectively says:
:: "If there are no parameters sent to runit.bat, call CatGT
:: with the parameters hard coded here, else, pass all of the
:: parameters through to CatGT."
:: F:\Data\Acute_neuropixels\DANA_NAc\Rat_B_unknown\52722_DANA_7000uM_bank0_g0
:: COWEN: docs suggest no -gblcar if extracting LFP - I agree since this is done in Kilosort anyway (kilosort also filters). HOWEVER, -gfix works better with gblcar so that's a thing. Perhaps best to run two separate versions of this code - one for just APs and one for LFPs (with no -gblcar)
::  -xa for digital inputs
:: To make sure a bad channel does not corrupt the -gblcar, then do as follows. SHould be blanked out though.
::  -chnexcl={0;285}
@echo off
@setlocal enableextensions
@cd /d "%~dp0"

set LOCALARGS=-dir=F:\Data\Acute_neuropixels\DANA_NAc\Rat_B_unknown -run=52722_DANA_7000uM_bank0 -g=0 -t=0,570 ^
-prb_fld -t_miss_ok ^
-ni -prb=0 -xa 2,.5

if [%1]==[] (set ARGS=%LOCALARGS%) else (set ARGS=%*)

%~dp0CatGT %ARGS%

