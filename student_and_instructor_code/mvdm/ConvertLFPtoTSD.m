function myLFP = ConvertLFPtoTSD(LFP, ch_in_list)
% function myLFP = ConvertLFPtoTSD(LFP, ch_in_list)
%
% convert NPXL LFP as loaded by obj.ReadBinVolts to vandermeerlab tsd data
% structure

myLFP = tsd;
myLFP.tvec = LFP.tvec;
myLFP.data = LFP.data(ch_in_list, :);
for iLFP = 1:size(myLFP.data, 1)
    myLFP.cfg.hdr{iLFP}.Fs = LFP.sFreq;
end