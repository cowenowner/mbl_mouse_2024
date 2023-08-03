function [muaSDF,muaS] = getMUA(varargin)
% function [muaSDF,muaS] = getMUA(varargin)
%
% outputs multiunit spike density function
%
% varargins:
% target_dir = pwd;
% sigma = 0.025; % sd for spike convolution
% dt = 0.01; % output resolution in seconds
%
% note gaussian convolution window is limited to 1s
%
% MvdM 10

target_dir = pwd;
sigma = 0.025; % sd for spike convolution
dt = 0.01; % output resolution in seconds

extract_varargin;

cd(target_dir);
sd = mtarlinit('.','allCells',1,'fullTime',1,'SET_MinSpks',10); sd.S = sd.fullS; % should standardize

% get start time from some CSC
csc_fn = sd.ExpKeys.goodSWR{1};
disp(sprintf('Loading CSC %s to get times...',csc_fn));
csc = LoadCSC(csc_fn,'ConvertToS',1);
t0 = StartTime(csc); t1 = EndTime(csc);
clear csc;

% initialize
t0 = dt.*(floor(t0./dt));
t1 = dt.*(ceil(t1./dt));
tvec = t0:dt:t1; tvec = tvec';
mua = zeros(size(tvec));

% remove interneurons
celltype = TwoTcellTypes(sd);
sd.S = sd.S(celltype == 5);

% add spikes
for iC = 1:length(sd.S)
   spks = Data(Restrict(sd.S{iC},t0,t1)); 
   spk_bin = ceil((spks-t0)./dt);
   mua(spk_bin) = mua(spk_bin)+1;
end

% convolve
flt = gausswin(round(1/dt),round(1/sigma)); % 1s gaussian window with SD sigma
flt = flt./sum(flt);
mua = conv2(mua,flt,'same');

muaSDF = tsd(tvec,mua);

% get spikes
muaS = [];
for iC = 1:length(sd.S)
    muaS = cat(1,muaS,Data(sd.S{iC}));
end
muaS = ts(sort(muaS,'ascend'));