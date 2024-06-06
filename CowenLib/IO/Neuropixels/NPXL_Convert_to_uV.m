function data = NPXL_Convert_to_uV(data,meta)
% Convert to uV...
data = double(data);
gain=str2double(meta.imChan0apGain);
imax=str2double(meta.imMaxInt);
vmax=str2double(meta.imAiRangeMax);
data = 1e6*data*(vmax/(imax*gain));