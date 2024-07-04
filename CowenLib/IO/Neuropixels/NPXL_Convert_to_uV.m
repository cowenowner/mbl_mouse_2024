function data = NPXL_Convert_to_uV(data,meta)
% Convert to uV...
data = double(data);
try
    gain=str2double(meta.imChan0apGain);
catch
    gain=str2double(meta.imChan0lfGain);

end
imax=str2double(meta.imMaxInt);
vmax=str2double(meta.imAiRangeMax);
data = 1e6*data*(vmax/(imax*gain));