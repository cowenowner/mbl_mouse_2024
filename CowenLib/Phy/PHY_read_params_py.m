function [INFO] = PHY_read_params_py(params_dir)
if nargin < 1
    params_dir = pwd;
end
fname = fullfile(params_dir,'params.py');
fp = fopen(fname,'r');
a = textscan(fp,'%s=%s\n');
fclose(fp);
%INFO = a;
INFO.sFreq = str2double(a{2}{5});