%% See SPEC_cross_fq_coupling_comod_dupre2017(signal,sFreq,low_range, method) for more detials
% this shows how to implement without my intermediate script and install.
% It uses direct calls to python in matlab rather than going round about.
% it is much more elegant. Much smaller code.
% Requires: 1: use anaconda and get a working version of pactools - verify.
% find the path to that environment and the associated python.exe code.
% put it in using pyenv (I think you only need to do this once).
% then follow the calls just like you would for the real function in
% python.
% you then need to use double() to convert back to matlab variables.
% Currently it only works for python v3.10 or lower.
% SEE https://pactools.github.io/generated/pactools.Comodulogram.html#pactools.Comodulogram
pyenv(Version="C:\ProgramData\Anaconda3\envs\Phy\python.exe")
py.importlib.import_module('pactools')

signal2 = py.pactools.simulate_pac(n_points=40000, fs=500, high_fq=100, low_fq=2, low_fq_width=1, noise_level=3);
signal1 = double(signal2); % just to show that you can just pass in a vector.

out = py.pactools.Comodulogram(fs = 500,low_fq_range = 1:1:10,low_fq_width = 1,method = 'tort', n_surrogates = -1, progress_bar = false, n_jobs = 1);
v = out.fit(signal1);
CM = double(v.comod_)';
high_fq_range = double(v.high_fq_range);
low_fq_range = double(v.low_fq_range);

figure
imagesc(low_fq_range,high_fq_range,CM)
axis xy

