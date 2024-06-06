% how to install python.
% do the matlab tutorial. FIRST INSTALL a version - has to be compatible -
% matlab will like an older version.
% but it still does not register. so, find the executable 
% 
% pyenv(Version="C:\Users\cowen\AppData\Local\Programs\Python\Python310\python.exe")
pyenv(Version="C:\ProgramData\Anaconda3\envs\Phy\python.exe")
% pyenv(Version="C:\Users\cowen\AppData\Local\Programs\Python\Python310\python.exe")
% 
% ans = 
% 
%   PythonEnvironment with properties:
% 
%           Version: "3.10"
%        Executable: "C:\Users\cowen\AppData\Local\Programs\Python\Python310\python.exe"
%           Library: "C:\Users\cowen\AppData\Local\Programs\Python\Python310\python310.dll"
%              Home: "C:\Users\cowen\AppData\Local\Programs\Python\Python310"
%            Status: NotLoaded
%     ExecutionMode: InProcess

% pyversion
% 
%        version: '3.10'
%     executable: 'C:\Users\cowen\AppData\Local\Programs\Python\Python310\python.exe'
%        library: 'C:\Users\cowen\AppData\Local\Programs\Python\Python310\python310.dll'
%           home: 'C:\Users\cowen\AppData\Local\Programs\Python\Python310'
%       isloaded: 0


% You might have to restart a few times - 
% Whn you do, 
py.list({'a'})
% should work.
% pyrun("import pactools as pt") % don't do this - do the next thign...
py.importlib.import_module('pactools')
signal1 = int16(round(1000*randn(10000,1)));
signal1 = signal1(:)'; % make it a row vec

signal2 = py.pactools.simulate_pac(n_points=40000, fs=500, high_fq=100, low_fq=2, low_fq_width=1, noise_level=1);

%.fit(signal)
out = py.pactools.Comodulogram(fs = 500,low_fq_range = 1:1:10,low_fq_width = 1,method = 'tort', n_surrogates = -1, progress_bar = false, n_jobs = 1);
v = out.fit(signal2);
CM = double(v.comod_)';
high_fq_range = double(v.high_fq_range);
low_fq_range = double(v.low_fq_range);

figure
imagesc(low_fq_range,high_fq_range,CM)
axis xy
