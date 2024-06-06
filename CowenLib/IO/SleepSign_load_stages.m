function T = SleepSign_load_stages(fname)
if ~exist(fname,'file')
    fname
    error('Could not find stages file')
end
% opts = detectImportOptions(fname);
% opts = setvaropts(opts,'VariableOptions','InputFormat','MM/dd/uuuu HH:mm:ss')
% It is in MM/dd/yyyy format but not clear how to force this in readtable.
%'Format','%{dd/MM/uuuu HH:mm}D%f %f %f %f %f %f %f %f %f'
T = readtable(fname,'NumHeaderLines',17,'Delimiter', ',','Format','%{MM/dd/uuuu HH:mm:ss}D%f%f%f%f');
% it's in Day/Month/year
T.original_t_sec = T.Time.Day*24*60*60 + T.Time.Hour*60*60 + T.Time.Minute*60 + T.Time.Second;
