function OUT = Numbers_from_filename_Abhi(filename)

if iscell(filename)
    for ii = 1:length(filename)
        numbers = [];
        numbers = regexp(filename{ii}, '\d+', 'match');
        OUT(ii,:) = cellfun(@str2double, numbers{1});
    end
else
    error('input should be a cell conatining a string')
end
        


