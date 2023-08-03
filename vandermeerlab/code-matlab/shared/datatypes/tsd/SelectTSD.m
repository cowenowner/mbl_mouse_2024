function [tsd_out, idx] = SelectTSD(cfg_in, tsd_in, selectspec)
%SELECTTSD Specify time points to keep.
%   [tsd_out, idx] = SelectTSD(cfg_in, tsd_in, selectspec)
%
%   INPUTS:
%         cfg: config struct with fields controlling function behavior
%       tsd_in: tsd struct
%  selectspec: selection specifics, either:
%         - [nx1] double: logical array or indices specifying which
%                 time points to keep.
%         - string: string specifying which usr field to work with. If
%                 selectspec is a string, the config options cfg.operation,
%                 cfg.threshold, and cfg.str apply.
%
%
%   OUTPUTS
%      tsd_out: tsd struct with time points selected and all
%              corresponding same-length usr trimmed accordingly
%      idx: indices into original tsd_in.tvec of kept time points
%
%   CFG OPTIONS
%       cfg.operation = '>='; How to perform numerical selection, see
%                     cfg.threshold.
%                '>' - usr data > threshold
%               '>=' - usr data >= threshold
%                '<' - usr data < threshold
%               '<=' - usr data <= threshold
%                '=' - usr data = threshold
%       cfg.threshold = 0; Set a numerical threshold for keeping time points.
%                     This works on numerical usr contents, but can also be
%                     applied to strings as long as the first character is
%                     number-convertible:
%                     If your field contains strings and cfg.str is
%                     empty, SelectTSD assumes that the first character is a
%                     number (i.e. a rating) and thresholds based on this 
%                     number. An example would be '1, very good', for which
%                     SelectTSD considers the 1 only. 
%       cfg.str = ''; If your target usr field contains strings that 
%                     are NOT number-convertible, input the string you
%                     want to select by. If this is not empty, it overrides 
%                     numerical selection. Examples of non-number-convertible 
%                     strings might be 'good' or 'maybe' or 'poor'.
%       cfg.verbose = 1; Tell me how many time points came in, and how many
%                     went out.
%
% MvdM Aug 2023, based on aacarey's SelectIV()
%
% see also SelectIV, restrict
cfg_def.operation = '>=';
cfg_def.threshold = 0;
cfg_def.str = ''; % if this is not empty, it overrides numerical selection
cfg_def.verbose = 1;
mfun = mfilename;

if ~CheckTSD(tsd_in,mfun)
    error('tsd_in must be a tsd data type.')
end

% parse cfg parameters
cfg = ProcessConfig(cfg_def, cfg_in, mfun);

if isempty(tsd_in.tvec)
    if cfg.verbose
        fprintf('%s: tsd_in is empty, returning iv_in\n', mfun)
        tsd_out = History(tsd_in, mfun, cfg);
        idx = [];
        return
    end
end

% choose which thing to do
if islogical(selectspec) || isnumeric(selectspec)
    spec_type = 'log_or_num'; % specifying intervals to keep using a logical or numerical array
elseif ischar(selectspec)
    spec_type = 'string'; % specifying a string which corresponds to a usr field name
else
    error('selectspec must be a logical array, numeric array of indices, or a string specifying a usr field name.')
end

switch spec_type
    case 'string'  
        % make sure usr exists
        if ~isfield(tsd_in,'usr')
            error('tsd_in requires usr for this type of selection.')
        end
        % check that the field actually exists and that it's the right length
        if ischar(selectspec) && ~isfield(tsd_in.usr,selectspec)
            error([selectspec,' does not exist.'])
        elseif ischar(selectspec) && length(tsd_in.usr.(selectspec)) ~= length(tsd_in.tvec)
            error(['tsd_in.usr.',selectspec,' must have the same dimensions as tsd_in.tvec.'])
        end
        
        % if the field contains strings, get ratings in numerical form
        if isempty(cfg.str) && ~isnumeric(tsd_in.usr.(selectspec)(1))
            str_type = 'rating'; % something like '5, or delete'
            temp = nan(size(tsd_in.usr.(selectspec)));
            for ii = 1:length(temp)
                temp(ii,1) = str2double(tsd_in.usr.(selectspec){ii,1}(1)); % we assume the rating is the first character in the string
            end
        elseif isempty(cfg.str) && isnumeric(tsd_in.usr.(selectspec))
            str_type = 'rating'; % something like '5, or delete'
            temp = tsd_in.usr.(selectspec);
        elseif ~isempty(cfg.str)
            str_type = 'description'; % something like 'good'
            temp = tsd_in.usr.(selectspec);
        end
        
        % do the thing
        switch str_type
            case 'rating'
                switch cfg.operation
                    case '>'
                        keep =  temp > cfg.threshold;
                    case '>='
                        keep = temp >= cfg.threshold;
                    case '<'
                        keep = temp < cfg.threshold;
                    case '<='
                        keep = temp <= cfg.threshold;
                    case '='
                        keep = temp == cfg.threshold;
                    otherwise
                        error('Unrecognized cfg.operation')
                end
            case 'description'
                keep = nan(size(temp));
                for iStr = 1:length(temp)
                    keep(iStr) = strcmp(cfg.str,temp(iStr));
                end
        end
        
        keep = logical(keep);
        
    case 'log_or_num'
        keep = selectspec;
        % these config options do not apply in this case, so don't give
        % them a value in history since they did not affect the output
        cfg.operation = '';
        cfg.threshold = [];
        cfg.str = '';
end

tsd_out = tsd_in;
tsd_out.tvec = tsd_out.tvec(keep);
tsd_out.data = tsd_out.data(:, keep);

% also select data from other same-length usr fields
if isfield(tsd_out,'usr') && ~isempty(tsd_out.usr)
    fields = fieldnames(tsd_out.usr);
    for iField = 1:length(fields)
        tsd_out.usr.(fields{iField}) = tsd_out.usr.(fields{iField})(keep);
    end
end

% make idx output
if islogical(keep)
    idx = find(keep);
elseif isnumeric(keep)
    idx = keep; 
end

% talk to me
if cfg.verbose
    disp([mfun,': ',num2str(length(tsd_in.tvec)),' time points in, ',num2str(length(tsd_out.tvec)),' time points out.'])
end

% housekeeping
tsd_out.cfg.history.mfun = cat(1,tsd_in.cfg.history.mfun,mfun);
tsd_out.cfg.history.cfg = cat(1,tsd_in.cfg.history.cfg,{cfg});

end % of function