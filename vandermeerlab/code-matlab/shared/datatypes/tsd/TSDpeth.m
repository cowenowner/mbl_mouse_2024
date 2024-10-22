function out = tsdPETH(cfg_in,tsd_in,t_in)
% function out = tsdPETH(cfg_in,tsd_in,t_in)
%
% compute peri-event time histogram (average)
%
% INPUTS:
%
% tsd_in: tsd to compute average over
% t_in: times to average: can be raw timestamps, ts or iv data. If iv, cfg.window is ignored
%
% OUTPUTS:
%
% out: tsd with PETH
%
% cfg options:
%
% cfg_def.window = [-2 2]; % start and end times of window (in s)
% cfg_def.dt = []; % time step, specify this for 'interp' mode
% cfg_def.mode = 'raw'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
% cfg_def.interp_mode = 'linear';
%
% EXAMPLE USAGE:
%
% cd('R042-2013-08-18'); LoadExpKeys; LoadMetadata;
% cfg = []; cfg.fc = ExpKeys.goodSWR(1);
% lfp = LoadCSC(cfg);
% peth = TSDpeth([],lfp,metadata.taskvars.trial_iv.tstart(2:end));
% plot(peth);
%
% NOTES:
%
% 'raw' mode collects tsd samples that fall within the specified intervals.
% Because it does not interpolate, it can easily fail if there are gaps in
% the data or the diffs between samples are unequal.
%
% 'interp' interpolates on a fixed timebase, so guarantees the number of
% samples will be constant.
%
% multidimensional TSDs are not yet supported!
%
% MvdM 2017-08-16 initial version

cfg_def = [];
cfg_def.window = [-2 2]; % start and end times of window (in s)
cfg_def.dt = []; % time step, specify this for 'interp' mode
cfg_def.mode = 'raw'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
cfg_def.interp_mode = 'linear';

cfg = ProcessConfig2(cfg_def,cfg_in);

if ~CheckTSD(tsd_in)
    error('Incorrectly formed tsd input.');
end

nDim = sum(size(tsd_in.data) > 1);

if nDim > 2
    error('More than 2 dimensions in TSD input not yet supported.');
end

if ~isfield(t_in,'type') % assume raw times
    this_iv = iv(t_in + cfg.window(1), t_in + cfg.window(2));
else
    switch t_in.type
        case 'iv'
            
            if ~CheckIV(t_in) % this should check that there are no t_end before t_start, etc...
                error('Incorrectly formed iv input.');
            end
            
            this_iv = t_in;
               
        case 'ts'
            
            if length(t_in.t) ~= 1
                error('ts input must have exactly one .t cell');
            end
            
            this_iv = iv(t_in.t{1} + cfg.window(1), t_in.t{1} + cfg.window(2));
            
        otherwise
            
            error('t input must be raw times, iv or ts data')
        
    end
end

% check that iv's don't exceed bounds of data
early_PETH_idx = this_iv.tstart < tsd_in.tvec(1);
if sum(early_PETH_idx) > 0
    warning('%d preceding PETH(s) removed.', sum(early_PETH_idx));
    this_iv.tstart = this_iv.tstart(~early_PETH_idx);
    this_iv.tend = this_iv.tend(~early_PETH_idx);
end

late_PETH_idx = this_iv.tend > tsd_in.tvec(end);
if sum(late_PETH_idx) > 0
    warning('%d exceeding PETH(s) removed.', sum(late_PETH_idx));
    this_iv.tstart = this_iv.tstart(~late_PETH_idx);
    this_iv.tend = this_iv.tend(~late_PETH_idx);
end

nT = length(this_iv.tstart);



switch cfg.mode
    
    case 'raw'
        
        % make big idx matrix
        start_idx = nearest_idx3(this_iv.tstart,tsd_in.tvec);
        end_idx = nearest_idx3(this_iv.tend,tsd_in.tvec);
        
        if length(unique(end_idx-start_idx)) ~= 1 % unequal length trials
           error('Raw mode requires equal tsd samples for each trial.');
        end
        
        % 1-D
        switch nDim
            case 1
                for iT = nT:-1:1 % slow! could be vectorized
                    out_data(iT,:) = tsd_in.data(start_idx(iT):end_idx(iT));
                end

                % 2-D
            case 2
                for iT = nT:-1:1 % slow! could be vectorized
                    out_data(iT,:,:) = tsd_in.data(:, start_idx(iT):end_idx(iT))';
                end
        end
        out_tvec = tsd_in.tvec(start_idx(1):end_idx(1));
        out_tvec = out_tvec-nanmean(out_tvec); % this is an approximation -- depends on exact spacing of input tsd
                
    case 'interp'
        
        if isempty(cfg.dt)
            error('interp mode requires cfg.dt to be specified.');
        end
        
        for iT = nT:-1:1 % slow! should be vectorized to do one interp1 first, then reshape
            
            out_data(iT,:) = interp1(tsd_in.tvec,tsd_in.data,this_iv.tstart(iT):cfg.dt:this_iv.tend(iT),cfg.interp_mode); % need to generalize to 2-D
            
        end
        out_tvec = cfg.window(1):cfg.dt:cfg.window(2);
          
    otherwise
        
        error('unknown cfg.mode %s',cfg.mode);
end

% average and package
out = tsd;
out.data = sq(nanmean(out_data, 1)); if nDim == 2, out.data = out.data'; end
out.tvec = out_tvec;