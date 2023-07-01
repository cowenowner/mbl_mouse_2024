function OUT = Q1_Do_neurons_correlate_with_string_pulling_single_unit()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do individual neurons correlate with features of string pulling.
% Calculate a metric for individual units. See the ensemble version for an
% ensemble-level metric.
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT = [];
GP = LK_Globals;
PLOT_IT = true;
Behavior_sFreq = 100;
neuron_quality_threshold = 2;
small_binsize_ms = 5;
params = 25:50:200; % smooth
params = 50; % smooth
PETH_win_ms = 300;

vbls = {'Right_speed' 'Left_speed' 'Nose_speed' 'Left_acc' 'Right_acc'  'Right_x_to_nose' 'Left_x_to_nose' 'Rot_speed' 'Rot_acc' 'IMU_speed' 'Right_y_d1' 'Left_y_d1'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SES = LK_Session_Info();
OUT.SES = SES;
OUT.aborted = false;
OUT.smooth_win_ms = params;
OUT.variables = vbls;
OUT.small_binsize_ms = small_binsize_ms;
OUT.Behavior_sFreq = Behavior_sFreq;
OUT.neuron_quality_threshold = neuron_quality_threshold;

% Estimate the individual pull times.
PULL = LK_Determine_Individual_Pull_Times();
% Do PETHs on these data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and create a magnificent table of motion data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[MT, string_pull_intervals_uSec,~,~,~,~,EVENTS] = LK_Combine_All_String_Pull_Motion_To_Table(Behavior_sFreq, false);
epoch_st_ed_uSec = [string_pull_intervals_uSec(1) - 2e6 string_pull_intervals_uSec(end) + 2e6 ]; % add a little padding.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load spikes and restrict to the times of interest.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SP,TS_uS] = LK_Load_Spikes(neuron_quality_threshold,epoch_st_ed_uSec);
OUT.n_neurons = length(SP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start with classic PETHs aligned to key movement events.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[OUT.PTH_L_st_up,OUT.PTH_x_axis,~,OUT.PTH_L_st_up_p] = PETH_raster_array(TS_uS,EVENTS.Left_start_up_t_uS, 4, PETH_win_ms, PETH_win_ms);
[OUT.PTH_R_st_up,OUT.PTH_x_axis,~,OUT.PTH_R_st_up_p] = PETH_raster_array(TS_uS,EVENTS.Right_start_up_t_uS, 4, PETH_win_ms, PETH_win_ms);

[OUT.PTH_L_st_down,OUT.PTH_x_axis,~,OUT.PTH_L_st_down_p] = PETH_raster_array(TS_uS,EVENTS.Left_start_down_t_uS, 4, PETH_win_ms, PETH_win_ms);
[OUT.PTH_R_st_down,OUT.PTH_x_axis,~,OUT.PTH_R_st_down_p] = PETH_raster_array(TS_uS,EVENTS.Right_start_down_t_uS, 4, PETH_win_ms, PETH_win_ms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a space of rising or falling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert each pull epoch into a space from 0 to 1 that spans the start and
% end of a pull. We could also restrict other analyses such as responses to
% pull speed to just these intervals.
smooth_win_ms = 50;
smooth_bins = round(smooth_win_ms/small_binsize_ms);

[NT,~,NT_t_uS] = Bin_and_smooth_ts_array(TS_uS, small_binsize_ms*1e3, smooth_bins);
res = 100;
PE.Left.up.spk = nan(res,Cols(NT),length(EVENTS.Left_start_up_t_uS));

SE = [];
SE{1} = [EVENTS.Left_start_up_t_uS(:) EVENTS.Left_end_up_t_uS(:)];
nm{1} = 'Lup';
SE{2} = [EVENTS.Right_start_up_t_uS(:) EVENTS.Right_end_up_t_uS(:)];
nm{2} = 'Rup';
SE{3} = [EVENTS.Left_start_down_t_uS(:) EVENTS.Left_end_down_t_uS(:)];
nm{3} = 'Ldown';
SE{4} = [EVENTS.Right_start_down_t_uS(:) EVENTS.Right_end_down_t_uS(:)];
nm{4} = 'Rdown';

% Angle of paw movement for each pull and the distance between the start and end of a pull...
% We can define bad pulls as pulls that did not travel far enough in
% distance.
clrs = lines(10);
MXY = [MT.t_uSec MT.Left_x MT.Left_y];

if PLOT_IT
    figure(201)
end

for ii = 1:length(SE)
    [OUT.(['Ang_' nm{ii}]),OUT.(['Ang_x_' nm{ii}]),OUT.(['Ang_y_' nm{ii}]),OUT.(['Ang_d_' nm{ii}])] = ...
        LK_reach_angle(MXY,SE{ii},'plot_it',PLOT_IT,'Color',clrs(ii,:));
end

th = 25;

if PLOT_IT
    for ii = 1:length(SE)
        figure(202)
        h = polarhistogram(OUT.(['Ang_' nm{ii}])(OUT.(['Ang_d_' nm{ii}])>th),30); h.FaceAlpha = .1; hold on
        h.DisplayStyle = 'stairs'; h.LineWidth = 4 ;
        figure(203)
        plot(OUT.(['Ang_x_' nm{ii}]),OUT.(['Ang_y_' nm{ii}]),'.'); hold on
        
    end
    figure(202);legend(nm); legend boxoff; figure(203);legend(nm); legend boxoff
end
% Occupancy normalized response fields. Very nice....
for ii = 1:length(SE)
    GIX = OUT.(['Ang_d_' nm{ii}])>th;
    [OUT.(['PF_' nm{ii}]),x_axis] = place_field_for_start_end_pull(TS_uS,SE{ii}(GIX,:),15,PLOT_IT);
end

new_x = linspace(0,1,res);

for iE = 1:length(EVENTS.Left_start_up_t_uS)
    % Determine things like the firing rates, velocities, and occupancies
    % for this space. (occupancy?) Time at each point would be the diff me
    % thinks.
    IX = NT_t_uS >= EVENTS.Left_start_up_t_uS(iE) & NT_t_uS < EVENTS.Left_end_up_t_uS(iE);
    % I might need more of a histogram here - just like place field
    x = linspace(0,1,sum(IX));
    if ~isempty(x)
    spk = interp1(x,NT(IX,:),new_x);
    PE.Left.up.spk(:,:,iE) = spk;
    % real time occupancy is just the mean speed at each point.
    occ = interp1(x(2:end)-.02,diff(NT_t_uS(IX)),new_x);
    PE.Left.up.occ(:,:,iE) = spk;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a reduced movement segment of just pulling times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MTp = MT(MT.Is_pulling,:);
% behavior_bin_size_ms = 1000/Behavior_sFreq;
cond = {'r','r_sh'};

OUT.r = nan(length(vbls),length(params),length(TS_uS));
OUT.r_sh = nan(length(vbls),length(params),length(TS_uS));

for iCond = 1:2
    for iNrn = 1:length(TS_uS)
        
        for iPar = 1:length(params)
            % smooth_win_ms = 200;
            smooth_win_ms = params(iPar);
            smooth_bins = round(smooth_win_ms/small_binsize_ms);
            switch cond{iCond}
                case 'r'
                    [NT,~,NT_t_uS] = Bin_and_smooth_ts_array(TS_uS, small_binsize_ms*1e3, smooth_bins);
                case 'r_sh'
                    [NT,~,NT_t_uS] = Bin_and_smooth_ts_array(TS_uS, small_binsize_ms*1e3, smooth_bins, true);
                otherwise
                    error('wtf')
            end
            NT_MT = interp1(NT_t_uS,NT,MTp.t_uSec,'nearest');
            % This is where I would split into two test sets.
            % To be significant, there must be a correlation in BOTH
            % subsets.
            pull_IDs = unique(MTp.Pull_count(MTp.Pull_count>0));
            gp1 = pull_IDs(1:2:end); gp2 = pull_IDs(2:2:end);
            IX1 = ismember(MTp.Pull_count,gp1);
            IX2 = ismember(MTp.Pull_count,gp2);
            
            x = NT_MT(:,iNrn);
            GIX = x>0; % Just select bins when the neuron actually fires.
            if sum(GIX) < 20 % need at least 20 active bins. Otherwise, forget it.
                continue
            end
            % Compute metrics for each neuron...
            % PETH - but there is no event to lock onto other than bout onset reverse correlation?, space field?
            % All of this could get corrpted if we restrcit all analysis to
            % just the pull intervals
            
            for iY = 1:length(vbls)
                y = MTp.(vbls{iY});
                if sum(~isnan(y(GIX & IX1))) < 10 || sum(~isnan(y(GIX & IX2))) < 10
                    continue
                end
                [~,pt(1)] = corr(x(GIX & IX1),y(GIX & IX1),'Type','Spearman','Rows','complete');
                [~,pt(2)] = corr(x(GIX & IX2),y(GIX & IX2),'Type','Spearman','Rows','complete');
                [r,p] = corr(x(GIX),y(GIX),'Type','Spearman','Rows','complete');
                if min(pt) > 0.05
                    % Only significant if BOTH sets show significance.
                    r = nan;
                end
                OUT.(cond{iCond})(iY,iPar,iNrn) = r;
                OUT.([cond{iCond} '_p'])(iY,iPar,iNrn) = p;
                
                if PLOT_IT && p < 0.01
                    figure(1)
                    clf
                    scatter(x,y)
                    ylabel(vbls{iY}); xlabel('nrn activity')
                    lsline
                    title(sprintf('%s nrn %d r=%2.3f p=%2.3f',vbls{iY},iNrn, r,p))
                    disp('')
                end
                
                %                 lm = fitlm(NT_MT(:,iNrn),y,'linear'); % 'RobustOpts','on');
                %                 OUT.(cond{iCond})(iY,iPar) = lm.Rsquared.Adjusted;
                %                 OUT.([cond{iCond} '_p'])(iY,iPar) = lm.coefTest;
            end
        end
    end
end

if PLOT_IT
    for iNrn = 1:size(OUT.r,3)
        figure
        imagesc(params,1:length(vbls) ,squeeze(OUT.r(:,:,iNrn)))
        set(gca,'YTickLabel',vbls)
        colorbar
    end
end



[OUT.PTHpk,OUT.PTH_x_axis] = PETH_raster_array(TS_uS,PULL.Peak_uS, 4, PETH_win_ms, PETH_win_ms);
[OUT.PTHst,OUT.PTH_x_axis] = PETH_raster_array(TS_uS,PULL.Start_uS,4, PETH_win_ms, PETH_win_ms);
% for ii = 1:length(OUT.PTHpk);figure;imagesc(OUT.PTHst{ii}>0); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike triggered averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rate Maps: spiking relative to nose.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return

M = [PAW.Time_uSec PAW.Right_x PAW.Right_y PAW.Left_x PAW.Left_y PAW.Right_speed PAW.Left_speed ];
M = M(PW.GIX,:);
for iS = 1:length(SP)
    t_sec = Restrict(SP(iS).t_uS/1e6,move_intervals_sec);
    if length(t_sec) > 40 && nanmin(p(:,iS)) < 0.01
        SpPaw = interp1(M(:,1),M(:,2:end),t_sec*1e6,'nearest');
        figure
        subplot(2,2,1)
        plot(M(:,2),M(:,3),'k.','MarkerSize',2)
        hold on
        plot(SpPaw(:,1),SpPaw(:,2),'ro','MarkerSize',4)
        title('Right')
        
        subplot(2,2,2)
        plot(M(:,4),M(:,5),'k.','MarkerSize',2)
        hold on
        plot(SpPaw(:,3),SpPaw(:,4),'ro','MarkerSize',4)
        title('Left')
        sgtitle(['nrn ' num2str(iS)])
        subplot(2,2,3)
        Plot_placefield(t_sec*1e6,M(:,1:3),[],100,11);
        subplot(2,2,4)
        Plot_placefield(t_sec*1e6,M(:,[1 4:5]),[],100,11);
    end
end

M = [IMU.t_uS IMU.data_V(:,1:2) IMU.absjerk IMU.absjerkpc ];
M = M(IMU.GIX,:);
for iS = 1:length(SP)
    t_sec = Restrict(SP(iS).t_uS/1e6,move_intervals_sec);
    if length(t_sec) > 40 && nanmin(p(:,iS)) < 0.01
        SpPaw = interp1(M(:,1),M(:,2:end),t_sec*1e6,'nearest');
        figure
        subplot(2,2,1)
        plot(M(:,2),M(:,3),'b.','MarkerSize',2)
        hold on
        plot(SpPaw(:,1),SpPaw(:,2),'ro')
        title('PC1 and 2')
        
        subplot(2,2,2)
        plot(M(:,4),M(:,5),'b.','MarkerSize',2)
        hold on
        plot(SpPaw(:,3),SpPaw(:,4),'ro')
        title('jerk and pc jerk')
        sgtitle(['nrn ' num2str(iS)])
        subplot(2,2,3)
        %         Plot_placefield(t_sec*1e6,M(:,1:3),[],100,11);
        Plot_placefield(t_sec*1e6,M(:,1:3),[],20);
        subplot(2,2,4)
        %         Plot_placefield(t_sec*1e6,M(:,[1 4:5]),[],100,11);
        Plot_placefield(t_sec*1e6,M(:,[1 4:5]),[],20);
    end
end

% Top down x y location.
M = [POS.Time_uS POS.COM_xy POS.Speed POS.COM_xy(:,1)];
M = M(POS.GIX,:);
for iS = 1:length(SP)
    t_sec = Restrict(SP(iS).t_uS/1e6,move_intervals_sec);
    if length(t_sec) > 40 && nanmin(p(:,iS)) < 0.01
        SpPaw = interp1(M(:,1),M(:,2:end),t_sec*1e6,'nearest');
        figure
        subplot(2,2,1)
        plot(M(:,2),M(:,3),'k.','MarkerSize',2)
        hold on
        plot(SpPaw(:,1),SpPaw(:,2),'ro','MarkerSize',4)
        title('POS xy')
        
        subplot(2,2,2)
        plot(M(:,4),M(:,5),'k.','MarkerSize',2)
        hold on
        plot(SpPaw(:,3),SpPaw(:,4),'ro','MarkerSize',4)
        title('speed by x')
        sgtitle(['nrn ' num2str(iS)])
        subplot(2,2,3)
        Plot_placefield(t_sec*1e6,M(:,1:3),[],100,11);
        subplot(2,2,4)
        Plot_placefield(t_sec*1e6,M(:,[1 4:5]),[],100,11);
    end
end
end % end function


function [PF,x_axis,TC,OCC] = place_field_for_start_end_pull(TS_uS,se,nbins,plot_it)
% I need to create a new definition of 'space' that spans the start and end
% of each pull.
POS = [];
for iE = 1:Rows(se)
    t = se(iE,1):1000:se(iE,2);
    d = linspace(0,1,length(t));
    POS = [POS;t(:) d(:)];
end


TS_uSr = Restrict(TS_uS,se);
% [PF,x_axis,TC,OCC] = Plot_placefield_1D(TS_uSr,POS,se,15,[],true);
[PF,x_axis,TC,OCC] = Plot_placefield_1D_by_trial(TS_uSr,POS,se,nbins,[],plot_it);

end



% string_start_end_sec = [string_pull_intervals_uSec(1) string_pull_intervals_uSec(end)]/1e6;
% move_intervals_uSec = find_intervals([ROT.t_uSec(ROT.GIX) ROT.Speed(ROT.GIX) > 4 & ROT.Speed(ROT.GIX) < 300],.5);
% pull_times_sec = ROT.Sec(diff(ROT.Speed) > 30)';
% pull_times_sec = Restrict(pull_times_sec,string_pull_intervals_uSec/1e6);
% neg_pull_times_sec = ROT.Sec(diff(ROT.Speed) < -30);
% neg_pull_times_sec = Restrict(neg_pull_times_sec,string_pull_intervals_uSec/1e6);