function OUT = Q3_single_unit_responses_simplified()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do individual neurons correlate with features of string pulling.
% Calculate a metric for individual units. See the ensemble version for an
% ensemble-level metric.
% Cowen 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sessions with spikes 
% rat 315, 24 26 27 28 29 30 31 (crap) 32 36
%
OUT = [];
GP = LK_Globals;
PLOT_IT = true;
Behavior_sFreq = 100;
neuron_quality_threshold = 2;
small_binsize_ms = 5;
PETH_win_ms = 1500; bin_and_inc_size_msec = 20;

params = 50; % smooth
close all
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
%
% FCC = load('Filtered_Time_Stamped_Coordinates_Corrected_Ori');
% GSI = load('good_string_pull_intervals_uSec.mat');

% Estimate the individual pull times.
% PULL = LK_Determine_Individual_Pull_Times();
% Do PETHs on these data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and create a magnificent table of motion data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MT, string_pull_intervals_uSec,~,~,~,~,EVENTS] = LK_Combine_All_String_Pull_Motion_To_Table(Behavior_sFreq, false);
epoch_st_ed_uSec = [string_pull_intervals_uSec(1) - 2e6 string_pull_intervals_uSec(end) + 2e6 ]; % add a little padding.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load spikes and restrict to the times of interest.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SP, TS_uS] = LK_Load_Spikes(neuron_quality_threshold,epoch_st_ed_uSec);
OUT.n_neurons = length(SP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the spikes on top of each trajectory.
%%
i = 1;
COND(i).sted = [EVENTS.Left_start_up_t_uS EVENTS.Left_end_up_t_uS]; 
COND(i).paw = 'Left';
COND(i).dir = 'up';
COND(i).titstr = [COND(i).paw ' ' COND(i).dir];
i = i + 1;

COND(i).sted = [EVENTS.Right_start_up_t_uS EVENTS.Right_end_up_t_uS]; 
COND(i).paw = 'Right';
COND(i).dir = 'up';
COND(i).titstr = [COND(i).paw ' ' COND(i).dir];
i = i + 1;

COND(i).sted = [EVENTS.Left_start_down_t_uS EVENTS.Left_end_down_t_uS]; 
COND(i).paw = 'Left';
COND(i).dir = 'down';
COND(i).titstr = [COND(i).paw ' ' COND(i).dir];
i = i + 1;

COND(i).sted = [EVENTS.Right_start_down_t_uS EVENTS.Right_end_down_t_uS]; 
COND(i).paw = 'Right';
COND(i).dir = 'down';
COND(i).titstr = [COND(i).paw ' ' COND(i).dir];


for iN = 1:length(SP)
    figure(iN)
    h = [];
    found_spikes = false;

    for iC = 1:length(COND)
        xstr = [COND(iC).paw '_x'];
        ystr = [COND(iC).paw '_y'];
        h(iC) = subplot(2,2,iC);
        % First plot the trajectories
        for ii = 1:Rows(COND(iC).sted)
            IX = MT.t_uSec > COND(iC).sted(ii,1) & MT.t_uSec < COND(iC).sted(ii,2);
            if sum(IX) > 10 && ~isnan(sum(MT.(xstr)(IX)))
                % Find the positions that most closely match these spikes.
                mnx = mean(MT.(xstr)(IX));
                mny = mean(MT.(ystr)(IX));

                %                 mnx = 0; mny = 0;

                plot(MT.(xstr)(IX)-mnx,MT.(ystr)(IX)-mny,'-','Color',[.88 .88 .88])
                hold on
            end
        end
        % now plot the spikes
        for ii = 1:Rows(COND(iC).sted)
            t_uS = Restrict(SP(iN).t_uS,COND(iC).sted(ii,:));
            IX = MT.t_uSec > COND(iC).sted(ii,1) & MT.t_uSec < COND(iC).sted(ii,2);
            
            if sum(IX) > 10 && length(t_uS) > 2 && ~isnan(sum(MT.(xstr)(IX)))
                found_spikes = true;
                % Find the positions that most closely match these spikes.
                sx = interp1(MT.t_uSec(IX),MT.(xstr)(IX), t_uS, 'spline');
                sy = interp1(MT.t_uSec(IX),MT.(ystr)(IX), t_uS, 'spline');
                mnx = mean(MT.(xstr)(IX));
                mny = mean(MT.(ystr)(IX));
                %                 mnx = 0; mny = 0;

                plot(sx-mnx,sy-mny,'.','Color',[.0 .0 .0],'MarkerSize',8)
            end
        end
        axis tight
        title(sprintf('%s', COND(iC).titstr),'Interpreter','none')
        xlabel('x');ylabel('y')
    end
    sgtitle({[ SP(iN).fname ' Neuron ' num2str(iN)] pwd},'FontSize',8,'Interpreter','none')
    if ~found_spikes
        close
    else 
        %         equalize_axes(h);
    end
end
% Look to see if neurons respond to start of bout.
for iN = 1:length(SP)

    figure(iN+100)
    PETH_raster_array(SP(iN).t_uS,string_pull_intervals_uSec(:,1), bin_and_inc_size_msec, PETH_win_ms, PETH_win_ms);
    sgtitle({[ SP(iN).fname ' Neuron ' num2str(iN)] pwd},'FontSize',8,'Interpreter','none')

end
%%
disp(['done'])

% Do the neurons map onto the start or end of a bout.

