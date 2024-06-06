function [hb ,he]= Error_bars_Abhi(data, error_type, color_grp, color_rgb,plot_points,plot_error_bars)
%[e,m,err] =
%function Error_bars(varargin)
% Plot a bar graph of the matrices in varargin of their means and sem and stds.
% INPUT : varargin of vectors for calculating means and stds and sems on
%         assumes that all vargins only one column.
% OUTPUT: Bar graph with the SEM as a red error whisker.
%         reference to the handle of the error bars

% cowen
% changed by Abhi 2023
if nargin < 2 || isempty(error_type)
    error_type = 'Sem';
end
if nargin < 3
    color_grp = false;
end
if nargin < 4
    color_rgb = [.2 .2 .2];
end
if nargin < 5
    plot_points = true;
end
if nargin < 6
    plot_error_bars = true;
end


MarkerFaceAlpha = 0.2;

% count = 1;
% % find commands and creat an argument cell array that just has the data.
% for ii = 1:length(varargin)
%     if ischar(varargin{ii})
%         switch varargin{ii}
%             case 'std'
%                 error_type = 'std';
%             case 'nanstd'
%                 error_type = 'nanstd';
%             case '1.95*std'
%                 error_type = '1.95*std';
%             otherwise
%                 error('Incorrect error type')
%         end
%     else
%         arg{count} = varargin{ii};
%         count = count + 1;
%     end
% end
%%%%%% the real code %%%%%
if iscell(data)
    if length(data) == 1
        m   = nanmean(data{1});
        sd  = nanstd(data{1});
        eval(['err = ' error_type '(data{1);']);

    else
        for ii = 1:length(data)
            m(ii) = nanmean(data{ii});
            sd(ii) = nanstd(data{ii});
            eval(['err(ii) = ' error_type '(data{ii});']);
            err(ii) = feval(error_type,data{ii});
        end
    end   
else
    if Cols(data) == 1

        m   = nanmean(data);
        sd  = nanstd(data);
        eval(['err = ' error_type '(data);']);

    else
        for ii = 1:Cols(data)
            m(ii) = nanmean(data(:,ii));
            sd(ii) = nanstd(data(:,ii));
            eval(['err(ii) = ' error_type '(data(:,ii));']);
            err(ii) = feval(error_type,data(:,ii));
        end
    end
end


hh = bar(m);

if color_grp && length(color_rgb) > 1
    hh.FaceColor = 'Flat';
    for ii = 1:length(m)
        hh.CData(ii,:) = color_rgb(ii,:);
    end
else
    set(hh,'FaceColor',color_rgb)
end
% hh.FaceColor = 'Flat';
% hh.CData(1,:) = GP.Colors.LDOPA;
% hh.CData(1,:) = GP.Colors.LID;
% hh.CData(2,:) = GP.Colors.SHAM;
% hh.CData(3,:) = GP.Colors.SHAM;
% hh.CData(4,:) = GP.Colors.SHAM;
% set(hh,'FaceColor',GP.Colors.Baseline)
% set(hh,'FaceColor',GP.Colors.Ketamine)


hold on
if plot_points
    if iscell(data)
        for ii = 1:length(data)
            pts = data{ii};
            x = repmat(ii,length(pts),1);
            x = x + (rand(size(x))-.5)*.2;
            % plot(x,pts,'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceAlpha',MarkerFaceAlpha)
            scatter(x, pts,5,[.2 .2 .2],'filled','MarkerFaceAlpha',MarkerFaceAlpha)
        end
    else
        for ii = 1:Cols(data)
            pts = data(:,ii);
            x = repmat(ii,length(pts),1);
            x = x + (rand(size(x))-.5)*.2;
            % plot(x,pts,'o','MarkerSize',5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceAlpha',MarkerFaceAlpha)
            scatter(x, pts,10,[.2 .2 .2],'filled','MarkerFaceAlpha',MarkerFaceAlpha)
        end
    end
end

%errorbar(m,sd,'+')
if plot_error_bars
    if 1

        errorb(1:length(m),m,err,'barwidth',2); e = [];
    else

        e = errorbar(m,err,'k','LineStyle','none');
    end
else
end

hold on
% set(e,'LineWidth',2)
% set(gca,'FontSize',12)
%set(e,'MarkerSize',0.1)
%hold on
%plot(m,'+k')
%set(e(2),'MarkerSize',.1);
%set(e(2),'Color','b');
%hold off
pubify_figure_axis

if nargout > 0
    hb = hh;
    he = e;
end