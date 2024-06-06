function O = Ensemble_analysis_and_comparison(Q, metrics, varargin)
% Compare ensemble recordings between two epochs in terms of ensemble
% metrics and generate nice plots.
% Q is a CELL ARRAY of TWO Q matrices.
% For now assumes you ONLY pass in 2 Q matrices.
% It is up to you to make the Q (ncells (column) ,time (row)) matrices.
% You determine the bin size.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLOT_IT = true;
x_label = [];

Extract_varargin

if nargin < 2
    metrics = {'corr' 'sparsity' 'n_effective_dimensions' 'state_similarity'}; % show everything;
end

O = [];
R = []; SparseQ = [];lat_explained = [];

for iM = 1:length(metrics)
    switch metrics{iM}
        case 'corr'
            % Assumes simultaneously recorded cells.
            for iQ = 1:length(Q)
                R{iQ} = Corr_upper_only(Q{iQ});
            end

            if PLOT_IT
                plot_R(Q,R,x_label);
            end
            O.R = R;

        case 'sparsity'
            % Can be applied to non-simultaneously recorded cells
            for iQ = 1:length(Q)
                QQ = Q{iQ} + abs(min(Q{iQ}(:))) + 1; % For some measures of sparsity, you need to have only positive values.
                SparseQ{iQ} = Sparseness_measures(QQ);
            end
            O.SparseQ = SparseQ;

            if PLOT_IT
                plot_Sparse(Q,SparseQ,x_label);
            end

        case 'n_effective_dimensions'
            % Assumes simultaneously recorded cells.
            lat_explained = []; neffdim = [];
            for iQ = 1:length(Q)
                warning off
                [pc,sc,lat,~,lat_explained(:,iQ)] = pca(Q{iQ},'Economy',false);
                neffdim_abbott(iQ) = n_effective_dimensions(Q{iQ});
                warning on
            end
            O.nEffDim.lat_explained = lat_explained;
            O.nEffDim.neffdim_abbott = neffdim_abbott;
            if PLOT_IT
                plot_n_effective_dimension(Q,O.nEffDim,x_label);
            end

        case 'state_similarity'
            % How similar are these two states?

            for iQ = 1:length(Q)
                R{iQ} = Corr_upper_only(Q{iQ});
            end

            [O.StatSim.r_sim,O.StatSim.p_sim ] = corr(R{1}.r,R{2}.r,'Rows','complete');
            O.StatSim.R = R;
            if PLOT_IT
                plot_state_similarity(Q,O.StatSim,x_label);
            end

        otherwise
            error('Unknown metric')
    end

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_R(Q,R,x_label)
figure
for iQ = 1:length(Q)

    subplot(4,length(Q),iQ);
    imagesc(x_label,[], Q{iQ}')
    colormap(turbo);
    colorbar('eastoutside')
    pubify_figure_axis
    ylabel('Neuron')
    if iQ == 1
        cl = clim;
    end
    clim(cl);

    subplot(4,length(Q),iQ + length(Q));
    imagesc(R{iQ}.Rfull)
    axis square
    xlabel('Neuron');ylabel('Neuron')
    pubify_figure_axis
    if iQ == 1
        cl2 = clim;
    end
    clim(cl2);

    subplot(4,length(Q),iQ + 2*length(Q));
    histogram(R{iQ}.r,80)
    xlabel('r');ylabel('Count')
    pubify_figure_axis
    if iQ == 1
        xl = xlim;
    end
    xlim(xl);

end
subplot(4,length(Q),4*length(Q)-1:4*length(Q));
histogram(R{2}.r-R{1}.r,80)
xlabel('r2-r1')
pubify_figure_axis
plot_vert_line_at_zero
set(gcf,"Position",[ 469   106   825   844])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Sparse(Q,Sparse,x_label)
for iQ = 1:length(Q)
    mnmx(iQ,:) = [min(Sparse{iQ}.CV) max(Sparse{iQ}.CV)];
end

yl = [min(mnmx(:,1)) max(mnmx(:,2))];

figure
for iQ = 1:length(Q)

    subplot(1,length(Q),iQ);
    imagesc(x_label,[], Q{iQ}')
    colormap(turbo);
    colorbar('eastoutside')
    pubify_figure_axis
    ylabel('Neuron')
    if iQ == 1
        cl = clim;
    end
    clim(cl);
    yyaxis right
    if isempty(x_label)
        plot(Sparse{iQ}.CV,'k','LineWidth',3)
    else
        plot(x_label,Sparse{iQ}.CV,'k','LineWidth',3)
    end
    % if iQ == 1
    % yl(i) = ylim;
    % end
    ylim(yl)

end
set(gcf,'Position',[282         395        1220         420])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_n_effective_dimension(Q,nEffDim,x_label)
figure
subplot(1,2,1)
bar(nEffDim.neffdim_abbott)
pubify_figure_axis
ylabel('N Effective Dim (Abbott)')
subplot(1,2,2)
plot(nEffDim.lat_explained)
ylabel('Latent explained variance')
pubify_figure_axis
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_state_similarity(Q,StateSim,x_label)
SIGIX = StateSim.R{1}.p < 0.25 | StateSim.R{2}.p < 0.25;
figure
subplot(1,2,1)
bar(StateSim.r_sim)
pubify_figure_axis
ylabel('r')
subplot(1,2,2)
scatter(StateSim.R{1}.r, StateSim.R{2}.r,3)
lsline
hold on
scatter(StateSim.R{1}.r(SIGIX), StateSim.R{2}.r(SIGIX),15,'r+')

lsline
ylabel('correlations in Q1 and Q2. red p<0.25 only')
xlabel('r Q1'); ylabel('r Q2')
pubify_figure_axis
end