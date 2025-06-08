
function hFig = plotClusterTrajectories(u_store_all, flat_reps, representatives, tsp_sim,filename, varargin)
%PLOTCLUSTERTRAJECTORIES  Visualise control inputs cluster by cluster.
%
%   hFig = plotClusterTrajectories(u_store, flat_reps, reps, tsp_sim, Name,Value)
%
%   INPUTS
%     u_store         [n_in × Nt × n_traj]  control trajectories
%     flat_reps       [n_in × 1]            node number for each u_store row
%     representatives [6 × 2]               node numbers per cluster
%     tsp_sim         [1 × Nt+1]            time stamps
%
%   OPTIONAL NAME–VALUE PAIRS
%     'PlotStyle' — 'lines' | 'envelope'     (default 'lines')
%     'Opacity'   — scalar in [0,1]          (default 0.15)
%     'ShowMean'  — logical                  (default true)
%     'Colors'    — 2×3 RGB matrix           (default blue/red)
%
%   OUTPUT
%     hFig — handle to created figure
%
%   Example (transparent lines):
%     plotClusterTrajectories(u_store, flat_reps, reps, tsp_sim,...
%                             'Opacity',0.1,'PlotStyle','lines');
%
%   Example (min–max envelope):
%     plotClusterTrajectories(u_store, flat_reps, reps, tsp_sim,...
%                             'PlotStyle','envelope');

% --------------------------------------------------------------------
p = inputParser;
addParameter(p,'PlotStyle','lines',@(s) any(strcmpi(s,{'lines','envelope'})));
addParameter(p,'Opacity',0.15,@(x) isnumeric(x)&&x>=0&&x<=1);
addParameter(p,'ShowMean',true,@islogical);
addParameter(p,'Colors',[0.11 0.53 0.97;   % node‑1 colour
                          0.88 0.28 0.25]); % node‑2 colour
parse(p,varargin{:});

plotStyle = lower(p.Results.PlotStyle);
alphaVal  = p.Results.Opacity;
showMean  = p.Results.ShowMean;
clr       = p.Results.Colors;

% sizes
[n_in, Nt, n_traj] = size(u_store_all);
assert(size(representatives,1)==6,  'Expected 6 clusters (3×2 grid).');
assert(size(representatives,2)==2,  'Each cluster must contain two nodes.');

t = tsp_sim(1:Nt).';                 % column vector for patch use
hFig = figure('Color','w','Name','Cluster‑wise control trajectories');

% helper for MATLAB < R2018a (no alpha support)
supportsAlpha = verLessThan('matlab','9.4')==0;

for c = 1:6
    % locate the two actuated rows for this cluster
    [~,idx] = ismember(representatives(c,:), flat_reps.');
    if any(idx==0)
        warning('Cluster %d: could not find both node rows in u_store.',c);
        continue;
    end

    % extract the trajectories: Nt × n_traj
    U1 = squeeze(u_store_all(idx(1),1:Nt,:));
    U2 = squeeze(u_store_all(idx(2),1:Nt,:));

    subplot(3,2,c);  hold on;  box on;  grid on;  set(gca,'Layer','top');

    switch plotStyle
        case 'lines'
            for k = 1:n_traj
                h1 = plot(t, U1(:,k), 'Color', clr(1,:), 'LineWidth',0.8);
                h2 = plot(t, U2(:,k), 'Color', clr(2,:), 'LineWidth',0.8);
                if supportsAlpha          % add transparency if available
                    h1.Color(4) = alphaVal;
                    h2.Color(4) = alphaVal;
                else                      % fallback: lighten colours
                    h1.Color = clr(1,:)*(1-alphaVal) + alphaVal;
                    h2.Color = clr(2,:)*(1-alphaVal) + alphaVal;
                end
            end

        case 'envelope'
            % choose envelope type: here min–max; change to ±std if desired
            lo1 = min(U1,[],2);  hi1 = max(U1,[],2);
            lo2 = min(U2,[],2);  hi2 = max(U2,[],2);

            % patches (translucent)
            p1 = patch([t; flipud(t)], [lo1; flipud(hi1)], clr(1,:), ...
                       'EdgeColor','none');
            p2 = patch([t; flipud(t)], [lo2; flipud(hi2)], clr(2,:), ...
                       'EdgeColor','none');
            if supportsAlpha
                p1.FaceAlpha = alphaVal;
                p2.FaceAlpha = alphaVal;
            else
                p1.FaceColor = clr(1,:)*(1-alphaVal) + alphaVal;
                p2.FaceColor = clr(2,:)*(1-alphaVal) + alphaVal;
            end
    end

    % mean lines (optional; clearer with transparency)
    if showMean
        hMean1 = plot(t, mean(U1,2), 'Color', clr(1,:), 'LineWidth',1.6);
        hMean2 = plot(t, mean(U2,2), 'Color', clr(2,:), 'LineWidth',1.6);
    end

    title(sprintf('Cluster %d [K]',c),'FontWeight','normal');
    xlabel('Time [s]');  ylabel('u [pu]');
    ylim tight;
    % if c==1
    %     hRuns  = [h1 h2];           % first trajectory of each node is fine
    %     hMeans = [hMean1 hMean2];   % handles returned when you plot the means
    %     if plotStyle=="lines"
    %         if showMean
    %             lg = legend([hRuns hMeans], ...
    %                         {'Node 1 (runs)','Node 2 (runs)', ...
    %                          'Node 1 (mean)','Node 2 (mean)'}, ...
    %                         'Location','best','FontSize',8);
    %         else
    %             lg = legend({'Node 1 (runs)','Node 2 (runs)'});
    %         end
    %     else  % envelope
    %         if showMean
    %             lg = legend({'Node 1 (envelope)','Node 2 (envelope)',...
    %                          'Node 1 (mean)','Node 2 (mean)'});
    %         else
    %             lg = legend({'Node 1 (envelope)','Node 2 (envelope)'});
    %         end
    %     end
    %     set(lg,'Location','best','FontSize',8);
    %     % Thicken all lines appearing in the legend
    %     set(hRuns,'LineWidth',10); 
    % end
end
print(gcf, filename, '-depsc2', '-painters');
end