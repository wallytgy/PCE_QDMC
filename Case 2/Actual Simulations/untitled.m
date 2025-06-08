% ------------ choose one cluster (e.g. cluster #1) ---------------
clusterID   = 1;                                % <‑‑ change as needed
nodes       = representatives(clusterID,:);     % two node numbers
[~,idx]     = ismember(nodes, flat_reps.');

% quick sanity
fprintf('Rows used for cluster %d: %d and %d\n', clusterID, idx(1), idx(2));

% ------------ plotting --------------------------------------------
Nt = length(tsp_sim)-1;  t = tsp_sim(1:Nt).';

figure; ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on');

for k = 1:size(u_store,3)
    plot(ax, t, u_store(idx(1),1:Nt,k), 'Color',[0 0.45 0.74 0.2]);
    plot(ax, t, u_store(idx(2),1:Nt,k), 'Color',[0.85 0.33 0.10 0.2]);
end

title(ax, sprintf('Cluster %d – overlay diagnostic',clusterID));
xlabel(ax,'Time [s]'); ylabel(ax,'u [pu]');

% ------------ how many lines ended up on the axis? -----------------
nLines = numel(ax.Children);
fprintf('Axis now contains %d graphic objects (should be %d).\n',...
        nLines, 2*size(u_store,3));