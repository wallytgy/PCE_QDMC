% ------------ choose one cluster (e.g. cluster #1) ---------------
clusterID   = 1;                                % <‑‑ change as needed
nodes       = representatives(clusterID,:);     % two node numbers
[~,idx]     = ismember(nodes, flat_reps.');

% quick sanity
fprintf('Rows used for cluster %d: %d and %d\n', clusterID, idx(1), idx(2));

% ------------ plotting --------------------------------------------
Nt = length(tsp_sim)-1;  t = tsp_sim(1:Nt).';

figure; 
hold on
for k = 1:size(u_store_all,3)
    plot(t, u_store_all(idx(1),1:Nt,k),'r');
    plot(t, u_store_all(idx(2),1:Nt,k)),'b';

end
title(sprintf('Cluster %d – overlay diagnostic',clusterID));
xlabel('Time [s]'); ylabel('u [pu]');

% ------------ how many lines ended up on the axis? -----------------
nLines = numel(ax.Children);
fprintf('Axis now contains %d graphic objects (should be %d).\n',...
        nLines, 2*size(u_store,3));