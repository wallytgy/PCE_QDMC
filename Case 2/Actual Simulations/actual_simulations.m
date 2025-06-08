clear all
close all
%% define simulation options
rng(42);

case118data = pglib_opf_case118_ieee
branch = case118data.branch
edges = branch(:,1:2)
xij = branch(:,4)
G = graph(edges(:,1),edges(:,2),xij)


%% get clusters

% Step 1: Get weighted Laplacian
Gsimp = simplify(G)
Adj = adjacency(Gsimp, 'weighted');  % Weighted adjacency matrix
D = diag(sum(Adj, 2));           % Degree matrix
L = D - Adj;    

% Step 2: Choose number of clusters
k = 8;  % adjust as needed

% Step 3: Compute first k eigenvectors of Laplacian
[V, ~] = eigs(L, k, 'smallestabs'); 

% Step 4: Cluster rows of eigenvector matrix
idx = kmeans(V, k); 


min_size   = 3;                       % threshold
cluster_sz = accumarray(idx,1);       % size of each cluster
k          = numel(cluster_sz);   


cluster_sizes = accumarray(idx, 1);  % how many nodes per cluster


% Copy the original labels
min_size = 3;  % Threshold for what counts as 'too small'

new_idx = idx;                        % start with the original labels

% pre‑compute the “largest” cluster label for the fallback case
[~, largestC] = max(cluster_sz);

for c = 1:k
    if cluster_sz(c) <= min_size                 % --- tiny cluster ---
        nodes_c = find(idx == c);                % nodes in that cluster

        for v = nodes_c'                         % iterate over those nodes
            nbrs      = neighbors(G, v);         % <-- FUNCTION call
            nbrLabels = idx(nbrs);               % labels of neighbours
            nbrLabels(nbrLabels == c) = [];      % throw away tiny‑cluster label

            if ~isempty(nbrLabels)
                new_idx(v) = mode(nbrLabels);    % majority vote
            else
                new_idx(v) = largestC;           % fallback
            end
        end
    end
end

figure;
[~,~,new_idx] = unique(new_idx);
idx = new_idx 
k = max(idx)
deg = degree(G);  % degree of all nodes

k = max(new_idx);
representatives = zeros(k, 2);  % two per cluster

deg = degree(G);  % degree of all nodes

for c = 1:k
    nodes = find(new_idx == c);
    
    % Get degrees of these nodes
    [~, sorted_idx] = sort(deg(nodes), 'descend');
    
    % Select top 2 (or just all if fewer than 2 nodes)
    top2 = nodes(sorted_idx(1:min(2, end)));
    reps = zeros(1, 2); reps(1:numel(top2)) = top2;
    representatives(c, :) = reps;
    
    fprintf('Cluster %d → Nodes %s\n', c, mat2str(reps(reps > 0)));
end
flat_reps = nonzeros(representatives(:));
p = plot(G, 'NodeCData', new_idx, 'Layout', 'force', 'MarkerSize', 5);
highlight(p, flat_reps, 'Marker', 'o', 'MarkerSize', 9);
title('IEEE 118-bus Network');

print(gcf, 'graph_network', '-depsc2', '-vector');



%%
%%%%%%%%%%%%%%%%%% OPEN LOOP SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
mc_samples = 1e2;
pce_order = 4 ;
delta = 0.1;
n_t = 1001;
n        = max(edges(:));        % number of nodes
breakpts = [0 1] * 5;
tspan = breakpts(1):(breakpts(end)-breakpts(1))/(n_t-1):breakpts(end);
Uvals = zeros(n,1)
Uvals(flat_reps) = 1
u_fun = @(t) piecewise_u(t, breakpts, Uvals);
x_responses_MC = zeros(mc_samples, n);
Bij = 1./xij;
Bij = Bij / mean(abs(Bij));
Vref = 1;
M = 1;
eta1 = 1
eta2 = 1
[A_nom,B_nom]  = build_AB(edges, Bij, Vref, M, eta2);

%% MC open loop for various values of eta1 and Bij
for i = 1:mc_samples
    i
    theta_1_sample = unifrnd(-1,1);
    theta_2_sample = unifrnd(-1,1);
    A = (1 + delta * 3^(1/2) * theta_1_sample) * A_nom - eta1 * eye(n);
    B = (1 + delta * 3^(1/2) * theta_2_sample) * B_nom;
    ode    = @(t,x)  A*x + B * u_fun(t);
    [t,x]  = ode15s(ode, tspan, zeros(n,1));
    x_responses_MC(i,:) = x(end,:);
end
%% PCE open loop for various values of eta1 and Bij
syms theta1 theta2
PCE_polynomials = generate_polynomials(pce_order, theta1, theta2);
N = length(PCE_polynomials)
%%%% construct \mathcal{A}_{PCE}
S = zeros(N);
for i = 1:N
    for j = 1:N
        poly = PCE_polynomials(i) * PCE_polynomials(j);
        S(i,j) = int(int(3^(1/2)*theta1 * poly * 1/2  * 1/2, theta1, [-1,1]), theta2, [-1,1]);
    end
end

A_PCE = kron(eye(N) + delta * S, A_nom) - eta1 * kron(eye(N), eye(n)) ; 
%%%% construct \mathcal{B}_PCE
v = zeros(N,1)
for i = 1:N
    v(i) = int(int((1 + delta * 3^(1/2) * theta2) * PCE_polynomials(i) * 1/2  * 1/2, theta1, [-1,1]), theta2, [-1,1]);
end

B_PCE = kron(v, B_nom)
ode_PCE    = @(t,X) A_PCE*X + B_PCE * u_fun(t);
[t,x_PCE]  = ode15s(ode_PCE, tspan, zeros(N * n,1));
x_PCE = reshape(x_PCE, [n_t, n, N]);

n_t = 1001

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CLOSED LOOP SIMULATION %%%%%%%%%%%%%%%%%%%%%%
time_horizon = 20 
breakpts = [0 1] * time_horizon;
tspan = breakpts(1):(breakpts(end)-breakpts(1))/(n_t-1):breakpts(end);
step_respones = zeros(length(flat_reps), n_t, k , N);
for i = 1:length(flat_reps)
    Uvals = zeros(n,1)
    Uvals(flat_reps(i)) = 1
    u_fun = @(t) piecewise_u(t, breakpts, Uvals);
    ode_PCE    = @(t,X) A_PCE*X + B_PCE * u_fun(t);
    [t,x_PCE_responses]  = ode15s(ode_PCE, tspan, zeros(N* n,1));
    x_PCE_responses = reshape(x_PCE_responses, [n_t, n, N]);
    for c = 1:k
        % coeff‑wise cluster averages
        nodes = find(idx==c);
        output_c  = mean(x_PCE_responses(:,nodes,:), 2);            % n_t ×1   (mean part)
        step_respones(i,:,c,:) = output_c;
        
    end
end
[n1,n2,n3,n4] = size(step_respones);          % grab the sizes (optional)
y_responses = reshape(step_respones, n1, n2, n3*n4); 
y_responses = permute(y_responses, [3 1 2]);


%% Parameters to create step response
dt = 1; 
n_sr = 20 % Delta t (has to be identical when creating the step)
t_simulation = 400 ; % response model and when running closed loop control)
tsp_sr = 0:time_horizon/10:time_horizon;   % Time span to create the step response % Length of the step response model
ops = "ode45";   % ODE options
idx_sr = (1:n_sr) * (n_t-1) / n_sr + 1;
Su = y_responses(:,:,idx_sr); 

%%
n_in = length(flat_reps); 
n_out = k * N;
n_dis = 0;
Sd = [];
p = 20;        % Prediction horizon (number of time steps)
c = 5;         % Control horizon (number of time steps)
La =ones(n_in,1);         % Weight for input movement
Q = ones(n_out,1);
Q(1:k) = 10;
ctg = [];       % Cost to go

u_past = zeros(n_in,1);    % Past input (assume process starts at steady state)
y_past = zeros(n_out,1);     % Past output (assume process starts at steady state)
d_past = [];    % Past measured disturbance
u_int = [];     % Integrating inputs
d_int = [];     % Integrating measured disturbances


u_min = -100*ones(n_in,1);     % Minimum input
u_max = 100*ones(n_in,1); % Maximum input
D_u_min = -25*ones(n_in,1);  % Minimum input movement
D_u_max = 25*ones(n_in,1);   % Maximum input movement
u_dist = 0.5 * [zeros(1,50) ones(1,50) zeros(1,100) -1*ones(1,100) zeros(1,100)];
y_min = -2 *ones(n_out,1);   % Minimum output
y_max = 2*ones(n_out,1);    % Maximum output


soft = [];   % Both, y_max and y_min are soft constraints
w_eps = [];     % Weight for linear constraint violation
W_eps = [];     % Weight for quadratic constraint violation

n_x = 118;

tsp_sim = 0:dt:t_simulation;
n_traj = 100;
y_store_all = zeros(6, 401, n_traj);
u_store_all = zeros(n_in,length(tsp_sim)-1,n_traj);

%%
rng(42)
for iteration = 1:n_traj
iteration
theta_1_sample = unifrnd(-1,1);
theta_2_sample = unifrnd(-1,1);
% randomly generate setpoints
ysp = zeros(k* N,length(tsp_sim)-1);
N_interval = (length(tsp_sim)-1)/4;
ysp(1:k,1:N_interval) = (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,N_interval+1:2*N_interval) = - (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,2*N_interval+1:3*N_interval) = (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,3*N_interval+1:4*N_interval) = - (-1 + 2 * ones(k,1)) * ones(1,N_interval);
udist = [];


% Create instance of the controller
Controller = QDMC_controller_soft(n_in,n_out,n_dis,Su,Sd,p,c,La,Q,ctg,...
                               u_min,u_max,D_u_min,D_u_max,y_min,y_max,...
                               u_past,y_past,d_past,...
                               u_int,d_int,soft,w_eps,W_eps);

y_actual_store = zeros(k,length(tsp_sim));
y_store = zeros(n_out,length(tsp_sim));     % Outputs are measured at t
u_store = zeros(n_in,length(tsp_sim)-1);
y = zeros(n_out,1); % Start the process from steady state (deviation variable)
% Write value for t=0 to storage
x_curr = zeros(n,1);

    %%% run closed loop
    for ii = 1:length(tsp_sim)-1
        % Current setpoint
        y = y_store(:,ii);
        ysp_curr = ysp(:,ii);
        % Current measured disturbance
        d = [];
        
        % Calculate controller output = process input
        u = Controller.run_controller(y,ysp_curr,d);
        u_proj = projection(u,flat_reps,n);
        breakpts = [0 dt];
        Uvals    = [u_proj];
        u_fun = @(t) piecewise_u(t, breakpts, Uvals);
        % Simulate one time step
        ode    = @(t,x)  (1 + delta * 3^(1/2) * theta_1_sample) * A*x + (1 + delta * 3^(1/2) * theta_2_sample) *B * u_fun(t);
        [t,x]  = ode45(ode, [0 dt], x_curr);
        x_curr = x(end,:);
        W = sparse(1:numel(idx), idx, 1, numel(idx), k);
        sums   = x_curr * W;                % n_t × K   sum of temps in each cluster
        counts = full(sum(W,1));       % 1   × K   vertices per cluster (all > 0)
        y_actual  = sums ./ counts;       % broadcast divide → averages
        y_store(1:k,ii+1) = y_actual;
        y_actual_store(:,ii+1) = y_actual;
        u_store(:,ii) = u;
    end
    y_store_all(:,:,iteration) = y_actual_store;
    u_store_all(:,:,iteration) = u_store;
end
%%
plot_trajectories(y_store_all, ysp,"varying_sp_nd_case_2_output.eps","r",false)
y_store_all_nd = y_store_all
%%
% 2) Shaded min–max envelopes plus mean trajectories
plotClusterTrajectories(u_store_all, flat_reps, representatives, tsp_sim, "varying_sp_nd_case_2_control_zerovar.eps",...
                        'PlotStyle','lines','Opacity',0.1,'ShowMean',true);



%%
rng(42)
u_dist = 10 * [zeros(n_in,50) ones(n_in,100) zeros(n_in,100) -1*ones(n_in,100) zeros(n_in,50)];
for iteration = 1:n_traj
iteration
theta_1_sample = unifrnd(-1,1);
theta_2_sample = unifrnd(-1,1);
% randomly generate setpoints
ysp = zeros(k* N,length(tsp_sim)-1);
N_interval = (length(tsp_sim)-1)/4;
ysp(1:k,1:N_interval) = (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,N_interval+1:2*N_interval) = - (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,2*N_interval+1:3*N_interval) = (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,3*N_interval+1:4*N_interval) = - (-1 + 2 * ones(k,1)) * ones(1,N_interval);
udist = [];


% Create instance of the controller
Controller = QDMC_controller_soft(n_in,n_out,n_dis,Su,Sd,p,c,La,Q,ctg,...
                               u_min,u_max,D_u_min,D_u_max,y_min,y_max,...
                               u_past,y_past,d_past,...
                               u_int,d_int,soft,w_eps,W_eps);

y_actual_store = zeros(k,length(tsp_sim));
y_store = zeros(n_out,length(tsp_sim));     % Outputs are measured at t
u_store = zeros(n_in,length(tsp_sim)-1);
y = zeros(n_out,1); % Start the process from steady state (deviation variable)
% Write value for t=0 to storage
x_curr = zeros(n,1);

    %%% run closed loop
    for ii = 1:length(tsp_sim)-1
        % Current setpoint
        y = y_store(:,ii);
        ysp_curr = ysp(:,ii);
        % Current measured disturbance
        d = [];
        
        % Calculate controller output = process input
        u = Controller.run_controller(y,ysp_curr,d);
        u = u + u_dist(:,ii);
        u_proj = projection(u,flat_reps,n);
        breakpts = [0 dt];
        Uvals    = [u_proj];
        u_fun = @(t) piecewise_u(t, breakpts, Uvals);
        % Simulate one time step
        ode    = @(t,x)  (1 + delta * 3^(1/2) * theta_1_sample) * A*x + (1 + delta * 3^(1/2) * theta_2_sample) *B * u_fun(t);
        [t,x]  = ode45(ode, [0 dt], x_curr);
        x_curr = x(end,:);
        W = sparse(1:numel(idx), idx, 1, numel(idx), k);
        sums   = x_curr * W;                % n_t × K   sum of temps in each cluster
        counts = full(sum(W,1));       % 1   × K   vertices per cluster (all > 0)
        y_actual  = sums ./ counts;       % broadcast divide → averages
        y_store(1:k,ii+1) = y_actual;
        y_actual_store(:,ii+1) = y_actual;
        u_store(:,ii) = u;
    end
    y_store_all(:,:,iteration) = y_actual_store;
    u_store_all(:,:,iteration) = u_store;
end
%%
plot_trajectories(y_store_all, ysp,"varying_sp_id_case_2_output.eps","g"  ,true)
y_store_all_id = y_store_all
%%
% 2) Shaded min–max envelopes plus mean trajectories
plotClusterTrajectories(u_store_all, flat_reps, representatives, tsp_sim, "varying_sp_id_case_2_control_zerovar.eps",...
                        'PlotStyle','lines','Opacity',0.1,'ShowMean',true);





%%
y_dist = 0.5 * [zeros(k,50) ones(k,100) zeros(k,100) -1*ones(k,100) zeros(k,50)];
rng(42)
for iteration = 1:n_traj
iteration
theta_1_sample = unifrnd(-1,1);
theta_2_sample = unifrnd(-1,1);
% randomly generate setpoints
ysp = zeros(k* N,length(tsp_sim)-1);
N_interval = (length(tsp_sim)-1)/4;
ysp(1:k,1:N_interval) = (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,N_interval+1:2*N_interval) = - (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,2*N_interval+1:3*N_interval) = (-1 + 2 * ones(k,1)) * ones(1,N_interval);
ysp(1:k,3*N_interval+1:4*N_interval) = - (-1 + 2 * ones(k,1)) * ones(1,N_interval);
udist = [];


% Create instance of the controller
Controller = QDMC_controller_soft(n_in,n_out,n_dis,Su,Sd,p,c,La,Q,ctg,...
                               u_min,u_max,D_u_min,D_u_max,y_min,y_max,...
                               u_past,y_past,d_past,...
                               u_int,d_int,soft,w_eps,W_eps);

y_actual_store = zeros(k,length(tsp_sim));
y_store = zeros(n_out,length(tsp_sim));     % Outputs are measured at t
u_store = zeros(n_in,length(tsp_sim)-1);
y = zeros(n_out,1); % Start the process from steady state (deviation variable)
% Write value for t=0 to storage
x_curr = zeros(n,1);

    %%% run closed loop
    for ii = 1:length(tsp_sim)-1
        % Current setpoint
        y = y_store(:,ii);
        ysp_curr = ysp(:,ii);
        % Current measured disturbance
        d = [];
        
        % Calculate controller output = process input
        u = Controller.run_controller(y,ysp_curr,d);
        u_proj = projection(u,flat_reps,n);
        breakpts = [0 dt];
        Uvals    = [u_proj];
        u_fun = @(t) piecewise_u(t, breakpts, Uvals);
        % Simulate one time step
        ode    = @(t,x)  (1 + delta * 3^(1/2) * theta_1_sample) * A*x + (1 + delta * 3^(1/2) * theta_2_sample) *B * u_fun(t);
        [t,x]  = ode45(ode, [0 dt], x_curr);
        x_curr = x(end,:);
        W = sparse(1:numel(idx), idx, 1, numel(idx), k);
        sums   = x_curr * W;                % n_t × K   sum of temps in each cluster
        counts = full(sum(W,1));       % 1   × K   vertices per cluster (all > 0)
        y_actual  = sums ./ counts;
        y_actual = y_actual';
        y_actual = y_actual + y_dist(:,ii);
        y_store(1:k,ii+1) = y_actual;
        y_actual_store(:,ii+1) = y_actual;
        u_store(:,ii) = u;
    end
    y_store_all(:,:,iteration) = y_actual_store;
    u_store_all(:,:,iteration) = u_store;
end
%%
plot_trajectories(y_store_all, ysp ,"varying_sp_od_case_2_output.eps", "magenta", true)
y_store_all_od = y_store_all
% 2) Shaded min–max envelopes plus mean trajectories
plotClusterTrajectories(u_store_all, flat_reps, representatives, tsp_sim,"varying_sp_od_case_2_control_zerovar.eps",...
                        'PlotStyle','lines','Opacity',0.1,'ShowMean',true);



%%

%%%%%%%%%%%% compute all IAE and width


IAE_nd = computeAvgIAE(y_store_all_nd(:,2:end,:), ysp(1:k,:), dt)
width_nd = computeEnvelopeWidth(y_store_all_nd)
IAE_id = computeAvgIAE(y_store_all_id(:,2:end,:), ysp(1:k,:), dt)
width_id = computeEnvelopeWidth(y_store_all_id)
IAE_od = computeAvgIAE(y_store_all_od(:,2:end,:), ysp(1:k,:), dt)
width_od =computeEnvelopeWidth(y_store_all_od)




























function u = piecewise_u(t, breakpts, Uvals)
    idx = discretize(t, breakpts, 'IncludedEdge','left');  % which interval?
    idx(isnan(idx)) = size(Uvals,2);                       % t == breakpts(end)
    u   = Uvals(:, idx);
end

function poly = generate_polynomials(pce_order, theta1, theta2)
% Return every Legendre-tensor product   Pn(theta1) · Pm(theta2)
% whose total degree n+m is ≤ pce_order.

    poly = sym([]);                 % column vector we’ll keep appending to

    for n = 0:pce_order
        Pn = sqrt(2*n+1) * legendreP(n, theta1);     % ϕn(θ₁)

        for m = 0:(pce_order - n)                    % only combos with n+m ≤ pce_order
            Pm = sqrt(2*m+1) * legendreP(m, theta2); % ϕm(θ₂)

            poly(end+1,1) = expand(Pn * Pm);         % store as a column vector
        end
    end
end

function u_proj = projection(u,flat_reps,n)
u_proj = zeros(n,1);
    for i = 1:length(flat_reps)
        u_proj(flat_reps(i)) = u(i);
    end
end


function plot_trajectories(y_store_all, ysp, filename, col , faults)
    [xxx xxxx iterations] = size(y_store_all);
    t  = 0:400;                      % time vector (adjust if needed)
         % base colour for the trajectories
    
    figure;
    tl = tiledlayout(3, 2, 'TileSpacing','compact', 'Padding','compact');
    
    for r = 1:6
        nexttile; hold on;
    
        % overlay 100 trajectories for this row
        for it = 1:iterations
            plot(t, squeeze(y_store_all(r, :, it)), ...
                 'Color',col);    % same hue, very light (α ≈ 0.15)
        end

        plot(t(1:end-1), squeeze(ysp(r, :)), 'Color',[0 0 1]);
        xlabel('t');  ylabel(sprintf('Cluster %d [K]', r));
        ylim([-2,2])
        if faults
            xregion(50,150,'DisplayName',"First Input Fault","FaceColor",'r')
            xregion(250,350,'DisplayName',"Second Input Fault","FaceColor",'r')
        end
    end

    print(gcf, filename, '-depsc2', '-vector'); 
end




function plot_ieee118_clusters(G, new_idx, representatives, filename)
% plot_ieee118_clusters  Plot IEEE-118 network coloured by cluster
%
%   plot_ieee118_clusters(G, new_idx, representatives, filename)
%
%   G               graph/digraph object
%   new_idx         cluster label for each node (vector)
%   representatives vector of cluster-representative node IDs (zeros allowed)
%   filename        (optional) name for EPS file *without* extension
%
% Example:
%   plot_ieee118_clusters(G, new_idx, reps, 'ieee118_clustered')

    figure('Color','w');

    % -------- core plot --------------------------------------------------
    p = plot(G, ...
        'Layout',     'force', ...
        'NodeCData',  new_idx, ...
        'MarkerSize', 6, ...
        'EdgeAlpha',  0.30, ...
        'NodeLabel', {}, ...
        'LineWidth',  0.5);

    colormap(jet(max(new_idx)));
    axis off
    set(gca,'FontSize',12)

    % -------- highlight representatives ----------------------------------
    flat_reps = nonzeros(representatives(:));
    % enlarge marker + dark edge, but DON'T change NodeColor
    highlight(p, flat_reps, ...
        'Marker',           'o',  ...
        'MarkerSize',       9,    ...
        'MarkerEdgeColor',  'k',  ...
        'LineWidth',        1);

    title('IEEE 118-Bus Network', 'FontSize',14,'FontWeight','bold');

    % -------- optional EPS export ----------------------------------------
    if nargin >= 4 && ~isempty(filename)
        print(gcf, filename, '-depsc2', '-painters');  % vector EPS
    end
end