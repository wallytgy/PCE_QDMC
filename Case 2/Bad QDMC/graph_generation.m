clear all
% make the graph and plot the graph to see the topology
case118data = pglib_opf_case118_ieee
branch = case118data.branch
edges = branch(:,1:2)
x = branch(:,4)
G = graph(edges(:,1),edges(:,2),x)
figure(1)
plot(G, ...
    'MarkerSize', 3, ...
    'LineWidth', 0.5, ...
    'Layout', 'force'); 

% calculate the impedences
Bij = 1./x
Bij = Bij / mean(abs(Bij))
Vref = 1
M = 1
eta1 = 1
[A,B] = build_AB(edges, Bij, Vref, M, eta1)



%%

% Step 1: Get weighted Laplacian
Gsimp = simplify(G)
Adj = adjacency(Gsimp, 'weighted');  % Weighted adjacency matrix
D = diag(sum(Adj, 2));           % Degree matrix
L = D - Adj;    

% Step 2: Choose number of clusters
k = 10;  % adjust as needed

% Step 3: Compute first k eigenvectors of Laplacian
[V, ~] = eigs(L, k, 'smallestabs'); 

% Step 4: Cluster rows of eigenvector matrix
idx = kmeans(V, k); 

% Step 5: Plot
figure(2)
plot(G, 'Layout', 'force', 'NodeCData', idx, 'MarkerSize', 5);

