function [A,B] = build_AB(edges, Bij, Vref, M, eta2, actuators)
% BUILD_AB  Construct A and B for  \dot x = A x + B u
%
% edges     : m-by-2   edge list  [i  j] (undirected)
% Bij       : m-by-1   susceptances (same order as edges)
% Vref, M   : scalars  (for k_ij = Bij * Vref / M^2)
% eta1      : scalar   input gain
% actuators : vector   node indices that have control inputs
%                         (default: all nodes)
%
% Returns
%   A : n-by-n
%   B : n-by-|actuators|

% --------- constants & sizes -----------------------------------------
n     = max(edges(:));          % number of buses/nodes
if nargin < 6 || isempty(actuators)
    actuators = 1:n;            % default: actuator on every node
end
m     = size(edges,1);

% --------- edge weights  k_ij  ---------------------------------------
k     = Bij .* (Vref / M^2);     % m-by-1 vector

% --------- sparse Laplacian  L  --------------------------------------
i_idx = edges(:,1); j_idx = edges(:,2);

L = sparse( ...
    [i_idx; j_idx], ...
    [j_idx; i_idx], ...
    [-k;    -k], ...
    n, n);

diag_vals = accumarray([i_idx; j_idx], [k; k], [n 1]);
L = L + spdiags(diag_vals, 0, n, n);

% --------- A and B ----------------------------------------------------
A = -L;

B = eta2 * sparse(1:n, actuators, 1, n, numel(actuators));

end