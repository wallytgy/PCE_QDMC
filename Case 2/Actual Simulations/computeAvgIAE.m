function avg_IAE = computeAvgIAE(y_store, y_sp , dt)
% Computes the average integral absolute error (IAE) over all clusters and realizations
%
% Inputs:
%   y_store : [k x n_t x n_iter] array of output trajectories
%   y_sp    : [k x n_t] array of setpoints
%
% Output:
%   avg_IAE : scalar average IAE across all clusters and realizations

[k, n_t, n_iter] = size(y_store);
IAE_total = 0;

for i = 1:n_iter
    y_i = y_store(:,:,i);                % [k x n_t]
    abs_err = abs(y_i - y_sp);           % element-wise absolute error
    IAE_i = trapz(1:n_t, abs_err, 2) * dt;    % integrate over time (dim 2) per cluster
    IAE_total = IAE_total + sum(IAE_i);  % sum across clusters and accumulate
end

avg_IAE = IAE_total / (k * n_iter);      % average over all clusters and realizations
end