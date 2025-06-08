function avgWidth = computeEnvelopeWidth(y_store)
% computeEnvelopeWidth  Computes the average final-time envelope width
%
%   avgWidth = computeEnvelopeWidth(y_store)
%
% Inputs
%   y_store : [k × n_t × n_iter] output trajectories
%
% Output
%   avgWidth : scalar average width across all clusters at final time
%
% -------------------------------------------------------------------------

[k, n_t, n_iter] = size(y_store);
final_outputs = squeeze(y_store(:, end, :));  % [k × n_iter]

% Compute max - min across realizations for each cluster
widths = max(final_outputs, [], 2) - min(final_outputs, [], 2);  % [k × 1]

% Average across all clusters
avgWidth = mean(widths);
end