function width = computeEnvelopeWidth_single(y_store)
% computeEnvelopeWidth_single  Computes the envelope width at final time
%                              for a single-cluster case
%
%   width = computeEnvelopeWidth_single(y_store)
%
% Inputs
%   y_store : [n_t × n_iter] output trajectories
%
% Output
%   width   : scalar envelope width at final time across realizations
%
% -------------------------------------------------------------------------

% Extract outputs at final time
final_outputs = y_store(end, :);  % [1 × n_iter]

% Compute envelope width
width = max(final_outputs) - min(final_outputs);
end