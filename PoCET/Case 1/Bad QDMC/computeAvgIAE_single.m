function avgIAE = computeAvgIAE_single(y_store, y_sp, dt)
% computeAvgIAE_single  Average IAE for a single‑cluster case
%
%   avgIAE = computeAvgIAE_single(y_store, y_sp)
%   avgIAE = computeAvgIAE_single(y_store, y_sp, dt)
%
% Inputs
%   y_store : [n_t × n_iter]  output trajectories (one cluster)
%   y_sp    : [n_t × 1]       set‑point trajectory
%   dt      : (optional)      sample time (default = 1)
%
% Output
%   avgIAE  : scalar          mean IAE across all realizations
%
% -------------------------------------------------------------------------

if nargin < 3
    dt = 1;                      % default sample period
end

% |error| for every time sample and realization
absErr = abs(y_store - y_sp);    % implicit expansion (R2016b+)

% Integrate |error| over time for each realization
IAE = trapz(absErr, 1) * dt;     % 1×n_iter row vector

% Average across realizations
avgIAE = mean(IAE);
end