classdef QDMC_controller_soft < handle
% Quadratic Dynamic Matrix Control for MIMO System with soft constraints
%
% This function creates and runs a quadratic dynamic matrix control
% algorithm for MIMO systems. See the note on "Model Predictive Control and
% State Estimation" by Enso Ikonen (2013).
% It determines the input changes based on the plant model (step response),
% current measurement, and current disturbance deviations w.r.t. to the
% preceeding step.
%
% Copied in great parts from Joel Paulson and Ali Mesbah, 3-13-2013
%
% Revised and improved by Matthias von Andrian 05/2016
%
% Transfered to classdef mode by Matthias von Andrian 06/2017
%
% Soft constraints added by Matthias von Andrian 06/2018
% a-n-d-r-i-a-n[at]m-i-t.e-d-u  -> remove the "-" signs
% For questions and to report bugs, please email me.
%
% For accademic use, the rules of properly citing other people's work apply


properties (GetAccess = public, SetAccess = public)
% Stored parameters, assigned in the initialization function
% QDMC_controller and used in the run function run_controller.
% With full access outside of class
    u_past      % Control input implemented in the previous time step,
                % vector [n_in by 1]
    d_past      % Measured disturbance in the previous time step,
                % vector [n_dis by 1]
    ops         % Options for quadratic program solver,
                % optimoptions('quadprog',...)
    F           % Free response -> Output, if no further controler action
                % was taken, vector [n*n_out by 1]
end % End of public properties

properties (GetAccess = public, SetAccess = protected)
% Can not be overriten
	n_in        % Number of input variables, scalar
    n_out       % Number of output variables, scalar
    n_dis       % Number of measured disturbances, scalar
    n           % Number of step response coefficients, scalar
    Su          % Step response coefficients input, matrix,
                % [n_out by n_in by n] -> 3D matrix
    Sd          % Step response coefficients disturbances, matrix,
                % [n_out by n_dis by n] -> 3D matrix
    hu          % Slope between the last two input step response
                % coefficients (assumed to be constant, assumed to be 0 for
                % nonintegrating pairings) matrix [n_out by n_in]
    hd          % Slope between the last two disturbance step response
                % coefficients (assumed to be constant, assumed to be 0 for
                % nonintegrating pairings) matrix [n_out by n_dis]
    p           % Prediction horizon, scalar
    c           % Control horizon, scalar
    u_min       % Minimum input, vector [n_in by 1]
    u_max       % Maximum input, vector [n_in by 1]
    D_u_min     % Minimum input change, vector [n_in by 1]
    D_u_max     % Maximum input change, vector [n_in by 1]
    y_min       % Minimum output, vector [n_out by 1]
    y_max       % Maximum output, vector [n_out by 1]
    soft        % Index of soft constraints [n_out by 2]
                % First column for ymax, second column for ymin
    T           % Matrix to extract free response within prediction horizon
                % binary matrix [p by n]
    G           % Dynamic matrix, combination of lower triangular matrices
                % with elements of Su, matrix [p*n_out by c*n_in]
    Wu          % Diagonal representation of the control weight parameters
    Wy          % Diagonal representation of the output weight parameters
    w_eps       % Vector of weights for slack variables for soft output
                % constraints in objective function (linear term)
    W_eps       % Diagonal matrix of weights for slack variables for soft 
                % output constraints in objective function (quadratic term)            
    H           % Hessian of the quadratic program
    Hsym        % Symmetric Hessian of the quadratic program
    M           % Matrix to update free response, binary matrix [n by n]
                % For soft constraits, decision variable becomes
                % x = [Delta_u; eps_max; eps_min]
    A_u         % Matrix to handle inequality constraints 
                % A_u * Delta_u <= b
                % -> has to be extended, if output constraints are softened
    A_eps       % Matrix to handle inequality constaints if output 
                % constraints are softened [A_u A_eps] * x <= b
    A           % Concatenated [A_u A_eps]
    lb_u        % Lower bound of Delta_u
                % -> has to be extended, if output constraints are softened
    lb_eps      % Lower bound of slack variables epsilon_max, epsilon_min
                % -> lower bound of x: [lb_u; lb_eps]
    ub_u        % Upper bound of Delta_u
                % -> has to be extended, if output constraints are softened
    ub_eps      % Upper bound of slack variables epsilon_max, epsilon_min
                % -> upper bound of x: [ub_u; ub_eps]
    lb          % Concatenated [lb_u; lb_eps];
    ub          % Concatenated [ub_u; ub_eps];
end % End of protected properties


methods

% ======================================================================= %
% Initialization                                           Initialization %
% ======================================================================= %
function obj = QDMC_controller_soft(n_in,n_out,n_dis,Su,Sd,p,c,La,Q,ctg,...
                               u_min,u_max,D_u_min,D_u_max,y_min,y_max,...
                               u_past,y_past,d_past,u_int,d_int,soft,w_eps,W_eps)
% This function checks for consistant inputs and brings them into the
% required format and creates structures needed to run the QDMC controller.
%
% Inputs (* marks optional input, leave empty [] if not applicable)
%   n_in        Number of input variables, scalar
%   n_out       Number of output variables, scalar
% * n_dis *     Number of measured disturbances, scalar
%   Su          Step response coefficients input, matrix,
%               [n_out by n_in by n] -> 3D matrix
% * Sd *        Step response coefficients disturbances, matrix,
%               [n_out by n_dis by n] -> 3D matrix
%   p           Prediction horizon, scalar
%   c           Control horizon, scalar
% * La *        Vector of control weight parameters, 
%               i.e., J = ||r-y|| + La*eye*||du|| where du is the input
%               change, [n_in by 1]
% * Q *         Vector of output weight parameters, vector [n_out by 1]
% * ctg *       Vector of cost to go parameters (additional weight for the
%               last predicted output error, vector [n_out by 1]
% * u_min *     Minimum input, vector [n_in by 1]
% * u_max *     Maximum input, vector [n_in by 1]
% * D_u_min *   Minimum input change, vector [n_in by 1]
% * D_u_max *   Maximum input change, vector [n_in by 1]
% * y_min *     Minimum output, vector [n_out by 1]
% * y_max *     Maximum output, vector [n_out by 1]
% * u_past *    Control input implemented in the previous time step,
%               vector [n_in by 1]
% * y_past *    Past output in the previous time step, vector [n_out by 1]
% * d_past *    Measured disturbance in the previous time step,
%               vector [n_dis by 1]
% * u_int *     Is Input integrating on specific output?
%               Binary matrix [n_out by n_in]
%               If not supplied, it will be set to 0
% * d_int *     Is disturbance integrating on specific output?
%               Binary matrix [n_out by n_dis]
%               If not supplied, it will be set to 0
% * soft *      Are output constraints soft constraints?
%               Binary matrix [n_out by 2], First column for max, second
%               for min. If not supplied, it will be set to 0
% * w_eps *     Weight for slack variable epsilon in objective function
%               J = ... + w_eps' * epsilon, vector [n_out by 1]
% * W_eps *     Weight for slack variable epsilon in objective function
%               J = ... + epsilon' + W_eps' * epsilon, vector [n_out by 1]

% Check input =========================================================== %
if nargin < 24
    error('Not enough input arguments')
end

% Number of input, output and measured disturbance variables ------------ %
if isempty(n_in) || ~isscalar(n_in)
    error('Please specify the number of inputs n_in as scalar')
end
if isempty(n_out) || ~isscalar(n_out)
    error('Please specify the number of outputs n_out as scalar')
end
if isempty(n_dis)
    n_dis = 0;
elseif ~isscalar(n_dis)
    error(['Please specify the number of measured disturbances n_dis ',...
           'as scalar'])
end

% Prediction and control horizon ---------------------------------------- %
if isempty(p) || ~isscalar(p)
    error('Please specify the prediction horizon p as scalar')
end
if isempty(c) || ~isscalar(c)
    error('Please specify the control horizon c as scalar')
end
if c >= p
    error(['Please specify the control horizon c smaller than the ',...
           'prediction horizon p'])
end

% Constraints ----------------------------------------------------------- %
if ~isempty(u_min) && (size(u_min,1) ~= n_in || size(u_min,2) ~= 1)
    error('u_min has to be a column vector with n_in entries')
end
if ~isempty(u_max) && (size(u_max,1) ~= n_in || size(u_max,2) ~= 1)
    error('u_max has to be a column vector with n_in entries')
end
if ~isempty(D_u_min) && (size(D_u_min,1) ~= n_in || size(D_u_min,2) ~= 1)
    error('D_u_min has to be a column vector with n_in entries')
end
if ~isempty(D_u_max) && (size(D_u_max,1) ~= n_in || size(D_u_max,2) ~= 1)
    error('D_u_max has to be a column vector with n_in entries')
end
if ~isempty(y_min) && (size(y_min,1) ~= n_out || size(y_min,2) ~= 1)
    error('y_min has to be a column vector with n_out entries')
end
if ~isempty(y_max) && (size(y_max,1) ~= n_out || size(y_max,2) ~= 1)
    error('y_max has to be a column vector with n_out entries')
end


% Past measurements ----------------------------------------------------- %
if ~isempty(u_past) && (size(u_past,1) ~= n_in || size(u_past,2) ~= 1)
    error('u_past has to be a column vector with n_in entries')
end
if isempty(u_past)
    u_past = zeros(n_in,1); % Assume steady state
end
if ~isempty(y_past) && (size(y_past,1) ~= n_out || size(y_past,2) ~= 1)
    error('y_past has to be a column vector with n_out entries')
end
if isempty(y_past)
    y_past = zeros(n_out,1); % Assume steady state
end
if ~isempty(d_past) && (size(d_past,1) ~= n_dis || size(d_past,2) ~= 1)
    error('d_past has to be a column vector with n_dis entries')
end
if isempty(d_past) && n_dis > 0
    d_past = zeros(n_dis,1); % Assume steady state
end

% Integrating inputs ---------------------------------------------------- %
if ~isempty(u_int) && (size(u_int,1) ~= n_out || size(u_int,2) ~= n_in)
    error('u_int has to be a matrix with dimensions [n_out by n_in]')
end
if isempty(u_int)
    u_int = zeros(n_out,n_in);
end
u_int = ~~u_int; % Make logical/binary
if ~isempty(d_int) && (size(d_int,1) ~= n_out || size(d_int,2) ~= n_dis)
    error('d_int has to be a matrix with dimensions [n_out by n_dis]')
end
if isempty(u_int)
    d_int = zeros(n_out,n_dis);
end
d_int = ~~d_int; % Make logical/binary

% Soft constraints ------------------------------------------------------ %
if ~isempty(soft) && (size(soft,1) ~= n_out || size(soft,2) ~= 2)
    error('soft has to be a matrix with dimensions [n_out by 2]')
end
if isempty(soft)
    soft = zeros(n_out,2);
end
if isempty(y_max)
    soft(:,1) = 0;
end
if isempty(y_min)
    soft(:,2) = 0;
end
soft = ~~soft; % Make logical/binary

% Weights --------------------------------------------------------------- %
if ~isempty(La) && (size(La,1) ~= n_in || size(La,2) ~= 1)
    error('La has to be a column vector with n_in entries')
end
if isempty(La)
    La = ones(n_in,1);
end
if ~isempty(Q) && length(Q(:)) ~= n_out
    error('Q has to be a vector with n_out entries')
end
if isempty(Q)
    Q = ones(n_out,1);
end
if ~isempty(ctg) && length(ctg(:)) ~= n_out
    error('ctg has to be a vector with n_out entries')
end
if isempty(ctg)
    ctg = zeros(n_out,1);
end
if sum(sum(soft)) > 0 && isempty(w_eps)
    warning('w_eps will be set to 100 * Q')
    w_eps = 100 * Q;
end
if sum(sum(soft)) > 0 && isempty(W_eps)
    warning('W_eps will be set to 100 * Q')
    W_eps = 100 * Q;
end
if ~isempty(w_eps) && (size(w_eps,1) ~= n_out || size(w_eps,2) ~= 1)
    error('w_eps has to be a vector with dimensions [n_out by 1]')
end
if ~isempty(W_eps) && (size(W_eps,1) ~= n_out || size(W_eps,2) ~= 1)
    error('W_eps has to be a vector with dimensions [n_out by 1]')
end
if sum(sum(soft)) == 0 && ~isempty(w_eps)
    warning('no soft constraints are considered, but w_eps is defined')
end
if sum(sum(soft)) == 0 && ~isempty(W_eps)
    warning('no soft constraints are considered, but W_eps is defined')
end

% Step response matrices ------------------------------------------------ %
[d1, d2, ~] = size(Su);
if d1 ~= n_out || d2 ~= n_in
    error(['Step response coefficients for input Su have to be in a ',...
           'matrix with dimensions [n_out by n_in by n]'])
end
if n_dis > 0
    [d1, d2, ~] = size(Sd);
    if d1 ~= n_out || d2 ~= n_dis
        error(['Step response coefficients for measured disturbances ',...
               'Sd have to be in a matrix with dimensions ',...
               '[n_out by n_dis by n]'])
    end
else
    Sd = [];
end
% Make sure, number of step response coefficients is equal for input
% and disturbance
if n_dis > 0
    [~, ~, d3u] = size(Su);
    [~, ~, d3d] = size(Sd);
    if d3u > d3d
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fill in Su so it becomes the same length as Sd
    elseif d3u < d3d
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fill in Sd so it becomes the same length as Su
    end
end
% Extend step response coefficients, if shorter than prediction horizon
n = size(Su,3);
if n <= p
    % Copy last value for non integrating, keep constant slope for
    % integrating
    % Input - output
    hu = (Su(:,:,n) - Su(:,:,n-1)) .* u_int;
    Su_extended = repmat(Su(:,:,n),1,1,p-n+1) + ...
        (repmat(hu,1,1,p-n+1) .* ...
         repmat(permute(1:p-n+1,[1,3,2]),n_out,n_in,1));
    Su = cat(3,Su,Su_extended);
    % Disturbance - output
    if n_dis > 0
        hd = (Sd(:,:,n) - Sd(:,:,n-1)) .* d_int;
        Sd_extended = repmat(Sd(:,:,n),1,1,p-n+1) + ...
            (repmat(hd,1,1,p-n+1) .* ...
             repmat(permute(1:p-n+1,[1,3,2]),n_out,n_dis,1));
        Sd = cat(3,Sd,Sd_extended);
    end
    n = p + 1;
end


% Set parameters (properties) =========================================== %
% Number of variables --------------------------------------------------- n
obj.n_in  = n_in;
obj.n_out = n_out;
obj.n_dis = n_dis;
obj.n = n;

% Step response matrices -------------------------------------------- Su Sn
obj.Su = Su;
obj.Sd = Sd;

% Difference between last step response coefficients ---------------- hu hd
obj.hu = (Su(:,:,n) - Su(:,:,n-1)) .* u_int;
if n_dis > 0
    obj.hd = (Sd(:,:,n) - Sd(:,:,n-1)) .* d_int;
else
    obj.hd = [];
end

% Horizons ------------------------------------------------------------ p c
obj.p = p;
obj.c = c;

% Constraints ------------------------------------------------- Constraints
obj.u_min = u_min;
obj.u_max = u_max;
obj.D_u_min = D_u_min;
obj.D_u_max = D_u_max;
obj.y_min = y_min;
obj.y_max = y_max;

% Soft constraints index --------------------------------------------- soft
obj.soft = soft;

% Past measurements ----------------------------------------- u_past d_past
obj.u_past = u_past;
obj.d_past = d_past;

% Determine matrix T ---------------------------------------------------- T
obj.T = ...
    sparse(p,n) + [sparse(p,1) spdiags(ones(p,1),0,p,p) sparse(p,n-p-1)];

% Determine matrix G ---------------------------------------------------- G
obj.G = zeros(p*n_out,c*n_in);
for ii = 1:n_out
    for jj = 1:n_in
        Su_part = permute(Su(ii,jj,:),[3,1,2]);
        obj.G((ii-1)*p+1:ii*p, (jj-1)*c+1:jj*c) = ...
            toeplitz(Su_part(1:p),Su_part(1,1)*eye(1,c));
    end
end

% Determine initial free response F ------------------------------------- F
obj.F = reshape(repmat(y_past,1,obj.n)',[],1);

% Determine input weighting matrix Wu ---------------------------------- Wu
obj.Wu = spdiags(...
          reshape(...
           repmat(La,1,c)',...
          n_in*c,1),...
         0,n_in*c,n_in*c);

% Determine output weighting matrix Wy --------------------------------- Wy
obj.Wy = spdiags(...
          reshape(...
           repmat(Q,1,p)' + [zeros(p-1,n_out);ctg'],...
          n_out*p,1),...
         0,n_out*p,n_out*p);

% Determine weighting vector and matrix w_eps and W_eps ------ w_eps, W_eps
obj.w_eps = [];
if sum(sum(soft)) > 0
    if ~isempty(y_max) && sum(soft(:,1)) > 0
        obj.w_eps = [obj.w_eps; reshape(repmat(w_eps .* soft(:,1),1,p)',[],1)];
    end
    if ~isempty(y_min) && sum(soft(:,2)) > 0
        obj.w_eps = [obj.w_eps; reshape(repmat(w_eps .* soft(:,2),1,p)',[],1)];
    end
end
if sum(sum(soft)) == 0
    obj.W_eps = [];
else
    W_eps_diag = [];
    if ~isempty(y_max) && sum(soft(:,1)) > 0
        W_eps_diag = [W_eps_diag; reshape(repmat(W_eps .* soft(:,1),1,p)',[],1)];
    end
    if ~isempty(y_min) && sum(soft(:,2)) > 0
        W_eps_diag = [W_eps_diag; reshape(repmat(W_eps .* soft(:,2),1,p)',[],1)];
    end
    obj.W_eps = spdiags(W_eps_diag,0,length(W_eps_diag),length(W_eps_diag));
end

% Calculate Hessian H --------------------------------------------------- H
obj.H = obj.G' * obj.Wy * obj.G + obj.Wu;
obj.Hsym = (obj.H+obj.H')/2;  % Quadprog requires symmetric hessian
    
% Determine matrix M to update free response ---------------------------- M
obj.M = spdiags(ones(n,1),1,n,n); obj.M(n,n) = 1;

% Set options for quadprog ops ---------------------------------------- ops
obj.ops = optimoptions('Quadprog'...
                      ,'Display','off'...
                      ,'Algorithm','interior-point-convex'...
                      ,'MaxIterations',1e4...
                      );

% Set matrix for constraint handling A_u ------------------------------ A_u
% Constraints: A_u * Delta_u <= b
% b is dependent on future delta u and is determined in run_controller
A_u = [];
% Use lower and upper bounds instead
% if ~isempty(obj.D_u_max)
%     A1max = eye(obj.n_in*obj.c);
%     A = [A; A1max];
% end
% if ~isempty(obj.D_u_min)
%     A1min = -eye(obj.n_in*obj.c);
%     A = [A; A1min];
% end

if ~isempty(obj.u_max)
    A2 = tril(ones(obj.c));
    A2max = [];
    for jj = 1:obj.n_in
        A2max = blkdiag(A2max,A2);
    end
    A_u = [A_u; A2max];
end
if ~isempty(obj.u_min)
    A2 = tril(ones(obj.c));
    A2min = [];
    for jj = 1:obj.n_in
        A2min = blkdiag(A2min,-A2);
    end        
    A_u = [A_u; A2min];
end

if ~isempty(obj.y_max)
    A3max = obj.G;
    A_u = [A_u; A3max];
end
if ~isempty(obj.y_min)
    A3min = -obj.G;
    A_u = [A_u; A3min];
end
obj.A_u = A_u;

% Set matrix for constraint handling A_eps -------------------------- A_eps
% Subtract slack variable epsilon from inequality constraints to soften:
% A_u * Delta_u - epsilon <= b  ->  [A_u A_eps] * [Delta_u; epsilon] <= b
% b is dependent on future delta u and is determined in run_controller
A_eps = [];

if ~isempty(obj.u_max)
    A2max = [];
    if ~isempty(obj.y_max) && sum(soft(:,1)) > 0
        Amax = sparse(n_in*c, n_out*p);
        A2max = [A2max Amax];
    end
    if ~isempty(obj.y_min) && sum(soft(:,2)) > 0
        Amin = sparse(n_in*c, n_out*p);
        A2max = [A2max Amin];
    end
    A_eps = [A_eps; A2max];
end
if ~isempty(obj.u_min)
    A2min = [];
    if ~isempty(obj.y_max) && sum(soft(:,1)) > 0
        Amax = sparse(n_in*c, n_out*p);
        A2min = [A2min Amax];
    end
    if ~isempty(obj.y_min) && sum(soft(:,2)) > 0
        Amin = sparse(n_in*c, n_out*p);
        A2min = [A2min Amin];
    end
    A_eps = [A_eps; A2min];
end

if ~isempty(obj.y_max)
    A3max = [];
    if sum(obj.soft(:,1)) > 0
        A3max = [A3max -spdiags(...
                         reshape(...
                          repmat(obj.soft(:,1),1,obj.p)',...
                         [],1),...
                        0,obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    if ~isempty(obj.y_min) && sum(obj.soft(:,2)) > 0
        A3max = [A3max sparse(obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    A_eps = [A_eps; A3max];
end
if ~isempty(obj.y_min)
    A3min = [];
    if ~isempty(obj.y_max) && sum(obj.soft(:,1)) > 0
        A3min = [A3min sparse(obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    if sum(obj.soft(:,2)) > 0
        A3min = [A3min -spdiags(...
                         reshape(...
                          repmat(obj.soft(:,2),1,obj.p)',...
                         [],1),...
                        0,obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    A_eps = [A_eps; A3min];
end
obj.A_eps = A_eps;

% Concatenate A --------------------------------------------------------- A
obj.A = [obj.A_u obj.A_eps];

% Set lower and upper bounds for Delta_u ----------------------- lb_u, ub_u
if ~isempty(D_u_min)
    obj.lb_u = reshape(repmat(D_u_min,1,c)',c*n_in,1);
else
    obj.lb_u = -inf(n_in*c,1);
end
if ~isempty(D_u_max)
    obj.ub_u = reshape(repmat(D_u_max,1,c)',c*n_in,1);
else
    obj.ub_u = inf(n_in*c,1);
end

% Set lower and upper bounds for slack variable epsilon ---- lb_eps, ub_eps
obj.lb_eps = [];
obj.ub_eps = [];
if ~isempty(y_max) && sum(soft(:,1)) > 0
    obj.lb_eps = [obj.lb_eps; sparse(n_out*p,1)];
    obj.ub_eps = [obj.ub_eps; inf(n_out*p,1)];
end
if ~isempty(y_min) && sum(soft(:,2)) > 0
    obj.lb_eps = [obj.lb_eps; sparse(n_out*p,1)];
    obj.ub_eps = [obj.ub_eps; inf(n_out*p,1)];
end

% Concatenate lb and ub --------------------------------------------- lb ub
obj.lb = [obj.lb_u; obj.lb_eps];
obj.ub = [obj.ub_u; obj.ub_eps];

end % End of initialization function QDMC_controller



% ======================================================================= %
% update_constraints                                   update_constraints %
% ======================================================================= %
function obj = update_constraints(obj,u_min,u_max,D_u_min,D_u_max,y_min,y_max)
% Update the matrix A_u to handle the constraints A_u * Delta_u <= b
% and the matrix A_eps to handle the constraints 
% [A_u A_eps] * [Delta_u; epsilon] <= b

% Check input
if ~isempty(u_min) && length(u_min(:)) ~= obj.n_in
    error('u_min has to be a vector with n_in entries')
end
if ~isempty(u_max) && length(u_max(:)) ~= obj.n_in
    error('u_max has to be a vector with n_in entries')
end
if ~isempty(D_u_min) && length(D_u_min(:)) ~= obj.n_in
    error('D_u_min has to be a vector with n_in entries')
end
if ~isempty(D_u_max) && length(D_u_max(:)) ~= obj.n_in
    error('D_u_max has to be a vector with n_in entries')
end
if ~isempty(y_min) && length(y_min(:)) ~= obj.n_out
    error('y_min has to be a vector with n_out entries')
end
if ~isempty(y_max) && length(y_max(:)) ~= obj.n_out
    error('y_max has to be a vector with n_out entries')
end

% Update constraints
obj.u_min = u_min;
obj.u_max = u_max;
obj.D_u_min = D_u_min;
obj.D_u_max = D_u_max;
obj.y_min = y_min;
obj.y_max = y_max;

% Set matrix for constraint handling A_u 
% Constraints: A_u * Delta_u <= b
% b is dependent on future delta u and is determined in run_controller
a_u = []; % A_u triggers a warning about using the same name as obj.A_u
% Use lower and upper bounds instead
% if ~isempty(obj.D_u_max)
%     A1max = eye(obj.n_in*obj.c);
%     A = [A; A1max];
% end
% if ~isempty(obj.D_u_min)
%     A1min = -eye(obj.n_in*obj.c);
%     A = [A; A1min];
% end

if ~isempty(obj.u_max)
    A2 = tril(ones(obj.c));
    A2max = [];
    for jj = 1:obj.n_in
        A2max = blkdiag(A2max,A2);
    end
    a_u = [a_u; A2max];
end
if ~isempty(obj.u_min)
    A2 = tril(ones(obj.c));
    A2min = [];
    for jj = 1:obj.n_in
        A2min = blkdiag(A2min,-A2);
    end        
    a_u = [a_u; A2min];
end

if ~isempty(obj.y_max)
    A3max = obj.G;
    a_u = [a_u; A3max];
end
if ~isempty(obj.y_min)
    A3min = -obj.G;
    a_u = [a_u; A3min];
end
obj.A_u = a_u;

% Set matrix for constraint handling A_eps -------------------------- A_eps
% Subtract slack variable epsilon from inequality constraints to soften:
% A_u * Delta_u - epsilon <= b  ->  [A_u A_eps] * [Delta_u; epsilon] <= b
% b is dependent on future delta u and is determined in run_controller
a_eps = []; % A_eps triggers a warning about using the same name as obj.A_eps

if ~isempty(obj.u_max)
    A2max = [];
    if ~isempty(obj.y_max) && sum(obj.soft(:,1)) > 0
        Amax = sparse(obj.n_in*obj.c, obj.n_out*obj.p);
        A2max = [A2max Amax];
    end
    if ~isempty(obj.y_min) && sum(obj.soft(:,2)) > 0
        Amin = sparse(obj.n_in*obj.c, obj.n_out*obj.p);
        A2max = [A2max Amin];
    end
    a_eps = [a_eps; A2max];
end
if ~isempty(obj.u_min)
    A2min = [];
    if ~isempty(obj.y_max) && sum(obj.soft(:,1)) > 0
        Amax = sparse(obj.n_in*obj.c, obj.n_out*obj.p);
        A2min = [A2min Amax];
    end
    if ~isempty(obj.y_min) && sum(obj.soft(:,2)) > 0
        Amin = sparse(obj.n_in*obj.c, obj.n_out*obj.p);
        A2min = [A2min Amin];
    end
    a_eps = [a_eps; A2min];
end

if ~isempty(obj.y_max)
    A3max = [];
    if sum(obj.soft(:,1)) > 0
        A3max = [A3max -spdiags(...
                         reshape(...
                          repmat(obj.soft(:,1),1,obj.p)',...
                         [],1),...
                        0,obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    if ~isempty(obj.y_min) && sum(obj.soft(:,2)) > 0
        A3max = [A3max sparse(obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    a_eps = [a_eps; A3max];
end
if ~isempty(obj.y_min)
    A3min = [];
    if ~isempty(obj.y_max) && sum(obj.soft(:,1)) > 0
        A3min = [A3min sparse(obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    if sum(obj.soft(:,2)) > 0
        A3min = [A3min -spdiags(...
                         reshape(...
                          repmat(obj.soft(:,2),1,obj.p)',...
                         [],1),...
                        0,obj.n_out*obj.p,obj.n_out*obj.p)];
    end
    a_eps = [a_eps; A3min];
end
obj.A_eps = a_eps;

% Concatenate A --------------------------------------------------------- A
obj.A = [obj.A_u obj.A_eps];

% Set lower and upper bounds for Delta_u ----------------------- lb_u, ub_u
if ~isempty(D_u_min)
    obj.lb_u = reshape(repmat(D_u_min,1,obj.c)',obj.c*obj.n_in,1);
else
    obj.lb_u = -inf(obj.n_in*obj.c,1);
end
if ~isempty(D_u_max)
    obj.ub_u = reshape(repmat(D_u_max,1,obj.c)',obj.c*obj.n_in,1);
else
    obj.ub_u = inf(obj.n_in*obj.c,1);
end

% Set lower and upper bounds for slack variable epsilon ---- lb_eps, ub_eps
obj.lb_eps = [];
obj.ub_eps = [];
if ~isempty(y_max) && sum(obj.soft(:,1)) > 0
    obj.lb_eps = [obj.lb_eps; sparse(obj.n_out*obj.p,1)];
    obj.ub_eps = [obj.ub_eps; inf(obj.n_out*obj.p,1)];
end
if ~isempty(y_min) && sum(obj.soft(:,2)) > 0
    obj.lb_eps = [obj.lb_eps; sparse(obj.n_out*obj.p,1)];
    obj.ub_eps = [obj.ub_eps; inf(obj.n_out*obj.p,1)];
end

% Concatenate lb and ub --------------------------------------------- lb ub
obj.lb = [obj.lb_u; obj.lb_eps];
obj.ub = [obj.ub_u; obj.ub_eps];

end % End of function update_constraints



% ======================================================================= %
% Run controller                                          Rund controller %
% ======================================================================= %
function [u, J] = run_controller(obj,y,y_sp,d)
%

% Check inputs 
if length(y(:)) ~= obj.n_out
    error('y has to be a vector with n_out entries')
end
if size(y_sp,1) ~= obj.n_out
    error(['Dimensions of y_sp do not match. 2 options for y_sp: ',...
           'vector [n_out by 1] -> Constant setpoint or ',...
           'matrix [n_out by p or more] -> setpoint changing over horizon'])
end
if ~isempty(d) && length(d(:)) ~= obj.n_dis
    error('d has to be a vector with n_dis entries')
end

% Predicted outputs due to past inputs and past and present disturbances  %
yp = zeros(obj.n_out*obj.p,1);
for ii = 1:obj.n_out
    % Calculate unmeasured disturbance bias for all outputs
    bias = y(ii) - obj.F((ii-1)*obj.n+1);
    
    % Calculate free response for all outputs
    TF = obj.T * obj.F((ii-1)*obj.n+1:ii*obj.n);
    
    % Calculate predicted disturbance
    if obj.n_dis > 0
        dist = reshape(permute(obj.Sd(ii,:,1:obj.p),...
                               [3,1,2]),...
                       obj.p,obj.n_dis,1) * ...
               (d - obj.d_past);
    end
    
    % Calculate outputs due to past and present terms
    if obj.n_dis > 0
        yp((ii-1)*obj.p+1:ii*obj.p) = TF + dist + bias;
    else
        yp((ii-1)*obj.p+1:ii*obj.p) = TF + bias;
    end
end

% Predicted deviation from setpoint ------------------------------------- %
% Bring setpoint to needed dimensions [n_out*p by 1]
if size(y_sp,2) == 1
    % Repeat constant setpoint for prediction horizon
    y_sp = repmat(y_sp,1,obj.p);
end
y_sp = reshape(y_sp(:,1:obj.p)',[],1);

er = y_sp - yp;

% Gradient -------------------------------------------------------------- f
f = - obj.G' * obj.Wy * er;

% Vector for constraint handling b -------------------------------------- b
% Constraints: A * [Delta_u; epsilon] <= b
b = [];
% User upper and lower bounds instead
% if ~isempty(obj.D_u_max)
%     b1max = reshape(repmat(obj.D_u_max,1,obj.c)',[],1);
%     b = [b; b1max];
% end
% if ~isempty(obj.D_u_min)
%     b1min = reshape(repmat(-obj.D_u_min,1,obj.c)',[],1);
%     b = [b; b1min];
% end

if ~isempty(obj.u_max)
    b2max = reshape(repmat(obj.u_max,1,obj.c)',[],1) - ...
            reshape(repmat(obj.u_past,1,obj.c)',[],1);
    b = [b; b2max];
end
if ~isempty(obj.u_min)
    b2min = reshape(repmat(-obj.u_min,1,obj.c)',[],1) - ...
            reshape(repmat(-obj.u_past,1,obj.c)',[],1);
    b = [b; b2min];
end

if ~isempty(obj.y_max)
    b3max = reshape(repmat(obj.y_max,1,obj.p)',[],1) - yp;
    b = [b; b3max];
end
if ~isempty(obj.y_min)
    b3min = -reshape(repmat(obj.y_min,1,obj.p)',[],1) + yp; ...
    b = [b; b3min];
end

% Future control moves -------------------------------------------------- %
if isempty(obj.A) % If unconstrained
    % Closed form solution for control moves
    x_opt = - obj.H \ f;
else % If constrained
    % Adjust H and f for soft constraints
    if sum(sum(obj.soft)) > 0
        hsym = [obj.Hsym sparse(size(obj.Hsym,1),size(obj.W_eps,2)); ...
                sparse(size(obj.W_eps,1),size(obj.Hsym,2)) obj.W_eps];
        f_soft = [f; 0.5 * obj.w_eps];
    else
        hsym = obj.Hsym;
        f_soft = f;
    end
    % Quadratic program
    [x_opt, ~, exitflag] = quadprog(hsym,f_soft,obj.A,b,[],[],obj.lb,obj.ub,[],obj.ops);
    
    % If unsuccessfull, try to soften ALL output constraints
    if exitflag ~= 1 && sum(sum(obj.soft)) < numel(obj.soft) &&...
            (~isempty(obj.y_min) || ~isempty(obj.y_max))
        warning(['Quadprog did not converge, exitflag: ',num2str(exitflag),'. Soften all output constraints'])
        % Fill in zeros in w_eps and W_eps with default value 100*Q
        q100 = [100 * diag(obj.Wy); 100 * diag(obj.Wy)];
        w_eps_full = obj.w_eps;
        if isempty(w_eps_full)
            w_eps_full = zeros(2*obj.n_out*obj.p,1);
        elseif length(w_eps_full) < length(q100)
            w_eps_full = [w_eps_full; w_eps_full];
        end
        w_eps_full(w_eps_full==0) = q100(w_eps_full==0);
        W_eps_full = diag(obj.W_eps);
        if isempty(W_eps_full)
            W_eps_full = zeros(2*obj.n_out*obj.p,1);
        elseif length(W_eps_full) < length(q100)
            W_eps_full = [W_eps_full; W_eps_full];
        end
        W_eps_full(W_eps_full==0) = q100(W_eps_full==0);
        W_eps_full = spdiags(W_eps_full,0,length(W_eps_full),length(W_eps_full));
        % Define Constraint matrix and bounds
        % depends on if there are ymax and ymin and umax and umin
        counter_u = 0;
        if ~isempty(obj.u_max)
            counter_u = counter_u + 1;
        end
        if ~isempty(obj.u_min)
            counter_u = counter_u + 1;
        end
        counter_y = 0;
        if ~isempty(obj.y_max)
            counter_y = counter_y + 1;
        end
        if ~isempty(obj.y_min)
            counter_y = counter_y + 1;
        end
        a = [obj.A_u ...
             [sparse(counter_u*obj.n_in*obj.c,counter_y*obj.n_out*obj.p); ...
          -speye(counter_y*obj.n_out*obj.p)]];
        Lb = [obj.lb_u; sparse(counter_y*obj.n_out*obj.p,1)];
        Ub = [obj.ub_u; inf(counter_y*obj.n_out*obj.p,1)];
        % Adjust H and f for soft constraints
        hsym = [obj.Hsym sparse(size(obj.Hsym,1),size(W_eps_full,2)); ...
            sparse(size(W_eps_full,1),size(obj.Hsym,2)) W_eps_full];
        f_soft = [f; 0.5 * w_eps_full];
        % Quadratic program
        [x_opt, ~, exitflag] = quadprog(hsym,f_soft,a,b,[],[],Lb,Ub,[],obj.ops);
    end
    
    % If unsuccessfull, use unconstrained solution
    if exitflag ~= 1
        warning(['Quadprog did not converge, exitflag: ',num2str(exitflag),' Unconstrained solution is used'])
        Delta_u = - obj.H \ f;
        % Make sure solutions are within bounds of u
        u_unc = Delta_u(1:obj.c:obj.n_in*obj.c) + obj.u_past; % From Delta_u to u
        if ~isempty(obj.u_max)
            u_unc = min([u_unc obj.u_max],[],2);    % smaller than u_max
        end
        if ~isempty(obj.u_min)  
            u_unc = max([u_unc obj.u_min],[],2);    % larger than u_min
        end
        Delta_u(1:obj.c:obj.n_in*obj.c) = u_unc - obj.u_past; % From u to Delta_u
        % Fill in zeros for epsilon
        x_opt = [Delta_u; zeros(length(obj.w_eps),1)];
    end
end

% First move to implement ----------------------------------------------- u
% Extract Delta_u and slack variable epsilon
Delta_u = x_opt(1:obj.n_in*obj.c);
if sum(sum(obj.soft)) > 0 || exist('w_eps_full','var')
    epsilon = x_opt(obj.n_in*obj.c+1:end);
else
    epsilon = zeros(2*obj.n_out*obj.p,1);
end

u = Delta_u(1:obj.c:obj.n_in*obj.c) + obj.u_past;
% Possibly add the possibility to round u to the desired amount of
% significant figures and use the rounded u for the free response update
% Alternatively, move the free response update to the beginning of the
% function, using the user provided u_past



% Calculate value of objective function (optional) ---------------------- J
if nargout == 2
    % Output part
    J1 = (yp + obj.G * Delta_u(1:obj.n_in*obj.c) - y_sp)' * obj.Wy * ...
         (yp + obj.G * Delta_u(1:obj.n_in*obj.c) - y_sp);
    % Input part
    J2 = Delta_u(1:obj.n_in*obj.c)' * obj.Wu * Delta_u(1:obj.n_in*obj.c);
    % Linear slack variable part
    if exist('w_eps_full','var')
        J3 = w_eps_full' * epsilon;
    elseif sum(sum(obj.soft)) > 0
        J3 = obj.w_eps' * epsilon;
    else
        J3 = 0;
    end
    
    % Quadratic slck variable part
    if exist('W_eps_full','var')
        J4 = epsilon' * W_eps_full * epsilon;
    elseif sum(sum(obj.soft)) > 0
        J4 = epsilon' * obj.W_eps * epsilon;
    else
        J4 = 0;
    end
    J = [J1; J2; J3; J4];
end


% Update free response for next iteration ------------------------------- %
for ii = 1:obj.n_out
    % Change index of F
    MF = obj.M * obj.F((ii-1)*obj.n+1:ii*obj.n);
    
    % Update last element with past input for integrating input
    lu = [zeros(obj.n-1,1); obj.hu(ii,:) * obj.u_past];
    
    % Update last element with past disturbance for integrating disturbance
    if obj.n_dis > 0
        ld = [zeros(obj.n-1,1); obj.hd(ii,:) * obj.d_past];
    end
    
    % Update all elements with current input change
    cu = permute(obj.Su(ii,:,:),[3,2,1]) * (u - obj.u_past);
    
    % Update all elements with current disturbance change
    if obj.n_dis > 0
        cd = permute(obj.Sd(ii,:,:),[3,2,1]) * (d - obj.d_past);
    end
    
    % Update free response
    if obj.n_dis > 0
        obj.F((ii-1)*obj.n+1:ii*obj.n) = MF + lu + ld + cu + cd;
    else
        obj.F((ii-1)*obj.n+1:ii*obj.n) = MF + lu + cu;
    end
    
end

% update u_past and d_past ---------------------------------------------- %
obj.u_past = u;
obj.d_past = d;

end % End of function QDMC_run


end % End of methods


end % End of classdef