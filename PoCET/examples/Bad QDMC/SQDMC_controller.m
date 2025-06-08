classdef SQDMC_controller < handle
% Stochastic Quadratic Dynamic Matrix Control for MIMO System
%
% This function creates and runs a stochastic quadratic dynamic matrix
% control algorithm for MIMO systems. See the note on "Model Predictive
% Control and State Estimation" by Enso Ikonen (2013) and the work on SQDMC
% by Joel Paulson. It determines the input changes based on the plant model
% (step response), current measurement, and current disturbance deviations
% w.r.t. the preceeding step.
%
% Copied in great parts from Joel Paulson, 1-23-2014
%
% Revised and improved by Matthias von Andrian 01/2018
%
% Transfered to classdef mode by Matthias von Andrian 01/2018


properties (GetAccess = public, SetAccess = public)
% Stored parameters, assigned in the initialization function
% SQDMC_controller and used in the run function run_controller.
% With full access outside of class
    u_past      % Control input implemented in the previous time step,
                % vector [n_in by 1]
    d_past      % Measured disturbance in the previous time step,
                % vector [n_dis by 1]
    ops         % Options for quadratic program solver,
                % optimoptions('quadprog',...)
end % End of public properties

properties (GetAccess = public, SetAccess = protected)
% Stored parameters, assigned in the initialization function
% SQDMC_controller and used in the run function run_controller.
% With no access outside of class
	n_in        % Number of input variables, scalar
    n_out       % Number of output variables, scalar
    n_dis       % Number of measured disturbances, scalar
    n           % Number of step response coefficients, scalar
    L           % Number of terms in the ploynomial chaos expansion (PCE)
    IP          % Inner product of basis function of the PCE with themselves,
                % IP(j) = <phi(j-1),phi(j-1)> since the index starts from 1 
                % while the nomenclature starts from 0, vector [L by 1]
    Su          % Step response coefficients input, matrix,
                % [n_out by L by n_in by n] -> 4D matrix
    Sd          % Step response coefficients disturbances, matrix,
                % [n_out by L by n_dis by n] -> 4D matrix
    hu          % Slope between the last two input step response
                % coefficients (assumed to be constant, assumed to be 0 for
                % nonintegrating pairings) matrix [n_out by L by n_in]
    hd          % Slope between the last two disturbance step response
                % coefficients (assumed to be constant, assumed to be 0 for
                % nonintegrating pairings) matrix [n_out by L by n_dis]
    p           % Prediction horizon, scalar
    c           % Control horizon, scalar
    u_min       % Minimum input, vector [n_in by 1]
    u_max       % Maximum input, vector [n_in by 1]
    D_u_min     % Minimum input change, vector [n_in by 1]
    D_u_max     % Maximum input change, vector [n_in by 1]
    y_min       % Minimum output, vector [n_out by 1]
    y_max       % Maximum output, vector [n_out by 1]
    T           % Matrix to extract free response within prediction horizon
                % binary matrix [p by n]
    G           % Dynamic matrix, combination of lower triangular matrices
                % with elements of Su, matrix [p*n_out by c*n_in by L]
    % G in cell form
    F           % Free response -> Output, if no further controler action
                % was taken, vector [n*n_out by L]
    Wu          % Diagonal representation of the control weight parameters
    Wy          % Diagonal representation of the output weight parameters
    Wvar        % Diagonal representation of the variance weight parameters
    H           % Hessian of the quadratic program
    Hsym        % Symmetric Hessian of the quadratic program
    M           % Matrix to update free response, binary matrix [n by n]
    A           % Matrix to handle inequality constraints A * Delta_u <= b
end % End of protected properties


methods
    
% ======================================================================= %
% Initialization                                           Initialization %
% ======================================================================= %
function obj = SQDMC_controller(n_in,n_out,L,IP,n_dis,Su,Sd,p,c,La,Q,ctg,Qvar,...
                               u_min,u_max,D_u_min,D_u_max,y_min,y_max,...
                               u_past,y_past,d_past,u_int,d_int)
% This function checks for consistant inputs and brings them into the
% required format and creates structures needed to run the SQDMC controller.
%
% Inputs (* marks optional input, leave empty [] if not applicable)
%   n_in        Number of input variables, scalar
%   n_out       Number of output variables, scalar
%   L           Number of terms in the ploynomial chaos expansion (PCE),
%               scalar
%   IP          Inner product of basis function of the PCE with themselves,
%               IP(j) = <phi(j-1),phi(j-1)> since the index starts from 1 
%               while the nomenclature starts from 0, vector [L by 1]
% * n_dis *     Number of measured disturbances, scalar
%   Su          Step response coefficients input, matrix,
%               [n_out by L by n_in by n] -> 4D matrix
% * Sd *        Step response coefficients disturbances, matrix,
%               [n_out by L by n_dis by n] -> 4D matrix
%   p           Prediction horizon, scalar
%   c           Control horizon, scalar
% * La *        Vector of control weight parameters, 
%               i.e., J = ||r-y|| + La*eye*||du|| where du is the input
%               change, [n_in by 1]
% * Q *         Vector of output weight parameters, vector [n_out by 1]
% * ctg *       Vector of cost to go parameters (additional weight for the
%               last predicted output error, vector [n_out by 1]
% * Qvar *      Vector of variance weight parameters, vector [n_out by 1]
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
%               Binary matrix [p.n_out by p.n_in]
%               If not supplied, it will be set to 0
% * d_int *     Is disturbance integrating on specific output?
%               Binary matrix [p.n_out by p.n_dis]
%               If not supplied, it will be set to 0

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
if isempty(L) || ~isscalar(L)
    error(['Please specify the number of terms in the ploynomial ',...
           'chaos expansion L as scalar'])
end
if ~isempty(IP) && (size(IP,1) ~= L || size(IP,2) ~= 1)
    error('IP has to be a column vector with L entries')
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
    error('y_min has to be a column vector with n_out entries') %%%%%%%%%%%%%%%% or n_out by p if constraints change over prediciton horizon
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
    error('y_past has to be a vector with dimensions [n_out by 1]')
end
if ~isempty(y_past) && (size(y_past,1) == n_out || size(y_past,2) == 1)
    y_past = [y_past zeros(n_out,L-1)]; % Assume steady state for all higher PCE coefficients
end
if isempty(y_past)
    y_past = zeros(n_out,L); % Assume steady state for all PCE coefficients
end
if ~isempty(d_past) && (size(d_past,1) ~= n_dis || size(d_past,2) ~= 1)
    error('d_past has to be a column vector with n_dis entries')
end
if isempty(d_past) && n_dis > 0
    d_past = zeros(n_dis,1); % Assume steady state
end

% Integrating inputs ---------------------------------------------------- %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% automate this with a cutoff value, if difference between last values is
% larger than cutoff (can be a relative value), mark as integrating.

% -> Inputs that give a stepresponse with constant nonzero slope at n
if ~isempty(u_int) && ...
   (size(u_int,1) ~= n_out || size(u_int,2) ~= L || size(u_int,3) ~= n_in)
    error('u_int has to be a matrix with dimensions [n_out by L by n_in]')
end
if isempty(u_int)
    u_int = zeros(n_out,L,n_in); % Assume non integrating inputs
end
if ~isempty(d_int) && ...
   (size(d_int,1) ~= n_out || size(d_int,2) ~= L || size(d_int,3) ~= n_dis)
    error('d_int has to be a matrix with dimensions [n_out by L by n_dis]')
end
if isempty(u_int)
    d_int = zeros(n_out,L,n_dis);
end

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
if ~isempty(Qvar) && length(Qvar(:)) ~= n_out
    error('Qvar has to be a vector with n_out entries')
end
if isempty(Qvar)
    Qvar = ones(n_out,1);
end

% Step response matrices ------------------------------------------------ %
[d1, ~, d3, ~] = size(Su);
if d1 ~= n_out || d3 ~= n_in
    error(['Step response coefficients for input Su have to be in a ',...
           'matrix with dimensions [n_out by L by n_in by n]'])
end
if n_dis > 0
    [d1, ~, d3, ~] = size(Sd);
    if d1 ~= n_out || d3 ~= n_dis
        error(['Step response coefficients for measured disturbances ',...
               'Sd have to be in a matrix with dimensions ',...
               '[n_out by L by n_dis by n]'])
    end
else
    Sd = [];
end
% Make sure, number of step response coefficients is equal for input
% and disturbance
if n_dis > 0
    [~, ~, ~, d4u] = size(Su);
    [~, ~, ~, d4d] = size(Sd);
    if d4u > d4d
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fill in Su so it becomes the same length as Sd
    elseif d4u < d4d
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fill in Sd so it becomes the same length as Su
    end
end
% Extend step response coefficients, if shorter than prediction horizon
n = size(Su,4);
if n <= p
    % Copy last value for non integrating, keep constant slope for
    % integrating
    % Input - output
    hu = (Su(:,:,:,n) - Su(:,:,:,n-1)) .* u_int;
    Su_extended = repmat(Su(:,:,:,n),1,1,1,p-n+1) + ...
        (repmat(hu,1,1,1,p-n+1) .* ...
         repmat(permute(1:p-n+1,[1,4,3,2]),n_out,L,n_in,1));
    Su = cat(4,Su,Su_extended);
    % Disturbance - output
    if n_dis > 0
        hd = (Sd(:,:,:,n) - Sd(:,:,:,n-1)) .* d_int;
        Sd_extended = repmat(Sd(:,:,:,n),1,1,1,p-n+1) + ...
            (repmat(hd,1,1,1,p-n+1) .* ...
             repmat(permute(1:p-n+1,[1,4,3,2]),n_out,L,n_dis,1));
        Sd = cat(4,Sd,Sd_extended);
    end
    n = p + 1;
end


% Set parameters (properties) =========================================== %
% Number of variables --------------------------------------------------- n
obj.n_in  = n_in;
obj.n_out = n_out;
obj.L = L;
obj.IP = IP;
obj.n_dis = n_dis;
obj.n = n;

% Step response matrices -------------------------------------------- Su Sd
obj.Su = Su;
obj.Sd = Sd;

% Difference between last step response coefficients ---------------- hu hd
obj.hu = (Su(:,:,:,n) - Su(:,:,:,n-1)) .* u_int;
if n_dis > 0
    obj.hd = (Sd(:,:,:,n) - Sd(:,:,:,n-1)) .* d_int;
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

% Past measurements ----------------------------------------- u_past d_past
obj.u_past = u_past;
obj.d_past = d_past;

% Determine matrix T ---------------------------------------------------- T
obj.T = ...
    sparse(p,n) + [sparse(p,1) spdiags(ones(p,1),0,p,p) sparse(p,n-p-1)];

% Determine matrix G ---------------------------------------------------- G
obj.G = zeros(p*n_out,c*n_in,L);
for ii = 1:n_out
    for jj = 1:n_in
        for kk = 1:L
            Su_part = permute(Su(ii,kk,jj,:),[4,1,2,3]);
            obj.G((ii-1)*p+1:ii*p, (jj-1)*c+1:jj*c,kk) = ...
                toeplitz(Su_part(1:p),Su_part(1,1)*eye(1,c));
        end
    end
end

% Determine initial free response F ------------------------------------- F
% Assume process starts from a steady state
obj.F = zeros(n_out*n,L);
for ii = 1:n_out
    obj.F((ii-1)*n+1:ii*n,:) = repmat(y_past(ii,:),n,1);
end

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

% Determine variance weighting matrix Wvar --------------------------- Wvar
obj.Wvar = spdiags(...
            reshape(...
             repmat(Qvar,1,p)',...
            n_out*p,1),...
           0,n_out*p,n_out*p);

% Calculate Hessian H --------------------------------------------------- H
obj.H = obj.G(:,:,1)' * obj.Wy * obj.G(:,:,1) + obj.Wu;
for ii = 2:L
    obj.H = obj.H + obj.G(:,:,ii)' * obj.Wvar * obj.G(:,:,ii) * obj.IP(ii);
end
obj.Hsym = (obj.H+obj.H')/2;  % Quadprog requires symmetric hessian
    
% Determine matrix M to update free response ---------------------------- M
obj.M = spdiags(ones(n,1),1,n,n); obj.M(n,n) = 1;

% Set options for quadprog ops ---------------------------------------- ops
obj.ops = optimoptions('Quadprog'...
                      ,'Display','off'...
                      ,'Algorithm','interior-point-convex'...
                      ,'MaxIterations',1e4);

% Set matrix for constraint handling A ---------------------------------- A
% Constraints: A * Delta_u <= b
% b is dependent on future delta u and is determined in run_controller
A = [];
if ~isempty(obj.D_u_max)
    A1max = eye(obj.n_in*obj.c);
    A = [A; A1max];
end
if ~isempty(obj.D_u_min)
    A1min = -eye(obj.n_in*obj.c);
    A = [A; A1min];
end

if ~isempty(obj.u_max)
    A2 = tril(ones(obj.c));
    A2max = [];
    for jj = 1:obj.n_in
        A2max = blkdiag(A2max,A2);
    end
    A = [A; A2max];
end
if ~isempty(obj.u_min)
    A2 = tril(ones(obj.c));
    A2min = [];
    for jj = 1:obj.n_in
        A2min = blkdiag(A2min,-A2);
    end        
    A = [A; A2min];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% can include tuning factor alpha here: y_min <= E[y] - alpha * Var[y]
if ~isempty(obj.y_max)
    A3max = obj.G(:,:,1);
    A = [A; A3max];
end
if ~isempty(obj.y_min)
    A3min = -obj.G(:,:,1);
    A = [A; A3min];
end
obj.A = A;

end % End of initialization function SQDMC_controller


% ======================================================================= %
% update_constraints                                   update_constraints %
% ======================================================================= %
function obj = update_constraints(obj,u_min,u_max,D_u_min,D_u_max,y_min,y_max)
% Update the matrix A to handle the constraints A * Delta_u <= b

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

% Set matrix for constraint handling A
a = []; % Lower case to avoid warning about obj.A and A
        % (although, here they are really the same)
if ~isempty(obj.D_u_max)
    A1max = eye(obj.n_in*obj.c);
    a = [a; A1max];
end
if ~isempty(obj.D_u_min)
    A1min = -eye(obj.n_in*obj.c);
    a = [a; A1min];
end

if ~isempty(obj.u_max)
    A2 = tril(ones(obj.c));
    A2max = [];
    for jj = 1:obj.n_in
        A2max = blkdiag(A2max,A2);
    end
    a = [a; A2max];
end
if ~isempty(obj.u_min)
    A2 = tril(ones(obj.c));
    A2min = [];
    for jj = 1:obj.n_in
        A2min = blkdiag(A2min,-A2);
    end        
    a = [a; A2min];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% can include tuning factor alpha here: y_min <= E[y] - alpha * Var[y]
if ~isempty(obj.y_max)
    A3max = obj.G(:,:,1);
    a = [a; A3max];
end
if ~isempty(obj.y_min)
    A3min = -obj.G(:,:,1);
    a = [a; A3min];
end
obj.A = a;

end % End of function update_constraints



% ======================================================================= %
% update_weights                                           update_weights %
% ======================================================================= %
function obj = update_weights(obj,La,Q,ctg,Qvar)
% Update the weight matrices

% Check input
if ~isempty(La) && (size(La,1) ~= obj.n_in || size(La,2) ~= 1)
    error('La has to be a column vector with n_in entries')
end
if isempty(La)
    La = ones(obj.n_in,1);
end
if ~isempty(Q) && length(Q(:)) ~= obj.n_out
    error('Q has to be a vector with n_out entries')
end
if isempty(Q)
    Q = ones(obj.n_out,1);
end
if ~isempty(ctg) && length(ctg(:)) ~= obj.n_out
    error('ctg has to be a vector with n_out entries')
end
if isempty(ctg)
    ctg = zeros(obj.n_out,1);
end
if ~isempty(Qvar) && length(Qvar(:)) ~= obj.n_out
    error('Qvar has to be a vector with n_out entries')
end
if isempty(Qvar)
    Qvar = ones(obj.n_out,1);
end

% Update weight matrices
% Determine input weighting matrix Wu ---------------------------------- Wu
obj.Wu = spdiags(...
          reshape(...
           repmat(La,1,obj.c)',...
          obj.n_in*obj.c,1),...
         0,obj.n_in*obj.c,obj.n_in*obj.c);

% Determine output weighting matrix Wy --------------------------------- Wy
obj.Wy = spdiags(...
          reshape(...
           repmat(Q,1,obj.p)' + [zeros(obj.p-1,obj.n_out);ctg'],...
          obj.n_out*obj.p,1),...
         0,obj.n_out*obj.p,obj.n_out*obj.p);

% Determine variance weighting matrix Wvar --------------------------- Wvar
obj.Wvar = spdiags(...
            reshape(...
             repmat(Qvar,1,obj.p)',...
            obj.n_out*obj.p,1),...
           0,obj.n_out*obj.p,obj.n_out*obj.p);

% Calculate Hessian H --------------------------------------------------- H
obj.H = obj.G(:,:,1)' * obj.Wy * obj.G(:,:,1) + obj.Wu;
for ii = 2:obj.L
    obj.H = obj.H + obj.G(:,:,ii)' * obj.Wvar * obj.G(:,:,ii) * obj.IP(ii);
end
obj.Hsym = (obj.H+obj.H')/2;  % Quadprog requires symmetric hessian


end % End of function update_weights



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
yp = zeros(obj.n_out*obj.p,obj.L);
for ii = 1:obj.n_out
    for jj = 1:obj.L
        % Calculate free response for all outputs
        TF = obj.T * obj.F((ii-1)*obj.n+1:ii*obj.n,jj);
    
        % Calculate predicted effect of disturbance
        if obj.n_dis > 0
            dist = reshape(permute(obj.Sd(ii,jj,:,1:obj.p),...
                                   [4,1,3,2]),...
                           obj.p,obj.n_dis,1) * ...
                   (d - obj.d_past);
        end
    
        % Calculate outputs due to past and present terms
        if obj.n_dis > 0
            yp((ii-1)*obj.p+1:ii*obj.p,jj) = TF + dist;
        else
            yp((ii-1)*obj.p+1:ii*obj.p,jj) = TF;
        end
    end
    
    % Calculate unmeasured disturbance bias for all outputs
    bias = y(ii) - obj.F((ii-1)*obj.n+1,1);
    
    % Calculate outputs due to past and present terms
    % (add bias to 0th PCE coefficient only)
    yp((ii-1)*obj.p+1:ii*obj.p,1) = yp((ii-1)*obj.p+1:ii*obj.p,1) + bias;
end

% Predicted expected deviation from setpoint ---------------------------- %
% Bring setpoint to needed dimensions [n_out*p by 1]
if size(y_sp,2) == 1
    % Repeat constant setpoint for prediction horizon
    y_sp = repmat(y_sp,1,obj.p);
end
y_sp = reshape(y_sp(:,1:obj.p)',[],1);

er = y_sp - yp(:,1);

% Gradient -------------------------------------------------------------- f
f = - obj.G(:,:,1)' * obj.Wy * er;
for ii = 2:obj.L
    f = f + obj.G(:,:,ii)' * obj.Wvar * yp(:,ii) * obj.IP(ii);
end

% Vector for constraint handling b -------------------------------------- b
% Constraints: A * Delta_u <= b
b = [];
if ~isempty(obj.D_u_max)
    b1max = reshape(repmat(obj.D_u_max,1,obj.c)',[],1);
    b = [b; b1max];
end
if ~isempty(obj.D_u_min)
    b1min = reshape(repmat(-obj.D_u_min,1,obj.c)',[],1);
    b = [b; b1min];
end

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% can include tuning factor alpha here: y_min <= E[y] - alpha * Var[y]
if ~isempty(obj.y_max)
    b3max = reshape(repmat(obj.y_max,1,obj.p)',[],1) - yp(:,1);
    b = [b; b3max];
end
if ~isempty(obj.y_min)
    b3min = -reshape(repmat(obj.y_min,1,obj.p)',[],1) + yp(:,1); ...
    b = [b; b3min];
end

% Future control moves -------------------------------------------------- %
if isempty(obj.A) % If unconstrained
    % Closed form solution for control moves
    Delta_u = - obj.H \ f;
else % If constrained
    % Quadratic program
    [Delta_u, ~, exitflag] = quadprog(obj.Hsym,f,obj.A,b,[],[],[],[],[],obj.ops);
    if exitflag ~= 1
        warning(['Quadprog did not converge, exitflag: ',num2str(exitflag),' Unconstrained solution is used'])
        Delta_u = - obj.H \ f;
    end
end

% First move to implement ----------------------------------------------- u
u = Delta_u(1:obj.c:length(Delta_u)) + obj.u_past;
% Possibly add the possibility to round u to the desired amount of
% significant figures and use the rounded u for the free response update
% Alternatively, move the free response update to the beginning of the
% function, using the user provided u_past

% Calculate value of objective function (optional) ---------------------- J
if nargout == 2
    % Expected value part
    J1 = (yp(:,1) + obj.G(:,:,1) * Delta_u - y_sp)' * obj.Wy * ...
         (yp(:,1) + obj.G(:,:,1) * Delta_u - y_sp);
    % Variance part
    Var = zeros(obj.n_out,1);
    for ii = 2:obj.L
        Var = Var + obj.IP(ii) * ...
              (yp(:,ii) + obj.G(:,:,ii) * Delta_u)' * ...
              (yp(:,ii) + obj.G(:,:,ii) * Delta_u);
    end
    Wvar_out = diag(obj.Wvar); Wvar_out = full(Wvar_out(1:obj.p:end));
%     Wvar_out = 1;
    J2 = Wvar_out' * Var;
    % Input part
    J3 = Delta_u' * obj.Wu * Delta_u;
    J = [J1; J2; J3];
end

% Update free response for next iteration ------------------------------- %
for ii = 1:obj.n_out
    for jj = 1:obj.L
        % Change index of F
        MF = obj.M * obj.F((ii-1)*obj.n+1:ii*obj.n,jj);

        % Update last element with past input for integrating input
        lu = [zeros(obj.n-1,1); permute(obj.hu(ii,jj,:),[1,3,2]) * obj.u_past];

        % Update last element with past disturbance for integrating disturbance
        if obj.n_dis > 0
            ld = [zeros(obj.n-1,1); permute(obj.hd(ii,jj,:),[1,3,2]) * obj.d_past];
        end

        % Update all elements with current input change
        cu = permute(obj.Su(ii,jj,:,:),[4,3,2,1]) * (u - obj.u_past);

        % Update all elements with current disturbance change
        if obj.n_dis > 0
            cd = permute(obj.Sd(ii,jj,:,:),[4,3,2,1]) * (d - obj.d_past);
        end

        % Update free response
        if obj.n_dis > 0
            obj.F((ii-1)*obj.n+1:ii*obj.n,jj) = MF + lu + ld + cu + cd;
        else
            obj.F((ii-1)*obj.n+1:ii*obj.n,jj) = MF + lu + cu;
        end
        % Keep variance at 0
%         obj.F(:,2:obj.L) = zeros(obj.n*obj.n_out,obj.L-1);
    end
    
end

% update u_past and d_past ---------------------------------------------- %
obj.u_past = u;
obj.d_past = d;

end % End of function run_controller


end % End of methods


end % End of classdef