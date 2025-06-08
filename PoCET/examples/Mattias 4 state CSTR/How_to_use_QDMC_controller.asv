% This file shows how to use the file QDMC_controller_soft.m,
% a quadratic dynamic matrix controller with soft constraints.
%
% Both files created by Matthias von Andrian, andrian@mit.edu
% a-n-d-r-i-a-n[at]m-i-t.e-d-u  -> remove the "-" signs
% For questions and to report bugs, please email me.
%
% For accademic use, the rules of properly citing other people's work apply

%% Clean up ============================================================= %
clc; clear all; close all;

%% Define parameters and functions ====================================== %
% Parameters to create step response
dt = 0.1;       % Delta t (has to be identical when creating the step 
                % response model and when running closed loop control)
tsp_sr = 0:dt:10;   % Time span to create the step response
n_sr = length(tsp_sr) - 1;  % Length of the step response model
ops = odeset;   % ODE options

% Parameters to design QDMC (more details, including required dimensions
% can be found in the file QDMC_controller_soft.m under initialization)
n_in = 1;       % Number of inputs
n_out = 1;      % Number of outputs
n_dis = 0;      % Number of measured disturbances

Sd = [];        % Step resposne of measured disturbances
p = 100;        % Prediction horizon (number of time steps)
c = 50;         % Control horizon (number of time steps)

La = 1;         % Weight for input movement
Q = 25;         % Weight for output error
ctg = [];       % Cost to go

u_past = 0;     % Past input (assume process starts at steady state)
y_past = 0;     % Past output (assume process starts at steady state)
d_past = [];    % Past measured disturbance
u_int = [];     % Integrating inputs
d_int = [];     % Integrating measured disturbances

u_min = -1;     % Minimum input
u_max = 1.5;    % Maximum input
D_u_min = -0.1; % Minimum input movement
D_u_max = 0.1;  % Maximum input movement
y_min = -0.1;   % Minimum output
y_max = 1.1;    % Maximum output

soft = [1,1];   % Both, y_max and y_min are soft constraints
w_eps = 10;     % Weight for linear constraint violation
W_eps = 10;     % Weight for quadratic constraint violation

% Parameters to run closed loop simulation
ysp = [ones(1,5/dt), zeros(1,5/dt), zeros(1,10/dt)]; % Output setpoints
% Setpoint of 1 for the first 5 time units, then 0
udist = [zeros(1,10/dt), 0.2.*ones(1,5/dt), zeros(1,5/dt)]; 
% Input disturbance (unmeasured), for time unit 10 to 15
tsp_sim = 0:dt:20;  % Time to simulate the closed loop

% Functions
dxdt = @(t,x,u) -0.5 * x + u;   % This is the SISO state space system 
                                % considered in this example
n_x = 1;        % Number of states


%% Create step response model =========================================== %
% In this section, a step test of the input is performed and the
% resulting system answer is recorded and stored as a step response model

% This is shown in a more general setting, for the SISO case the for-loops
% are not needed

% The step response Sd for the measured disturbances is created in the same
% way

x0 = zeros(n_x,1);              % Start the process at steady state,
                                % expressed as deviation variable
Su = zeros(n_out,n_in,n_sr);    % Preallocation 
                                % of step resposne coefficient matrix
for ii = 1:n_in
    % Run step test on system
    u_step = zeros(n_in,1);
    u_step(ii) = 1;     % Perform unit step for each input
    [~,X] = ode45(@(t,x) dxdt(t,x,u_step),tsp_sr,x0,ops);
    
    % Record step response (do not record response for t=0)
    Su(:,ii,:) = permute(X(2:end,:),[2,3,1]);
end

% Plot step response model
figure
for ii = 1:n_in
    for jj = 1:n_out
        subplot(n_in,n_out,(ii-1)*n_in + jj)
        plot(tsp_sr(2:end),permute(Su(jj,ii,:),[3,2,1]),'*')
        xlabel('Time')
        ylabel(['y_',num2str(jj)])
        title(['Step of u_',num2str(ii)])
    end
end


%% Perform closed loop control ========================================== %
% In this section, a setpoint change and an input disturbance rejection is
% simulated

% Create instance of the controller
Controller = QDMC_controller_soft(n_in,n_out,n_dis,Su,Sd,p,c,La,Q,ctg,...
                               u_min,u_max,D_u_min,D_u_max,y_min,y_max,...
                               u_past,y_past,d_past,...
                               u_int,d_int,soft,w_eps,W_eps);

% Preallocate storage for inputs and outputs
x_store = zeros(n_x,length(tsp_sim));       % States are measured at t
y_store = zeros(n_out,length(tsp_sim));     % Outputs are measured at t
u_store = zeros(n_in,length(tsp_sim)-1);    % Inputs are applied over delta t
% Define states and outputs for t=0
x = zeros(n_x,1);   % Start the process from steady state (deviation variable)
y = zeros(n_out,1); % Start the process from steady state (deviation variable)
% Write value for t=0 to storage
x_store(:,1) = x;
y_store(:,1) = y;

% Run closed loop
for ii = 1:length(tsp_sim)-1
    % Current setpoint
    ysp_curr = ysp(:,ii);
    % Current (unmeasured) input disturbance
    udist_curr = udist(ii);
    % Current measured disturbance
    d = [];
    
    % Calculate controller output = process input
    u = Controller.run_controller(y,ysp_curr,d);
    
    % Add unmeasured input disturbance to controller signal
    u_process = u + udist_curr;
    
    % Simulate one time step
    [~,X] = ode45(@(t,x) dxdt(t,x,u_process),[0 dt],x,ops);
    x = X(end,:)';  % State at the end of the time step
    y = x;  % Process output
    
    % Store results
    x_store(:,ii+1) = x;
    y_store(:,ii+1) = y;
    u_store(:,ii) = u;
end

%% Plot results ========================================================= %
% Outputs
figure
for ii = 1:n_out
    subplot(n_out,1,ii)
    hold on
    % Plot constraints
    if ~isempty(y_min)
        plot([tsp_sim(1);tsp_sim(end)],[y_min(ii);y_min(ii)],'k')
    end
    if ~isempty(y_max)
        plot([tsp_sim(1);tsp_sim(end)],[y_max(ii);y_max(ii)],'k')
    end
    plot(tsp_sim(1:end-1),ysp(ii,:),'k--') % Plot setpoint
    plot(tsp_sim,y_store(ii,:)) % Plot closed loop output
    grid on
    xlabel('Time')
    ylabel(['y_',num2str(ii)])
end

% States
figure
for ii = 1:n_x
    subplot(n_x,1,ii)
    plot(tsp_sim,x_store(ii,:))
    grid on
    xlabel('Time')
    ylabel(['x_',num2str(ii)])
end

% Inputs
figure
for ii = 1:n_in
    subplot(n_in,1,ii)
    hold on
    % Plot constraints
    if ~isempty(u_min)
        plot([tsp_sim(1);tsp_sim(end)],[u_min(ii);u_min(ii)],'k')
    end
    if ~isempty(u_max)
        plot([tsp_sim(1);tsp_sim(end)],[u_max(ii);u_max(ii)],'k')
    end
    stairs(tsp_sim,[u_store(ii,:),u_store(ii,end)]) % repeat last element to finish stair plot
    grid on
    xlabel('Time')
    ylabel(['u_',num2str(ii)])
end