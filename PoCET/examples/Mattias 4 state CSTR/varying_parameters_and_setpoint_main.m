clc, clear all, close all
addpath(genpath('../../../PoCET'));
%%

q_1ss = 1.1
q_2 = 1
V = 5
k_3 = 3.264
k_4 = 0.01591
c_Ain = 1
c_Bin = 7.35





%%

states(1).name = 'c_A';
states(1).dist = 'none';
states(1).data = 0;
states(1).rhs = '(2.75+u)/5*0.3 -(2.75+u+1)/5*c_A -c_A*c_B*(k_1+k_2)';

states(2).name = 'c_B';
states(2).dist = 'none';
states(2).data = 3.5;
states(2).rhs = '1/5*7.35 - (2.75+u+1)/5*c_B - c_A*c_B*(k_1+k_2) - c_B*c_C*3.264 - c_B*c_D*0.01591';

states(3).name = 'c_C';
states(3).dist = 'none';
states(3).data = 0;
states(3).rhs = '-(2.75+u+1)/5*c_C + c_A*c_B*k_1 - c_B*c_C*3.264'

states(4).name = 'c_D';
states(4).dist = 'none';
states(4).data = 0.0025;
states(4).rhs = '-(2.75+u+1)/5*c_D + c_A*c_B*k_2 - c_B*c_D*0.01591';

%%
parameter_names = {'k_1','k_2'}
parameter_dist{1} = 'uniform'
parameter_dist{2} = 'uniform'

parameter_data =  {[0.2789,0.8927], [0.1894,0.9331]}

parameters = struct('name', parameter_names, 'dist', parameter_dist, 'data', parameter_data)

%%
inputs(1).name = 'u';
inputs(1).rhs  = 'piecewise(u_t, u_v, t)';
inputs(1).u_t = [-1]
inputs(1).u_v = [0]

%% define simulation options
simoptions.tspan = [0, 40];
simoptions.dt = 0.1;
simoptions.setup = odeset;
simoptions.solver = 'ode45';

mc_samples = 1e3;
pce_order = 4;


%% (1.c) compose the PCE system and write file for PCE-ODE and MC simulations
% mySystem = PoCETcompose(states, parameters, inputs, options);
sys = PoCETcompose(states,parameters,inputs,[],pce_order);
MomMats = PoCETmomentMatrices(sys,pce_order);
PoCETwriteFiles(sys,'my_ODE.m','my_OUT.m','my_MCODE.m','my_MCOUT.m')

%% (2.a) simulate system and compute moments
% run PCE simulation
results = PoCETsimGalerkin(sys,'my_ODE',[],simoptions);

% compute moments from simulation results and store also in results
results = PoCETcalcMoments(sys,MomMats,results);

% run Monte-Carlo simulations
samples = PoCETsample(sys,'variables',mc_samples);
mcresults = PoCETsimMonteCarlo(sys,'my_MCODE',[],samples,simoptions,'method','moments');


%%
% plot results of both simulations
figure(1) % central moments of state 1
subplot(2,2,1); plot(results.time,results.c_A.moments(1,:),'r',mcresults.time,mcresults.c_A.moments(1,:),'b'); title('Mean')
subplot(2,2,2); plot(results.time,results.c_A.moments(2,:),'r',mcresults.time,mcresults.c_A.moments(2,:),'b'); title('Variance')
subplot(2,2,3); plot(results.time(2:end),results.c_A.moments(3,2:end),'r',mcresults.time(2:end),mcresults.c_A.moments(3,2:end),'b'); title('Skewness')
subplot(2,2,4); plot(results.time(2:end),results.c_A.moments(4,2:end),'r',mcresults.time(2:end),mcresults.c_A.moments(4,2:end),'b'); title('Excess Kurtosis')

%%
figure(2) % central moments of state 2
subplot(2,2,1); plot(results.time,results.c_B.moments(1,:),'r',mcresults.time,mcresults.c_B.moments(1,:),'b'); title('Mean')
subplot(2,2,2); plot(results.time,results.c_B.moments(2,:),'r',mcresults.time,mcresults.c_B.moments(2,:),'b'); title('Variance')
subplot(2,2,3); plot(results.time(2:end),results.c_B.moments(3,2:end),'r',mcresults.time(2:end),mcresults.c_B.moments(3,2:end),'b'); title('Skewness')


%%

%%%% first construct the steady states by using asymptotic value of the
%%%% above simulation


c_A_pce_ss = results.c_A.pcvals(:,end);
c_B_pce_ss = results.c_B.pcvals(:,end);
c_C_pce_ss = results.c_C.pcvals(:,end);
c_D_pce_ss = results.c_D.pcvals(:,end);

X_ss = [c_A_pce_ss;c_B_pce_ss;c_C_pce_ss;c_D_pce_ss];


X_ss = fsolve(@(X) my_ODE(0,X,sys,inputs(1).u_t,inputs(1).u_v) , X_ss)
display('checking gradient is 0:')
norm(my_ODE(0,X_ss,sys,inputs(1).u_t,inputs(1).u_v))


u = 1
[t, X_response] = ode45( @(t,X) my_ODE(t,X,sys,-1,u) , 0:0.1:40 ,  X_ss );

figure(3)
for u = -2.5:0.5:2.5
    [t, X_response] = ode45( @(t,X) my_ODE(t,X,sys,-1,u) , 0:0.1:40 ,  X_ss )
    plot(t,X_response(:,31))
    hold on 
end
xlabel('Time/s')
ylabel('c_C/M')
title('Output response for various deviation in flow rate')
legend('-2.5','-2','-1.5','-1','-0.5','0','0.5','1','1.5','2','2.5')
%%%% lets control C


%% generate step response matrix
u = 1
[t, X_response] = ode45( @(t,X) my_ODE(t,X,sys,-1,u) , 0:0.1:40 ,  X_ss );
X_response = X_response';
[n_r , ~] = size(X_response);
% y = c_C, u = F
c_C_PCE = X_response((n_r/4) * 2 + 1: (n_r/4) * 3, :);
S_yu = c_C_PCE(:,6:5:101);
Su = reshape(S_yu, [15,1,20]);
for i = 1:20
    Su(:,1,i) = Su(:,1,i) - X_ss(31:45);
end

%% 
% Parameters to create step response
dt = 0.5;       % Delta t (has to be identical when creating the step 
t_simulation = 400 ;% response model and when running closed loop control)
tsp_sr = 0:dt:10;   % Time span to create the step response
n_sr = length(tsp_sr) - 1;  % Length of the step response model
ops = odeset;   % ODE options


% Parameters to design QDMC (more details, including required dimensions
% can be found in the file QDMC_controller_soft.m under initialization)
n_in = 1;       % Number of inputs
n_out = 15;      % Number of outputs
n_dis = 0;      % Number of measured disturbances

Sd = [];        % Step resposne of measured disturbances
p = 20;        % Prediction horizon (number of time steps)
c = 10;         % Control horizon (number of time steps)


La = 1;         % Weight for input movement
Q = (1./X_ss(31:45).^2);  % Weight for output error
for i = 2:15
    Q(i) = Q(i)* 0.05;
end

ctg = [];       % Cost to go

u_past = 0;     % Past input (assume process starts at steady state)
y_past = zeros(15,1);     % Past output (assume process starts at steady state)
d_past = [];    % Past measured disturbance
u_int = [];     % Integrating inputs
d_int = [];     % Integrating measured disturbances



u_min = -2.5;     % Minimum input
u_max = 2.5;    % Maximum input
D_u_min = -0.1; % Minimum input movement
D_u_max = 0.1;  % Maximum input movement
y_min = -ones(15,1);   % Minimum output
y_max = ones(15,1);    % Maximum output

soft = [];   % Both, y_max and y_min are soft constraints
w_eps = [];     % Weight for linear constraint violation
W_eps = [];     % Weight for quadratic constraint violation

% Parameters to run closed loop simulation
ysp = [Su(1,1,20) * ones(1,t_simulation/dt) ; zeros(14,t_simulation/dt)]; % Output setpoints
ysp(:,201:400) = -0.5 * ysp(:,201:400)
ysp(:,601:800) = -1 * ysp(:,601:800)
udist = []; 

tsp_sim = 0:dt:t_simulation;

%% Perform closed loop control for varying paramters========================================== %


%%% find steady state of current system
for iteration = 1:100
PAR.k_1 = samples.parameters(1,iteration);
PAR.k_2 = samples.parameters(2,iteration);
x_ss_trial =  [X_ss(1), X_ss(16), X_ss(31), X_ss(46)];
x_ss = fsolve(@(x) my_MCODE(0,x,PAR,inputs(1).u_t,inputs(1).u_v) , x_ss_trial);
x_ss = x_ss';
% In this section, a setpoint change and an input disturbance rejection is
% simulated
n_x = 4;

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
x_dev = zeros(n_x,1);   % Start the process from steady state (deviation variable)
y = zeros(n_out,1); % Start the process from steady state (deviation variable)
% Write value for t=0 to storage
x_store(:,1) = x_dev;
y_store(:,1) = y;

% Run closed loop
for ii = 1:length(tsp_sim)-1
    % Current setpoint
    ysp_curr = ysp(:,ii);
    % Current measured disturbance
    d = [];
    
    % Calculate controller output = process input
    u = Controller.run_controller(y,ysp_curr,d);
    
    % Simulate one time step
    x_start = x_ss + x_dev;
    [~,X] = ode45(@(t,x) my_MCODE(0,x,PAR,-1,u) , [0 dt] , x_start ,ops);
    x_dev = X(end,:)' - x_ss ; % State at the end of the time step
    c_C_dev = x_dev(3);
    y = [c_C_dev ; zeros(14,1)];  % Process output
    
    % Store results
    x_store(:,ii+1) = x_dev;
    y_store(:,ii+1) = y;
    u_store(:,ii) = u;
end

%% plot output response
figure(4)
plot(tsp_sim, y_store(1,:), 'r' , tsp_sim(2:end), ysp(1,:), 'b')
xlabel('Time/s')
ylabel('c_C/M')

hold on
end 







