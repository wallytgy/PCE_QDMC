clc, clear all, close all
addpath(genpath('../../../PoCET'));
%%

k_01 = 1.287*1e12
k_02 = 1.287*1e12
k_03 = 9.043*1e9
E_A1_over_R = 9758.3
E_A2_over_R = 9758.3
E_A3_over_R = 7704.0
DeltaH_AB = 4.2
DeltaH_BC = -11.0
DeltaH_AD = -41.85
rho = 0.9342
c_p = 3.01
c_pK = 2
A = 0.215
V_R = 10.01
m_k = 5.0
T_in = 130
k_W = 4032
dot_Q_k = -4500
c_A0 = 2.5 %% not mentioned in paper




%%

states(1).name = 'c_A';
states(1).dist = 'none';
states(1).data = 0.8;
states(1).rhs = 'F*(2.5 - c_A) - (1+0.1*b)*c_A *c_B - (1+0.1*a) *c_A^2';

states(2).name = 'c_B';
states(2).dist = 'none';
states(2).data = 0.5;
states(2).rhs = '-F *c_B +  (1+0.1*b)*c_A*c_B  -  c_B';

states(3).name = 'T_R';
states(3).dist = 'none';
states(3).data = 134.14;
states(3).rhs = 'F *(130 - T_R) + (T_K - T_R) - 100*((1+0.1*b)*c_A*c_B  + c_B + (1+0.1*a) *c_A^2)'

states(4).name = 'T_K';
states(4).dist = 'none';
states(4).data = 134;
states(4).rhs = '(100 + (T_R - T_K))';

%%
parameter_names = {'a','b'}
parameter_dist{1} = 'uniform'
parameter_dist{2} = 'uniform'

parameter_data =  {[-1,1], -[1,1]}

parameters = struct('name', parameter_names, 'dist', parameter_dist, 'data', parameter_data)

%%
inputs(1).name = 'F';
inputs(1).rhs  = 'piecewise(u_t, u_v, t)';
inputs(1).u_t = [-0.1]
inputs(1).u_v = [5]

%% define simulation options
simoptions.tspan = [0, 3];
simoptions.dt = 0.01;
simoptions.setup = odeset;
simoptions.solver = 'ode45';

mc_samples = 1e4;
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
subplot(2,2,4); plot(results.time(2:end),results.c_B.moments(4,2:end),'r',mcresults.time(2:end),mcresults.c_B.moments(4,2:end),'b'); title('Excess Kurtosis')