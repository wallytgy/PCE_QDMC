function dXdt = my_ODE(t,X,PoCETsys,u_t,u_v)
 M = PoCETsys.coeff_matrices;

 c_A = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 c_B = X(1*PoCETsys.pce.options.n_phi+1:2*PoCETsys.pce.options.n_phi);
 c_C = X(2*PoCETsys.pce.options.n_phi+1:3*PoCETsys.pce.options.n_phi);
 c_D = X(3*PoCETsys.pce.options.n_phi+1:4*PoCETsys.pce.options.n_phi);
 
 u = piecewise(u_t, u_v, t);
 
 c_Ac_B = ckron(c_A,c_B);
 c_Bc_C = ckron(c_B,c_C);
 c_Bc_D = ckron(c_B,c_D);
 
 ddt_c_A = - M.k_1_O2*c_Ac_B - M.k_2_O2*c_Ac_B -0.2*u*M.one_O1*c_A -0.75*M.one_O1*c_A + 0.06*u*M.one_O0 + 0.165*M.one_O0;
 ddt_c_B = - M.k_1_O2*c_Ac_B - M.k_2_O2*c_Ac_B -3.264*M.one_O2*c_Bc_C -0.01591*M.one_O2*c_Bc_D -0.2*u*M.one_O1*c_B -0.75*M.one_O1*c_B + 1.47*M.one_O0;
 ddt_c_C = M.k_1_O2*c_Ac_B -3.264*M.one_O2*c_Bc_C -0.2*u*M.one_O1*c_C -0.75*M.one_O1*c_C;
 ddt_c_D = M.k_2_O2*c_Ac_B -0.01591*M.one_O2*c_Bc_D -0.2*u*M.one_O1*c_D -0.75*M.one_O1*c_D;
 
 dXdt = [ddt_c_A; ddt_c_B; ddt_c_C; ddt_c_D];
end