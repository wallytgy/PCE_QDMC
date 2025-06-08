function dXdt = my_ODE(t,X,PoCETsys,u_t,u_v)
 M = PoCETsys.coeff_matrices;

 c_A = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 c_B = X(1*PoCETsys.pce.options.n_phi+1:2*PoCETsys.pce.options.n_phi);
 T_R = X(2*PoCETsys.pce.options.n_phi+1:3*PoCETsys.pce.options.n_phi);
 T_K = X(3*PoCETsys.pce.options.n_phi+1:4*PoCETsys.pce.options.n_phi);
 
 F = piecewise(u_t, u_v, t);
 
 
 ddt_c_A = *M.one_O0;
 ddt_c_B = *M.one_O0;
 ddt_T_R = *M.one_O0;
 ddt_T_K = *M.one_O0;
 
 dXdt = [ddt_c_A; ddt_c_B; ddt_T_R; ddt_T_K];
end