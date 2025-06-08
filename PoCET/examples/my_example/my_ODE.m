function dXdt = my_ODE(t,X,PoCETsys,u_t,u_v)
 M = PoCETsys.coeff_matrices;

 c_A = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 c_B = X(1*PoCETsys.pce.options.n_phi+1:2*PoCETsys.pce.options.n_phi);
 T_R = X(2*PoCETsys.pce.options.n_phi+1:3*PoCETsys.pce.options.n_phi);
 T_K = X(3*PoCETsys.pce.options.n_phi+1:4*PoCETsys.pce.options.n_phi);
 
 F = piecewise(u_t, u_v, t);
 
 c_Ac_A = ckron(c_A,c_A);
 c_Ac_B = ckron(c_A,c_B);
 
 ddt_c_A = - F*M.one_O1*c_A + 2.5*F*M.one_O0 -0.1*M.a_O2*c_Ac_A -0.1*M.b_O2*c_Ac_B - M.one_O2*c_Ac_A - M.one_O2*c_Ac_B;
 ddt_c_B = - F*M.one_O1*c_B + 0.1*M.b_O2*c_Ac_B + M.one_O2*c_Ac_B - M.one_O1*c_B;
 ddt_T_R = - F*M.one_O1*T_R + 130*F*M.one_O0 + M.one_O1*T_K - M.one_O1*T_R -10*M.a_O2*c_Ac_A -10*M.b_O2*c_Ac_B -100*M.one_O2*c_Ac_A -100*M.one_O2*c_Ac_B -100*M.one_O1*c_B;
 ddt_T_K = - M.one_O1*T_K + M.one_O1*T_R + 100*M.one_O0;
 
 dXdt = [ddt_c_A; ddt_c_B; ddt_T_R; ddt_T_K];
end