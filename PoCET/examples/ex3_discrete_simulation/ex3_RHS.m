function dXdt = ex3_RHS(t,X,PoCETsys,vars)
 M = PoCETsys.coeff_matrices;

 A = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 B = X(1*PoCETsys.pce.options.n_phi+1:2*PoCETsys.pce.options.n_phi);
 
 
 
 ddt_A = M.p_1_O1*A + 0.5*M.one_O1*B;
 ddt_B = M.p_2_O1*B + M.R_O0;
 
 dXdt = [ddt_A; ddt_B];
end