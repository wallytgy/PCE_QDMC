function dXdt = ex2_ODE(t,X,PoCETsys,ut1,uv1,uv2)
 M = PoCETsys.coeff_matrices;

 x1 = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 x2 = X(1*PoCETsys.pce.options.n_phi+1:2*PoCETsys.pce.options.n_phi);
 x3 = X(2*PoCETsys.pce.options.n_phi+1:3*PoCETsys.pce.options.n_phi);
 
 u1 = piecewise(ut1,uv1,t);
 u2 = -uv2;
 
 
 ddt_x1 = - M.R1_O1*x1 + u1*M.one_O0;
 ddt_x2 = M.R1_O1*x1 - M.R2_O1*x2 + u2*M.one_O0;
 ddt_x3 = M.R2_O1*x2 - M.R3_O1*x3;
 
 dXdt = [ddt_x1; ddt_x2; ddt_x3];
end