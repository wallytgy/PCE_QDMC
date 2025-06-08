function fout = ex2_OUT(t,X,PoCETsys,ut1,uv1,uv2)
 M = PoCETsys.coeff_matrices;

 x1 = X(0*PoCETsys.pce.options.n_phi+1:1*PoCETsys.pce.options.n_phi);
 x2 = X(1*PoCETsys.pce.options.n_phi+1:2*PoCETsys.pce.options.n_phi);
 x3 = X(2*PoCETsys.pce.options.n_phi+1:3*PoCETsys.pce.options.n_phi);
 
 u1 = piecewise(ut1,uv1,t);
 u2 = -uv2;
 
 
 y = M.one_O1*x3;
 
 fout = [y];
end