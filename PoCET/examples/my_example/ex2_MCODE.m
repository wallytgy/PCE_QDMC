function dXdt = ex2_MCODE(t,X,PAR,ut1,uv1,uv2)

 x1 = X(1);
 x2 = X(2);
 x3 = X(3);
 
 u1 = piecewise(ut1,uv1,t);
 u2 = -uv2;
 
 ddt_x1 = u1 - PAR.R1*x1;
 ddt_x2 = u2 + PAR.R1*x1 - PAR.R2*x2;
 ddt_x3 = PAR.R2*x2 - PAR.R3*x3;

 dXdt = [ddt_x1; ddt_x2; ddt_x3];
end