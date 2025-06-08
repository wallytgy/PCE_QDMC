function dXdt = ex5_ODE_org(t,X,PAR,ut1,uv1,ut2,uv2)

 x1 = X(1);
 x2 = X(2);
 x3 = X(3);
 
 u1 = piecewise(ut1,uv1,t);
 u2 = piecewise(ut2,uv2,t);
 
 ddt_x1 = (-PAR.c1*PAR.Sp*sign(x1-x3)*sqrt(2*PAR.g*abs(x1-x3))+u1)/PAR.A;
 ddt_x2 = (-PAR.c3*PAR.Sp*sign(x2-x3)*sqrt(2*PAR.g*abs(x2-x3))-PAR.c2*PAR.Sp*sqrt(2*PAR.g*x2)-PAR.Qf+u2)/PAR.A;
 ddt_x3 = (PAR.c1*PAR.Sp*sign(x1-x3)*sqrt(2*PAR.g*abs(x1-x3))-PAR.c3*PAR.Sp*sign(x3-x2)*sqrt(2*PAR.g*abs(x3-x2)))/PAR.A;

 dXdt = [ddt_x1; ddt_x2; ddt_x3];
end