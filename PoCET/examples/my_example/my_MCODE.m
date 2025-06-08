function dXdt = my_MCODE(t,X,PAR,u_t,u_v)

 c_A = X(1);
 c_B = X(2);
 T_R = X(3);
 T_K = X(4);
 
 F = piecewise(u_t, u_v, t);
 
 ddt_c_A = F*(2.5 - c_A) - (1+0.1*PAR.b)*c_A *c_B - (1+0.1*PAR.a) *c_A^2;
 ddt_c_B = -F *c_B +  (1+0.1*PAR.b)*c_A*c_B  -  c_B;
 ddt_T_R = F *(130 - T_R) + (T_K - T_R) - 100*((1+0.1*PAR.b)*c_A*c_B  + c_B + (1+0.1*PAR.a) *c_A^2);
 ddt_T_K = (100 + (T_R - T_K));

 dXdt = [ddt_c_A; ddt_c_B; ddt_T_R; ddt_T_K];
end