function dXdt = my_MCODE(t,X,PAR,u_t,u_v)

 c_A = X(1);
 c_B = X(2);
 c_C = X(3);
 c_D = X(4);
 
 u = piecewise(u_t, u_v, t);
 
 ddt_c_A = (2.75+u)/5*0.3 -(2.75+u+1)/5*c_A -c_A*c_B*(PAR.k_1+PAR.k_2);
 ddt_c_B = 1/5*7.35 - (2.75+u+1)/5*c_B - c_A*c_B*(PAR.k_1+PAR.k_2) - c_B*c_C*0.3264 - c_B*c_D*0.01591;
 ddt_c_C = -(2.75+u+1)/5*c_C + c_A*c_B*PAR.k_1 - c_B*c_C*0.3264;
 ddt_c_D = -(2.75+u+1)/5*c_D + c_A*c_B*PAR.k_2 - c_B*c_D*0.01591;

 dXdt = [ddt_c_A; ddt_c_B; ddt_c_C; ddt_c_D];
end