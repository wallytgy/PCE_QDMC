function dXdt = my_MCODE(t,X,PAR,u_t,u_v)

 c_A = X(1);
 c_B = X(2);
 T_R = X(3);
 T_K = X(4);
 
 F = piecewise(u_t, u_v, t);
 
 ddt_c_A = F  (PAR.c_A0 - c_A) - (1+0.1*PAR.b)  PAR.k_01  exp(-PAR.E_A1/( PAR.R *(T_R + 273.15) ))  c_A - PAR.k_03  exp(-(1+0.1  PAR.a)  PAR.E_A3/( PAR.R *(T_R + 273.15) ) )   c_A^2;
 ddt_c_B = -F  c_B +  (1+0.1*PAR.b)  PAR.k_01  exp(-PAR.E_A1/( PAR.R *(T_R + 273.15) ))  c_A - PAR.k_02   exp(-PAR.E_A2/( PAR.R *(T_R + 273.15) ))   c_B;
 ddt_T_R = F  (PAR.T_in - T_R) + PAR.k_W  PAR.A / ( PAR.rho  PAR.c_p  PAR.V_R )  (T_K - T_R) - 1 / ( PAR.rho  PAR.c_p )  ((1+0.1*PAR.b)  PAR.k_01  exp(-PAR.E_A1/( PAR.R *(T_R + 273.15) ))  c_A  PAR.DeltaH_AB + PAR.k_02   exp(-PAR.E_A2/( PAR.R *(T_R + 273.15) ))  c_B  PAR.DeltaH_BC + PAR.k_03  exp(-(1+0.1  PAR.a)  PAR.E_A3/( PAR.R *(T_R + 273.15) ) )  c_A^2  PAR.DeltaH_AD);
 ddt_T_K = 1/(PAR.m_K  PAR.c_pK)  (PAR.Q_K + PAR.k_W  PAR.A  (T_R - T_K));

 dXdt = [ddt_c_A; ddt_c_B; ddt_T_R; ddt_T_K];
end