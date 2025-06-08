function fout = ex2_MCOUT(t,X_in,PAR,ut1,uv1,uv2)

 x1 = X_in(1,:);
 x2 = X_in(2,:);
 x3 = X_in(3,:);
 
 u1 = piecewise(ut1,uv1,t);
 u2 = -uv2;
 
 y = x3;
 
 fout = [y];
end