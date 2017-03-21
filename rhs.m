function yp = rhs (time,y)
% Computes the right-hand side of equation given t and P(t)
% a = -1;
% B = 2;
% g_0 = -0.9;
% g_1= -0.0001;
% 
% yp = y*(a+(B*y)+(g_0+g_1*time)*y^2);
yp = -y + 2*y^2;

end



