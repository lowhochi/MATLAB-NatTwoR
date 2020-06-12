% fg_cubic_solver.m
syms lapPhi lapQ K real
syms t
f = -t^3 + (lapPhi)/sqrt(2)*t^2 -2*lapQ*t - sqrt(2)*K;
% solve a*t^3 + b*t^2 + c*t + d = 0
a = -1;
b = (lapPhi)/sqrt(2);
c = -2*lapQ;
d = -sqrt(2)*K;
% p q r
p = -b/(3*a);
q = p^3 + (b*c-3*a*d)/(6*a^2);
r = c/(3*a);
p = complex_simple3(p, t);
q = complex_simple3(q, t);
r = complex_simple3(r, t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find first solution by cubic formula
temp = q^2+(r-p^2)^3;
temp = complex_simple3(temp, t);
t1 = p +(q+temp^(1/2))^(1/3) + (q-temp^(1/2))^(1/3);
t1 = complex_simple3(t1, t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find other solutions
f2 = t^2 + (t1-lapPhi/sqrt(2))*t -sqrt(2)*K/t1;
a2 = 1;
b2 = t1-(lapPhi/sqrt(2));
c2 = -sqrt(2)*K/t1;

dtm2 = b2^2-4*a2*c2;
dtm2 = simplify(dtm2);
% t2 = (-b2+(b2^2-4*a2*c2)^(1/2))/(2*a2);
% t3 = (-b2-(b2^2-4*a2*c2)^(1/2))/(2*a2);
% t2 = complex_simple3(t2, t);
% t3 = complex_simple3(t3, t);
% 
% solnSet = [t1, t2, t3];






