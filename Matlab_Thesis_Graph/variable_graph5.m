% variable_graph5.m
syms x y z
Gamma = sym('Gamma',[3 3 3]);
Gamma(1,1,1) = 1/(x+y*z);
Gamma(1,1,2) = -z*(x+y*z);
Gamma(1,1,3) = -y*(x+y*z);
Gamma(1,2,1) = z/(x+y*z);
Gamma(1,2,2) = 0;
Gamma(1,2,3) = 0;
Gamma(1,3,1) = y/(x+y*z);
Gamma(1,3,2) = 0;
Gamma(1,3,3) = 0;

Gamma(2,1,1) = z/(x+y*z);
Gamma(2,1,2) = 0;
Gamma(2,1,3) = 0;
Gamma(2,2,1) = 0;
Gamma(2,2,2) = 0;
Gamma(2,2,3) = 0;
Gamma(2,3,1) = 0;
Gamma(2,3,2) = 0;
Gamma(2,3,3) = 0;

Gamma(3,1,1) = y/(x+y*z);
Gamma(3,1,2) = 0;
Gamma(3,1,3) = 0;
Gamma(3,2,1) = 0;
Gamma(3,2,2) = 0;
Gamma(3,2,3) = 0;
Gamma(3,3,1) = 0;
Gamma(3,3,2) = 0;
Gamma(3,3,3) = 0;