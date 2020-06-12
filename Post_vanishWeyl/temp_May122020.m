syms p y
syms c0 c1 c2
realVariable = [y, c0, c2];
MVar = [p, y, c0, c1, c2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_bottom = (y+c0)^2+(p-c1)*(conj(p)-conj(c1));
lambda = c2*(y+c0)^2/lambda_bottom^2;
a = (-p+c1)/(y+c0);
dlambda_y = complexdiff3(lambda, y, 0);
dlambda_p = complexdiff3(lambda, p ,0);
dlambda_conjp = complexdiff3(lambda, p ,1);
da_p = complexdiff3(a, p, 0);
da_y = complexdiff3(a, y, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testA = 0.5*dlambda_y-conj(da_p)*lambda-conj(a)*dlambda_conjp;
% testA = myRealVariableFun(testA, realVariable);
% testA = complex_simple3(testA, MVar);
% 
% testB = -lambda*conj(a)*conj(da_y)-conj(a)*dlambda_p...
%     -1/2*conj(a)^2*dlambda_y;
% testB = myRealVariableFun(testB, realVariable);
% testB = complex_simple3(testB, MVar);
% 
% testC = dlambda_p + 2*lambda*conj(a)*conj(da_p) +conj(a)^2*dlambda_conjp;
% testC = myRealVariableFun(testC, realVariable);
% testC = complex_simple3(testC, MVar);
% 
% testD = conj(a)*dlambda_y+dlambda_p*(1-a*conj(a))...
%     -lambda*(da_p*conj(a)-conj(da_y));
% testD = myRealVariableFun(testD, realVariable);
% testD = complex_simple3(testD, MVar);
% 
% testD_left = conj(a)*dlambda_y+dlambda_p*(1-a*conj(a));
% testD_left = myRealVariableFun(testD_left, realVariable);
% testD_left = complex_simple3(testD_left, MVar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test on q and its derivatives
q = c2/2*log((y+c0)^2+(p-c1)*(conj(p)-conj(c1)));
dq_p = complexdiff3(q, p, 0);
dq_conjp = conj(dq_p);
dq_y = complexdiff3(q, y, 0);
d2q_pp = complexdiff3(dq_p, p, 0);
d2q_pconjp = complexdiff3(dq_p, p, 1);
d2q_py = complexdiff3(dq_p, y, 0);
d2q_yy = complexdiff3(dq_y, y, 0);

test201 = d2q_py-conj(a)*lambda;
test201 = myRealVariableFun(test201, realVariable);
test201 = complex_simple3(test201, MVar);
%
test202 = d2q_pconjp - lambda/2;
test202 = myRealVariableFun(test202, realVariable);
test202 = complex_simple3(test202, MVar);
%
test203 = 4*d2q_pp*d2q_pconjp + d2q_py^2;
test203 = myRealVariableFun(test203, realVariable);
test203 = complex_simple3(test203, MVar);
%
test204 = 4*d2q_pconjp^2 +2*d2q_pconjp*d2q_yy -d2q_py*conj(d2q_py);
test204 = myRealVariableFun(test204, realVariable);
test204 = complex_simple3(test204, MVar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







