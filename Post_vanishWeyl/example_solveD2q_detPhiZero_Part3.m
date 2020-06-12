load('Data_solveD2q_detPhiZero.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppose d2q_pp=/=0;
syms t s %(t=/=0, |t|^2=1)
% assume(t, 'clear');
% assume(s, 'clear');
% 'clear all' to delete assumptions.
%assume(t~=0);
%assumeAlso((s^2-conj(s)^2*t^2)==0);
MVar02 = [p, y, d2q_pp, s, t];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2Q_conjpconjp = conj(t)^2*d2q_pp; 
d2Q_py = 2*s*t*d2Q_conjpconjp;
d2Q_conjpy = 2*s*d2Q_conjpconjp;
d2Q_pconjp = -s^2*d2Q_conjpconjp;
d2Q_yy = (-2*t+2*s^2)*d2Q_conjpconjp;

variableSet = [d2q_py, d2q_pconjp, d2q_yy, conj(d2q_pp), conj(d2q_py)];
subSet = [d2Q_py, d2Q_pconjp, d2Q_yy, d2Q_conjpconjp, d2Q_conjpy];
totalMat02 = sym('totalMat02', [9,6]);
for j=1:9
    for k=1:6
        temp = totalMat(j,k);
        temp = subs(temp, variableSet, subSet);
        temp = myRealVariableFun(temp, realVariable);
        totalMat02(j,k) = complex_simple3(temp, MVar02);
    end
end
clearvars j k temp
step01 = myRowOperation(totalMat02, 2, 1, 1, -2, MVar02);
step02 = myRowOperation(step01, 4, 1, 1, -2*t, MVar02);
step03 = myRowOperation(step02, 5, 1, 1, -2*s, MVar02);
step04 = myRowOperation(step03, 4, 2, 2, -t, MVar02);
step05 = myRowOperation(step04, 5, 2, 2, -s, MVar02);
step06 = myRowOperation(step05, 7, 2, 2, t*conj(t)^2, MVar02);
step07 = myRowOperation(step06, 8, 2, 2, s*t*conj(t)^2, MVar02);
step08 = myRowOperation(step07, 1, 4, 2, 1, MVar02);
% 3rd column
step09 = myRowOperation(step08, 4, 1, 3, -t, MVar02);
step10 = myRowOperation(step09, 5, 1, 3, -s, MVar02);
step11 = myRowOperation(step10, 3, 1, 6, -s, MVar02);
step12 = myRowOperation(step11, 7, 1, 6, -s*t*conj(t)^2, MVar02);
step13 = myRowOperation(step12, 8, 1, 6, -s^2*t*conj(t)^2, MVar02);
step14 = myRowOperation(step13, 9, 1, 6, -t, MVar02);
step15 = myRowOperation(step14, 1, 1, 6, -s, MVar02);
step16 = myRowOperation(step15, 2, 1, 6, -s, MVar02);
% 4th column; use t*conj(t)=1;
step16(7,4) = 0;
step16(8,4) = 0;
step16(9,4) = 0;
step17 = myRowOperation(step16, 8, 1, 7, -s, MVar02);
step18 = myRowOperation(step17, 7, 1, 4, -conj(t)^2, MVar02);
step18(4,5) = -4*d2q_pp*(s^2*conj(t)-1);
step18(5,5) = -4*d2q_pp*s*conj(t)*(s^2*conj(t)-1);
step19 = myRowOperation(step18, 5, 1, 4, -s*conj(t), MVar02);
step19(5,6) = 0;
step19(1,:) = step19(1,:)/conj(t)^2;
step19(2,2) = -4*t^2*d2q_pp; 
step19(2,4:6) = step19(2,4:6)/conj(t)^2;
step19(6,:) = step19(6,:)/conj(t)^2;
step19(4,6) = -2*d2q_pp*s^2*conj(t); 
% s = conj(s)*t
step19(4,5) = -4*d2q_pp*(s*conj(s)- 1);
step19(4,6) = -2*d2q_pp*s*conj(s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%