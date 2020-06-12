% example5_Part2_NatTwoR_CR_rmW.m
variable_example5_Part2_CR_rmW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, d2phi_conjpconjp, d2phi_conjpy;
    d2phi_py, d2phi_conjpy, d2phi_yy];
d2q = [d2q_pp, d2q_pconjp, d2q_py;
    d2q_pconjp, d2q_conjpconjp, d2q_conjpy;
    d2q_py, d2q_conjpy, d2q_yy];
Cmat = [2*d2q_conjpconjp, -2*d2q_pconjp-d2q_yy, d2q_conjpy;
    -2*d2q_pconjp-d2q_yy, 2*d2q_pp, d2q_py;
    d2q_conjpy, d2q_py, -2*d2q_pconjp];
% det_of_Cmat = det(Cmat);
% det_of_Cmat = complex_simple3(det_of_Cmat, MVar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% focus on <<det(d2phi)=0>>
totalMat = d2phi*Cmat; 
phiVec1 = [d2phi_py; d2phi_conjpy; d2phi_yy];
sysMat1 = [d2q_conjpy, d2q_py, -2*d2q_pconjp; %(3,3)
    2*d2q_conjpconjp, -d2q_yy-2*d2q_pconjp, d2q_conjpy; %(3,1)
    -d2q_yy-2*d2q_pconjp, 2*d2q_pp, d2q_py]; %(3,2)
% (A1) d2q_pp=/=0
step01 = myRowOperation(sysMat1, 1, 2*d2q_pp, 3, -d2q_py, MVar); 
step02 = myRowOperation(step01, 2, 2*d2q_pp, 3, d2q_yy+2*d2q_pconjp, MVar);
% (A12)-d2q_py^2-4*d2q_pp*d2q_pconjp =/=0
step03 =  myRowOperation(step02, 2, d2q_py^2+4*d2q_pp*d2q_pconjp,...
    1, step02(2,3), MVar);
step04 = myRowOperation(step03, 3, d2q_py^2+4*d2q_pp*d2q_pconjp,...
    1, step02(3,3), MVar);

% zeroTerm01 = step04(2,1);
% zeroTerm01 = expand(zeroTerm01);
% zeroTerm01 = complex_simple3(zeroTerm01, MVar);
step04(2,1)=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testZ01 = zeroTerm01 +2*d2q_pp*det_of_Cmat02;
% testZ01 = complex_simple3(testZ01, MVar);
% let d2phi_py = t*(d2q_py^2+4*d2q_pp*d2q_pconjp);
d2phi_pySoln = tFun*(d2q_py^2+4*d2q_pp*d2q_pconjp);
d2phi_conjpySoln = -step04(3,1)/step04(3,2)*d2phi_pySoln;
d2phi_yySoln = -step04(1,1)/step04(1,3)*d2phi_pySoln;
d2phi_conjpySoln = complex_simple3(d2phi_conjpySoln, MVar);
d2phi_yySoln = complex_simple3(d2phi_yySoln, MVar);

solnVec = [d2phi_pySoln; d2phi_conjpySoln; d2phi_yySoln];
tempValue = subs(d2phi_yySoln,tFun,1);
% normalize solnVec so that d2phi_conjpy = conj(d2phi_py);
solnVec2 = solnVec*conj(tempValue);
for j=1:3
    solnVec2(j)= myRealVariableFun(solnVec2(j),realVariable);
    solnVec2(j)=expand(solnVec2(j));
    solnVec2(j)=complex_simple3(solnVec2(j),MVar);
end
clearvars tempValue j
% forceZero
% forceZ1 = det_of_Cmat02;
% forceZ2 = subs(solnVec2(2)-conj(solnVec2(1)),t,1);
% forceZ2 = expand(forceZ2);
% forceZ2 = myRealVariableFun(forceZ2,realVariable);
% forceZ2 = complex_simple3(forceZ2,MVar);
% testfZ02 = forceZ2 -d2q_conjpy*(det_of_Cmat02);
% testfZ02 = complex_simple3(testfZ02);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
termA = subs(solnVec2(1),tFun,1);
termB = subs(solnVec2(3),tFun,1);
dtermA = df_example_pconjpy(termA,CVar,derivativeDict);
dtermB = df_example_pconjpy(termB,CVar,derivativeDict);
termC = dtermA(3)-dtermB(1);
termC = complex_simple3(termC,MVar);
equationOne = termA*dtFun_y -termB*dtFun_p +termC*tFun;
save('Data_example5_Part2.mat');
