% WeylFlat5_NatTwoR_CR_rmW.m
% Verify Theorem 7.3
load('DataWeylFlat4_NatTwoR_CR_rmW.mat');
variableSet = [Psi, dPsi_p, d2Psi_pp, K0, dK0_y,...
    Zeta1, dZeta1_p, d2Zeta1_pp, Zeta2, dZeta2_p];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soln (1.11)
% syms Psi0 % constant
% PsiSub = Psi0;
% dPsi_pSub = 0;
% d2Psi_ppSub = 0;
% Zeta1Sub =0;
% dZeta1_pSub=0;
% d2Zeta1_ppSub=0;
% Zeta2Sub = 0;
% dZeta2_pSub = 0;
% subSet = [PsiSub, dPsi_pSub, d2Psi_ppSub, K0, dK0_y,...
%     Zeta1Sub, dZeta1_pSub, d2Zeta1_ppSub, Zeta2Sub, dZeta2_pSub];
% MVarSet = [u, p, y, Psi0, K0, dK0_y];
% for number = 1:4
%     myTest = WflatSet03(number);
%     myTest = subs(myTest, variableSet, subSet);
%     myTest = complex_simple3(myTest, MVarSet);
%     disp(myTest);
% end
% clearvars myTest number MVarSet subSet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soln (1.12)
% syms Psi0% constant
% PsiSub = Psi0;
% dPsi_pSub = 0;
% d2Psi_ppSub = 0;
% Zeta1Sub =0;
% dZeta1_pSub=0;
% d2Zeta1_ppSub=0;
% Zeta2Sub = Zeta2;
% dZeta2_pSub = 0;
% K0Sub = K0;
% dK0_ySub = 0;
% subSet = [PsiSub, dPsi_pSub, d2Psi_ppSub, K0Sub, dK0_ySub,...
%     Zeta1Sub, dZeta1_pSub, d2Zeta1_ppSub, Zeta2Sub, dZeta2_pSub];
% MVarSet = [u, p, y, Psi0, Zeta2, K0];
% for number = 1:4
%     myTest = WflatSet03(number);
%     myTest = subs(myTest, variableSet, subSet);
%     myTest = complex_simple3(myTest, MVarSet);
%     disp(myTest);
% end
% clearvars myTest number MVarSet subSet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soln (1.21)
% syms Const1 Const2 Const3 %constant
% PsiSub = Const1*p +Const2;
% dPsi_pSub = Const1;
% d2Psi_ppSub = 0;
% Zeta1Sub =0;
% dZeta1_pSub=0;
% d2Zeta1_ppSub = 0;
% Zeta2Sub = 0;
% dZeta2_pSub = 0;
% K0Sub = Const3/y;
% dK0_ySub = -Const3/y^2;
% subSet = [PsiSub, dPsi_pSub, d2Psi_ppSub, K0Sub, dK0_ySub,...
%     Zeta1Sub, dZeta1_pSub, d2Zeta1_ppSub, Zeta2Sub, dZeta2_pSub];
% MVarSet = [u, p, y, Const1, Const2, Const3];
% for number = 1:4
%     myTest = WflatSet03(number);
%     myTest = subs(myTest, variableSet, subSet);
%     myTest = subs(myTest, conj(Const3), Const3);
%     myTest = complex_simple3(myTest, MVarSet);
%     disp(myTest);
% end
% clearvars myTest number MVarSet subSet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soln (1.22)
% syms Const1 Const2 Const3 Const4%constant
% PsiSub = Const1*p +Const2;
% dPsi_pSub = Const1;
% d2Psi_ppSub = 0;
% Zeta1Sub =0;
% dZeta1_pSub=0;
% d2Zeta1_ppSub = 0;
% Zeta2Sub = i*Const4*Const1;
% dZeta2_pSub = 0;
% K0Sub = Const3/(y-Const4);
% dK0_ySub = -Const3/(y-Const4)^2;
% subSet = [PsiSub, dPsi_pSub, d2Psi_ppSub, K0Sub, dK0_ySub,...
%     Zeta1Sub, dZeta1_pSub, d2Zeta1_ppSub, Zeta2Sub, dZeta2_pSub];
% MVarSet = [u, p, y, Const1, Const2, Const3, Const4];
% for number = 1:4
%     myTest = WflatSet03(number);
%     myTest = subs(myTest, variableSet, subSet);
%     myTest = subs(myTest, conj(Const3), Const3);
%     myTest = subs(myTest, conj(Const4), Const4);
%     myTest = complex_simple3(myTest, MVarSet);
%     disp(myTest);
% end
% clearvars myTest number MVarSet subSet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soln (2.1)
% syms constA1 constA2 constB %constant
% PsiSub = 0;
% dPsi_pSub = 0;
% d2Psi_ppSub = 0;
% Zeta1Sub = constA1*p +constA2;
% dZeta1_pSub = constA1;
% d2Zeta1_ppSub = 0;
% Zeta2Sub = -constA1^2/conj(constA1)*p...
%     -2*constA1*constA2/(conj(constA1))...
%     +constB*constA1^2;
% dZeta2_pSub = -constA1^2/conj(constA1);
% K0Sub = 0;
% dK0_ySub = 0;
% subSet = [PsiSub, dPsi_pSub, d2Psi_ppSub, K0Sub, dK0_ySub,...
%     Zeta1Sub, dZeta1_pSub, d2Zeta1_ppSub, Zeta2Sub, dZeta2_pSub];
% MVarSet = [u, p, y, constA1, constA2, constB];
% for number = 1:4
%     myTest = WflatSet03(number);
%     myTest = subs(myTest, variableSet, subSet);
%     myTest = subs(myTest, conj(constB), constB);
%     myTest = complex_simple3(myTest, MVarSet);
%     disp(myTest);
% end
% clearvars myTest number MVarSet subSet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soln (2.2)
syms Psi0 constR constB constA1 %constant
constA = constA1...
    +i*(conj(Psi0)^2*constB-Psi0^2*conj(constB))/(4*i*Psi0*conj(Psi0));
PsiSub = Psi0;
dPsi_pSub = 0;
d2Psi_ppSub = 0;
Zeta1Sub = i*constR*Psi0*p +constA;
dZeta1_pSub = i*constR*Psi0;
d2Zeta1_ppSub = 0;
Zeta2Sub = i*constR*Psi0^2/conj(Psi0)*p +constB;
dZeta2_pSub = i*constR*Psi0^2/conj(Psi0);
K0Sub = 0;
dK0_ySub = 0;
subSet = [PsiSub, dPsi_pSub, d2Psi_ppSub, K0Sub, dK0_ySub,...
    Zeta1Sub, dZeta1_pSub, d2Zeta1_ppSub, Zeta2Sub, dZeta2_pSub];
MVarSet = [u, p, y, Psi0 constR, constB, constA1];
for number = 1:4
    myTest = WflatSet03(number);
    myTest = subs(myTest, variableSet, subSet);
    myTest = subs(myTest, conj(constR), constR);
    myTest = subs(myTest, conj(constA1), constA1);
    myTest = complex_simple3(myTest, MVarSet);
    disp(myTest);
end
clearvars myTest number MVarSet subSet















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%