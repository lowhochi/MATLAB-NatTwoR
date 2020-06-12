% compare Chern(1,1,1,2) and conj(Chern(1,1,2,1)).
% b/f: Chern3_Sep24_CR_rmW.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Data_ChernThree_Section2_Sep24_2019.mat');
assumeAlso([x y z theta], 'real');
Ch1112 = Chern_in_aV(1,1,1,2);
Ch1121 = Chern_in_aV(1,1,2,1);
test = Ch1112 - conj(Ch1121);
test = subs(test, Y0, 1+u*conj(u));
test = complex_simple3(test, MVarCh22);
% Evaluate test by Gij_k and its derivatives.
% WeylCurv_Part5_NatTwo_CR_rmW.m
% variableGijk_NatTwo_CR_rmW.m
temp_Sep24_CR_rmW
% assumeAlso([G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3,...
%     G33_1, G33_2, dG12_3_vnormv, dG23_1_vnormv, dG31_2_vnormv,...
%     dG11_2_vnormv, dG11_3_vnormv, dG22_1_vnormv, dG22_3_vnormv,...
%     dG33_1_vnormv, dG33_2_vnormv],'real');

thetaSub = G12_3 + G23_1 + G31_2;
dtheta_vnormvSub = dG12_3_vnormv + dG23_1_vnormv + dG31_2_vnormv;
aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);

daVVec = df_NatTwo_updated_CR_rmW(aVSub, CVarTest, dCVarTest);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);
daV_vnormvSub = daVVec(3);
d2aV_conjuVec = df_NatTwo_updated_CR_rmW(daV_conjuSub, CVarTest, dCVarTest);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);

variableTest1 = [aV, daV_u, daV_conju, daV_vnormv, d2aV_conjuconjMu,...
    theta, dtheta_vnormv];
subSetTest1 = [aVSub, daV_uSub, daV_conjuSub, daV_vnormvSub, d2aV_conjuconjMuSub,...
    thetaSub, dtheta_vnormvSub];
test = subs(test, variableTest1, subSetTest1);

test01 = test;
save('DataOne_Sep26.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
load('DataOne_Sep26.mat');
assumeAlso([x y z theta], 'real');
assumeAlso(Gijk,'real');
assumeAlso(d2GijkE,'real');

% variableGijk_NatTwo_CR_rmW.m
variableMu = [dG12_3_Mu, dG23_1_Mu, dG31_2_Mu,...
    dG11_2_Mu, dG11_3_Mu, dG22_1_Mu, dG22_3_Mu, dG33_1_Mu, dG33_2_Mu];
variableConjMu = [dG12_3_conjMu, dG23_1_conjMu, dG31_2_conjMu,...
    dG11_2_conjMu, dG11_3_conjMu, dG22_1_conjMu, dG22_3_conjMu,...
    dG33_1_conjMu, dG33_2_conjMu];
variableVnormv = [dG12_3_vnormv, dG23_1_vnormv, dG31_2_vnormv,...
    dG11_2_vnormv, dG11_3_vnormv, dG22_1_vnormv, dG22_3_vnormv,...
    dG33_1_vnormv, dG33_2_vnormv];

subSetMu = [dG12_3_MuSub, dG23_1_MuSub, dG31_2_MuSub,...
    dG11_2_MuSub, dG11_3_MuSub, dG22_1_MuSub, dG22_3_MuSub,...
    dG33_1_MuSub, dG33_2_MuSub];
subSetConjMu = [dG12_3_conjMuSub, dG23_1_conjMuSub, dG31_2_conjMuSub,...
    dG11_2_conjMuSub, dG11_3_conjMuSub, dG22_1_conjMuSub, dG22_3_conjMuSub,...
    dG33_1_conjMuSub, dG33_2_conjMuSub];
subSetVnormv = [dG12_3_vnormvSub, dG23_1_vnormvSub, dG31_2_vnormvSub,...
    dG11_2_vnormvSub, dG11_3_vnormvSub, dG22_1_vnormvSub, dG22_3_vnormvSub,...
    dG33_1_vnormvSub, dG33_2_vnormvSub];

test = subs(test, variableMu, subSetMu);
test = subs(test, variableConjMu, subSetConjMu);
test = subs(test, variableVnormv, subSetVnormv);
test = complex_simple3(test, [u]);
save('DataTwo_Sep26.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% load('DataTemp_Sep26.mat');
% latex_Test = fopen('latex_Test1112_Sep26_2019.txt','w');
% fclose(latex_Test);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%