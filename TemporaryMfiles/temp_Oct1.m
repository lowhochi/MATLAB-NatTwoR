load('DataRicCurv2_NatTwoR_CR_rmW.mat');
% diffMatrix = sym('difference',[6 6]);
% for m=1:6
%     for n=1:6
%         if (m<n)
%             temp = RicCurv_in_aV(m,n) - RicCurv_in_aV(n,m);
%             temp = subs(temp, Y, u*conj(u)+1);
%             temp = complex_simple3(temp, MVarRicCurv2);   
%             diffMatrix(m,n) = temp;
%         else
%             diffMatrix(m,n) = 0;
%         end
%     end
% end

d12 = RicCurv_in_aV(1,2) - RicCurv_in_aV(2,1);
d15 = RicCurv_in_aV(1,5) - RicCurv_in_aV(5,1);
d25 = RicCurv_in_aV(2,5) - RicCurv_in_aV(5,2);
differenceSet = [d12,d15,d25];
% [aV, d2aV_conjuconjMu, d2aV_conjuvnormv, d2aV_uMu, d2aV_uu,...
%     d2aV_uvnormv, d3aV_uuMu, d4aV_uuconjuconju,...
%     daV_Mu, daV_conjMu, daV_conju, daV_u, daV_vnormv,...
%     dtheta_Mu, dtheta_vnormv, theta];

variable_NatTwoR_Chern2_CR_rmW

variableSet0 = [aM, bV, h11, T4];
subSet0 = [aMSubTwo, bVSubTwo, T4SubTwo, h11SubTwo];

variableTemp1 = [aV, theta];
variableTemp2 = [dtheta_Mu, dtheta_vnormv,...
    daV_Mu, daV_conjMu, daV_conju, daV_u, daV_vnormv, d2aV_uu];
variableTemp3 = [d2aV_conjuconjMu, d2aV_conjuvnormv, d2aV_uMu,...
     d2aV_uvnormv, d3aV_uuMu, d4aV_uuconjuconju];
 
thetaSub = G12_3 + G23_1 + G31_2;
dtheta_vnormvSub = dG12_3_vnormv + dG23_1_vnormv + dG31_2_vnormv;
dtheta_MuSub = dG12_3_Mu + dG23_1_Mu + dG31_2_Mu;

aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);

daVVec = df_NatTwo_MuSet_CR_rmW(aVSub, CVarG1, dCVarG1);
daV_MuSub = daVVec(1);
daV_conjMuSub = daVVec(2);
daV_vnormvSub = daVVec(3);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);

d2aV_uVec = df_NatTwo_MuSet_CR_rmW(daV_uSub, CVarG1, dCVarG1);
d2aV_uMuSub = d2aV_uVec(1);
d2aV_uvnormvSub = d2aV_uVec(3);
d2aV_uuSub = d2aV_uVec(4);

d2aV_conjuVec = df_NatTwo_MuSet_CR_rmW(daV_conjuSub, CVarG1, dCVarG1);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);
d2aV_conjuvnormvSub = d2aV_conjuVec(3);

d3aV_uuVec = df_NatTwo_MuSet_CR_rmW(d2aV_uuSub, CVarG1, dCVarG1);
d3aV_uuMuSub = d3aV_uuVec(1);
d3aV_uuconjuSub = d3aV_uuVec(5);
d4aV_uuconjuconjuSub = complexdiff3(d3aV_uuconjuSub,u,1);

% check that d12, d15 and d25 are zeros.
% variableTemp1 = [aV, theta];
% variableTemp2 = [dtheta_Mu, dtheta_vnormv,...
%     daV_Mu, daV_conjMu, daV_conju, daV_u, daV_vnormv, d2aV_uu];
% variableTemp3 = [d2aV_conjuconjMu, d2aV_conjuvnormv, d2aV_uMu,...
%      d2aV_uvnormv, d3aV_uuMu, d4aV_uuconjuconju];


subSetTemp1 = [aVSub, thetaSub];
subSetTemp2 = [dtheta_MuSub, dtheta_vnormvSub,...
    daV_MuSub, daV_conjMuSub, daV_conjuSub, daV_uSub, daV_vnormvSub,...
    d2aV_uuSub];
subSetTemp3 = [d2aV_conjuconjMuSub, d2aV_conjuvnormvSub, d2aV_uMuSub,...
    d2aV_uvnormvSub, d3aV_uuMuSub, d4aV_uuconjuconjuSub];

MVarTemp = [d2w_uMu, d2w_uconjMu, d2w_uu, d2w_uvnormv,... 
    d3w_uuMu, d3w_uuconjMu, d3w_uuu, dw_Mu, dw_conjMu,...
    dw_u, dw_vnormv, u, w];

for j=1:3
    temp = differenceSet(j);
    temp = subs(temp, Y, 1+u*conj(u));
    temp = subs(temp, variableTemp3, subSetTemp3);
    temp = subs(temp, variableSet0, subSet0);
    temp = subs(temp, variableTemp2, subSetTemp2);
    temp = subs(temp, variableTemp1, subSetTemp1);
    
    temp = subs(temp, variableSetMu, subSetMu);
    temp = subs(temp, variableSetConjMu, subSetConjMu);
    temp = subs(temp, variableSetVnormv, subSetVnormv);
    temp = subs(temp, variableSetBianchi1, subSetBianchi1);
    temp = complex_simple3(temp, MVarTemp);
    differenceSet(j) = temp;
end
