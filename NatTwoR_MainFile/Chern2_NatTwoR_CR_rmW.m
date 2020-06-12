%load('DatChern1_NatTwoR_CR_rmW.mat');
load('DataChern1_Part2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [aV, theta, dtheta_Mu, dtheta_vnormv, d2theta_MuconjMu,...
%   daV_Mu, daV_conjMu, daV_vnormv, daV_u, daV_conju,...
%   d2aV_conjuMu, d2aV_conjuconjMu, d2aV_uMu, d2aV_uconjMu, d2aV_uu,... 
%   d2aV_MuconjMu]; 

% wSet = [w, dw_conjMu, dw_vnormv, dw_u, d2w_uu, d2w_uconjMu, d2w_conjMuconjMu]
% other variables in Chern_in_aV: Y0, u;
variable_NatTwoR_Chern2_CR_rmW

variableSet0 = [aM, bV, h11, T4];
subSet0 = [aMSubTwo, bVSubTwo, T4SubTwo, h11SubTwo];
% substitute variable0 BEFORE [theta, aV, daV_u, daV_conju, d2aV_uu].

ChernTwo = sym('ChernTwo',[2 2 2 2]);

thetaSub = G12_3 + G23_1 + G31_2;
dtheta_vnormvSub = dG12_3_vnormv + dG23_1_vnormv + dG31_2_vnormv;
dtheta_MuSub = dG12_3_Mu + dG23_1_Mu + dG31_2_Mu;
d2theta_MuconjMuSub = d2G12_3_MuconjMu + d2G23_1_MuconjMu + d2G31_2_MuconjMu;

aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);

daVVec = df_NatTwo_MuSet_CR_rmW(aVSub, CVarG1, dCVarG1);
daV_MuSub = daVVec(1);
% symvar(daV_MuSub)=[dG11_2_Mu, dG11_3_Mu, dG12_3_Mu, dG22_1_Mu,...
%   dG22_3_Mu, dG23_1_Mu, dG31_2_Mu, dG33_1_Mu, dG33_2_Mu, u];
daV_conjMuSub = daVVec(2);
daV_vnormvSub = daVVec(3);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);

d2aV_uVec = df_NatTwo_MuSet_CR_rmW(daV_uSub, CVarG1, dCVarG1);
d2aV_uMuSub = d2aV_uVec(1);
d2aV_uconjMuSub = d2aV_uVec(2);
d2aV_uuSub = d2aV_uVec(4);

d2aV_conjuVec = df_NatTwo_MuSet_CR_rmW(daV_conjuSub, CVarG1, dCVarG1);
d2aV_conjuMuSub = d2aV_conjuVec(1);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);

d2aV_MuVec = df_NatTwo_MuSet_CR_rmW(daV_MuSub, CVarG1, dCVarG1);
d2aV_MuconjMuSub = d2aV_MuVec(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Substitution on Chern_in_aV
variableSetG1 = [theta, aV, daV_u, daV_conju, d2aV_uu];
subSetG1 = [thetaSub, aVSub, daV_uSub, daV_conjuSub, d2aV_uuSub]; %latest

variableSetG2 = [daV_Mu, daV_conjMu, daV_vnormv, dtheta_Mu, dtheta_vnormv,...
    d2theta_MuconjMu];
subSetG2 = [daV_MuSub, daV_conjMuSub, daV_vnormvSub, dtheta_MuSub,...
    dtheta_vnormvSub, d2theta_MuconjMuSub]; %middle

variableSetG3 = [d2aV_uMu, d2aV_uconjMu, d2aV_conjuMu, d2aV_conjuconjMu,...
    d2aV_MuconjMu];
subSetG3 =  [d2aV_uMuSub, d2aV_uconjMuSub, d2aV_conjuMuSub, d2aV_conjuconjMuSub,...
    d2aV_MuconjMuSub]; %first

for j=1:16
    m = indexChern(j,1);
    n = indexChern(j,2);
    k = indexChern(j,3);
    ll = indexChern(j,4);
    temp = Chern_in_aV(m,n,k,ll);
    temp = subs(temp, variableSetG3, subSetG3);
    temp = subs(temp, variableSetG2, subSetG2);
    temp = subs(temp, variableSet0, subSet0);
    temp = subs(temp, variableSetG1, subSetG1);
    ChernTwo(m,n,k,ll) = temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Substitution

for j=1:16
    m = indexChern(j,1);
    n = indexChern(j,2);
    k = indexChern(j,3);
    ll = indexChern(j,4);
    temp = ChernTwo(m,n,k,ll); 
    temp = subs(temp, Y0, u*conj(u)+1);
    temp = subs(temp, d2G12_3N, d2G12_3Vec);
    temp = subs(temp, d2G23_1N, d2G23_1Vec);
    temp = subs(temp, d2G31_2N, d2G31_2Vec);
    temp = subs(temp, d2G11_2N, d2G11_2Vec);
    temp = subs(temp, d2G11_3N, d2G11_3Vec);
    temp = subs(temp, d2G22_1N, d2G22_1Vec);
    temp = subs(temp, d2G22_3N, d2G22_3Vec);
    temp = subs(temp, d2G33_1N, d2G33_1Vec);
    temp = subs(temp, d2G33_2N, d2G33_2Vec);
    temp = subs(temp, variableSetMu, subSetMu);
    temp = subs(temp, variableSetConjMu, subSetConjMu);
    temp = subs(temp, variableSetVnormv, subSetVnormv);
    temp = subs(temp, Y0, u*conj(u)+1);
    % Bianchi Identities
    temp = subs(temp, variableSetBianchi2, subSetBianchi2);
    temp = subs(temp, variableSetBianchi1, subSetBianchi1);
    temp = complex_simple3(temp,[w, dw_conjMu, dw_vnormv, dw_u,...
        d2w_uu, d2w_uconjMu, d2w_conjMuconjMu, u]);
    ChernTwo(m,n,k,ll) = temp;
end

clearvars j m n k ll temp 
clearvars subSetG1 subSetG2 subSetG3 variableSetG1 variableSetG2
clearvars variableSetG3 variableSet0 subSet0
clearvars dG12_3Row dG23_1Row dG31_2Row dG11_2Row dG11_3Row
clearvars dG22_1Row dG22_3Row dG33_1Row dG33_2Row
clearvars dG12_3ERow dG23_1ERow dG31_2ERow dG11_2ERow dG11_3ERow
clearvars dG22_1ERow dG22_3ERow dG33_1ERow dG33_2ERow
clearvars d2G12_3 d2G23_1 d2G31_2 d2G11_2 d2G11_3 d2G22_1 d2G22_3 d2G33_1 d2G33_2
clearvars daMVec daVVec dh11Vec dT4Vec d3aV_uconjuVec d4aV_uconjuMuVec 
clearvars d4aV_uconjuconjMuVec d4aV_uconjuconjuVec d4aV_uconjuvnormvVec

clearvars d2aV_MuVec d2aV_conjuVec d2aV_uVec d2h11_MuVec d2h11_uVec

list_of_variables_Chern2 = who;
save('DataChern2_NatTwoR_CR_rmW.mat');
% save('DataTempC2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

