% WeylCurv_Part52_NatTwo_CR_rmW.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Data_WeylChern_Part5_Jun10.mat');
assumeAlso([x y z u1 u2 gamma], 'real');

clearvars dG12_3Row d2G12_3 
clearvars dG23_1Row d2G23_1
clearvars dG31_2Row d2G31_2
clearvars dG11_2Row d2G11_2
clearvars dG11_3Row d2G11_3
clearvars dG22_1Row d2G22_1
clearvars dG22_3Row d2G22_3
clearvars dG33_1Row d2G33_1
clearvars dG33_2Row d2G33_2
clearvars d3w_uuRow d4w_uu dwRow 
clearvars dCVar3 dCVar4 d2w d2wRow d2rho

symSetWeylF = [];
% Put w = f in this context.
for j=1:120
    tempSet = symvar(WeylF(j,5));
    symSetWeylF = union(symSetWeylF, tempSet);
end
clear tempSet
% df = df_NatTwo_WeylPart5_CR_rmW(f, CVar5, dCVar5)

% f = -i/2*mu1*(mu3*G11_2 + mu1*G12_3 - mu2*G11_3)...
%     -i/2*mu2*(-mu3*G22_1 + mu1*G22_3 + mu2*G23_1)...
%     -i/2*mu3*(mu3*G31_2 - mu1*G33_2 + mu2*G33_1);
% dfdu = complexdiff3(f,u,0);
% d2fdu2 = complexdiff3(dfdu,u,0);
% phiF = d2fdu2 - 6*conj(u)/(1+u*conj(u))*dfdu + 12*conj(u)^2/(1+u*conj(u))^2*f;

dwN = [dw_Mu, dw_conjMu, dw_vnormv, dw_u];
dfVec0 = df_NatTwo_WeylPart5_CR_rmW(f, CVar5, dCVar5);
dfVec = [dfVec0(1), dfVec0(2), dfVec0(3), dfdu];
%
d2wN = [d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_uu,...
    d2w_MuMu, d2w_MuconjMu, d2w_Muvnormv, d2w_conjMuconjMu, d2w_conjMuvnormv,...
    d2w_vnormvvnormv];
d2f_uVec = df_NatTwo_WeylPart5_CR_rmW(dfdu, CVar5, dCVar5);
d2f_MuVec = df_NatTwo_WeylPart5_CR_rmW(dfVec0(1), CVar5, dCVar5);
d2f_conjMuVec = df_NatTwo_WeylPart5_CR_rmW(dfVec0(2), CVar5, dCVar5);
d2f_vnormvVec = df_NatTwo_WeylPart5_CR_rmW(dfVec0(3), CVar5, dCVar5);

d2f_uMu = d2f_uVec(1);
d2f_uconjMu = d2f_uVec(2);
d2f_uvnormv = d2f_uVec(3);
d2f_uu = d2fdu2;
d2f_MuMu = d2f_MuVec(1);
d2f_MuconjMu = d2f_MuVec(2);
d2f_Muvnormv = d2f_MuVec(3);
d2f_conjMuconjMu = d2f_conjMuVec(2);
d2f_conjMuvnormv = d2f_conjMuVec(3);
d2f_vnormvvnormv = d2f_vnormvVec(3);
clear d2fdu2

d2fVec = [d2f_uMu, d2f_uconjMu, d2f_uvnormv, d2f_uu,...
    d2f_MuMu, d2f_MuconjMu, d2f_Muvnormv, d2f_conjMuconjMu, d2f_conjMuvnormv,...
    d2f_vnormvvnormv];
%
d3wN = [d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv, d3w_uuu...
    d3w_uMuMu, d3w_uMuconjMu, d3w_uMuvnormv, d3w_uconjMuconjMu,...
    d3w_uconjMuvnormv];
d3f_uuVec = df_NatTwo_WeylPart5_CR_rmW(d2f_uu, CVar5, dCVar5);
d3f_uMuVec = df_NatTwo_WeylPart5_CR_rmW(d2f_uMu, CVar5, dCVar5);
d3f_uconjMuVec = df_NatTwo_WeylPart5_CR_rmW(d2f_uconjMu, CVar5, dCVar5);

d3fVec = [d3f_uuVec(1), d3f_uuVec(2), d3f_uuVec(3), d3f_uuVec(4),...
    d3f_uMuVec(1), d3f_uMuVec(2), d3f_uMuVec(3),...
    d3f_uconjMuVec(2), d3f_uconjMuVec(3)];
%
d4wN = [d4w_uuuMu, d4w_uuuconjMu, d4w_uuuu, ...
    d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuconjMuconjMu];
d4f_uuuVec = df_NatTwo_WeylPart5_CR_rmW(d3f_uuVec(4), CVar5, dCVar5);
d4f_uuMuVec = df_NatTwo_WeylPart5_CR_rmW(d3f_uuVec(1), CVar5, dCVar5);
d4f_uuconjMuVec = df_NatTwo_WeylPart5_CR_rmW(d3f_uuVec(2), CVar5, dCVar5);

d4fVec = [d4f_uuuVec(1), d4f_uuuVec(2), d4f_uuuVec(4),...
    d4f_uuMuVec(1), d4f_uuMuVec(2), d4f_uuconjMuVec(2)];
%
MVarFive2 = [u, G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3, G33_1, G33_2,...
    dG12_3_Mu, dG12_3_conjMu, dG12_3_vnormv,...
    dG23_1_Mu, dG23_1_conjMu, dG23_1_vnormv,...
    dG31_2_Mu, dG31_2_conjMu, dG31_2_vnormv,...
    dG11_2_Mu, dG11_2_conjMu, dG11_2_vnormv,...
    dG11_3_Mu, dG11_3_conjMu, dG11_3_vnormv,...
    dG22_1_Mu, dG22_1_conjMu, dG22_1_vnormv,...
    dG22_3_Mu, dG22_3_conjMu, dG22_3_vnormv,...
    dG33_1_Mu, dG33_1_conjMu, dG33_1_vnormv,...
    dG33_2_Mu, dG33_2_conjMu, dG33_2_vnormv,...
    d2G12_3_MuMu, d2G12_3_MuconjMu, d2G12_3_Muvnormv,...
    d2G12_3_conjMuconjMu, d2G12_3_conjMuvnormv, d2G12_3_vnormvvnormv,...
    d2G23_1_MuMu, d2G23_1_MuconjMu, d2G23_1_Muvnormv,...
    d2G23_1_conjMuconjMu, d2G23_1_conjMuvnormv, d2G23_1_vnormvvnormv,...
    d2G31_2_MuMu, d2G31_2_MuconjMu, d2G31_2_Muvnormv,...
    d2G31_2_conjMuconjMu, d2G31_2_conjMuvnormv, d2G31_2_vnormvvnormv,...
    d2G11_2_MuMu, d2G11_2_MuconjMu, d2G11_2_Muvnormv,...
    d2G11_2_conjMuconjMu, d2G11_2_conjMuvnormv, d2G11_2_vnormvvnormv,...
    d2G11_3_MuMu, d2G11_3_MuconjMu, d2G11_3_Muvnormv,...
    d2G11_3_conjMuconjMu, d2G11_3_conjMuvnormv, d2G11_3_vnormvvnormv,...
    d2G22_1_MuMu, d2G22_1_MuconjMu, d2G22_1_Muvnormv,...
    d2G22_1_conjMuconjMu, d2G22_1_conjMuvnormv, d2G22_1_vnormvvnormv,...
    d2G22_3_MuMu, d2G22_3_MuconjMu, d2G22_3_Muvnormv,...
    d2G22_3_conjMuconjMu, d2G22_3_conjMuvnormv, d2G22_3_vnormvvnormv,...
    d2G33_1_MuMu, d2G33_1_MuconjMu, d2G33_1_Muvnormv,...
    d2G33_1_conjMuconjMu, d2G33_1_conjMuvnormv, d2G33_1_vnormvvnormv,...
    d2G33_2_MuMu, d2G33_2_MuconjMu, d2G33_2_Muvnormv,...
    d2G33_2_conjMuconjMu, d2G33_2_conjMuvnormv, d2G33_2_vnormvvnormv];
%
WeylFf = sym('WeylFf',[120,5]);
for j=1:120
    m = countSetTwo(j,1);
    n = countSetTwo(j,2);
    k = countSetTwo(j,3);
    ll = countSetTwo(j,4);
    temp = WeylF(j,5);
    temp = subs(temp, w, f);
    temp = subs(temp, dwN, dfVec);
    temp = subs(temp, d2wN, d2fVec);
    temp = subs(temp, d3wN, d3fVec);
    temp = subs(temp, d4wN, d4fVec);
%     temp = complex_simple3(temp, MVarFive2);
    WeylFf(j,:) = [m, n, k, ll, temp];
end

symSetWf = [];
for j=1:120
    tempSet = symvar(WeylFf(j,5));
    symSetWf = union(symSetWf, tempSet);
end
clearvars j m n k ll temp tempSet

save('DataZeroV2_WeylChern_Part52_Jun15.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


