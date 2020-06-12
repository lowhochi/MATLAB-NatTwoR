% check1_NatTwoR_CR_rmW.m
syms u
syms G12_3 G23_1 G31_2 G11_2 G11_3 G22_1 G22_3 G33_1 G33_2 real
syms dG12_3_Mu dG12_3_conjMu dG12_3_vnormv
syms dG23_1_Mu dG23_1_conjMu dG23_1_vnormv
syms dG31_2_Mu dG31_2_conjMu dG31_2_vnormv 
syms dG11_2_Mu dG11_2_conjMu dG11_2_vnormv 
syms dG11_3_Mu dG11_3_conjMu dG11_3_vnormv 
syms dG22_1_Mu dG22_1_conjMu dG22_1_vnormv 
syms dG22_3_Mu dG22_3_conjMu dG22_3_vnormv
syms dG33_1_Mu dG33_1_conjMu dG33_1_vnormv
syms dG33_2_Mu dG33_2_conjMu dG33_2_vnormv
y0 = 1+u*conj(u);
mu1 = u^2-1;
mu2 = 2*u;
mu3 = i*(u^2+1);
norm_square_of_mu = mu1*conj(mu1) +mu2*conj(mu2) +mu3*conj(mu3);
v1 = i*(mu2*conj(mu3)-mu3*conj(mu2));
v2 = i*(mu3*conj(mu1)-mu1*conj(mu3));
v3 = i*(mu1*conj(mu2)-mu2*conj(mu1));
norm_square_of_v = v1*v1 +v2*v2 +v3*v3;
norm_of_v = norm_square_of_v^(1/2);
v1normv = v1/norm_of_v;
v2normv = v2/norm_of_v;
v3normv = v3/norm_of_v;
v1normv = complex_simple3(v1normv,u);
v2normv = complex_simple3(v2normv,u);
v3normv = complex_simple3(v3normv,u);

dmu_du1 = (2*conj(u)/y0)*mu1 + 2*v1normv;
dmu_du2 = (2*conj(u)/y0)*mu2 + 2*v2normv;
dmu_du3 = (2*conj(u)/y0)*mu3 + 2*v3normv;
dmu_du1 = complex_simple3(dmu_du1, [u]);
dmu_du2 = complex_simple3(dmu_du2, [u]);
dmu_du3 = complex_simple3(dmu_du3, [u]);

lie_e1_e2 = [-G11_2; G22_1; G12_3+G23_1];
lie_e2_e3 = [G23_1+G31_2; -G22_3; G33_2];
lie_e1_e3 = [-G11_3; -G12_3-G31_2; G33_1];

lie_mu_conjmu = (mu1*conj(mu2)-mu2*conj(mu1))*lie_e1_e2...
    +(mu2*conj(mu3)-mu3*conj(mu2))*lie_e2_e3...
    +(mu1*conj(mu3)-mu3*conj(mu1))*lie_e1_e3;
lie_mu_vnormv = (mu1*v2normv -mu2*v1normv)*lie_e1_e2...
    +(mu2*v3normv -mu3*v2normv)*lie_e2_e3...
    +(mu1*v3normv -mu3*v1normv)*lie_e1_e3;

CVar = [u, G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3, G33_1, G33_2];
derivativeDict.u = [0;0;0;1;0];
derivativeDict.G12_3 = [dG12_3_Mu; dG12_3_conjMu; dG12_3_vnormv; 0; 0];
derivativeDict.G23_1 = [dG23_1_Mu; dG23_1_conjMu; dG23_1_vnormv; 0; 0];
derivativeDict.G31_2 = [dG31_2_Mu; dG31_2_conjMu; dG31_2_vnormv; 0; 0];
derivativeDict.G11_2 = [dG11_2_Mu; dG11_2_conjMu; dG11_2_vnormv; 0; 0];
derivativeDict.G11_3 = [dG11_3_Mu; dG11_3_conjMu; dG11_3_vnormv; 0; 0];
derivativeDict.G22_1 = [dG22_1_Mu; dG22_1_conjMu; dG22_1_vnormv; 0; 0];
derivativeDict.G22_3 = [dG22_3_Mu; dG22_3_conjMu; dG22_3_vnormv; 0; 0];
derivativeDict.G33_1 = [dG33_1_Mu; dG33_1_conjMu; dG33_1_vnormv; 0; 0];
derivativeDict.G33_2 = [dG33_2_Mu; dG33_2_conjMu; dG33_2_vnormv; 0; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aVtest = i/(2*y0^2)*( (mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    + (mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    + (mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    -mu3*conj(mu1)*G11_2 +mu2*conj(mu1)*G11_3 +mu3*conj(mu2)*G22_1...
    -mu1*conj(mu2)*G22_3 -mu2*conj(mu3)*G33_1 + mu1*conj(mu3)*G33_2);

aV = (1/norm_square_of_mu)*(conj(mu1)*lie_mu_vnormv(1)...
    +conj(mu2)*lie_mu_vnormv(2) +conj(mu3)*lie_mu_vnormv(3));
daVVec = df_check1_NatTwoR_CR_rmW(aV, CVar, derivativeDict);
% difference_aV = aV - aVtest;
% difference_aV = complex_simple3(difference_aV, [u]);

aM = (1/norm_square_of_mu)*(conj(mu1)*lie_mu_conjmu(1)...
    +conj(mu2)*lie_mu_conjmu(2) +conj(mu3)*lie_mu_conjmu(3));
aMtest = -y0^2*daVVec(4);
% difference_aM = aM - aMtest;
% difference_aM = complex_simple3(difference_aM, [u]);
daV_u = complex_simple3(daVVec(4), u);
daV_conju = complex_simple3(daVVec(5),u);
d2aV_uVec = df_check1_NatTwoR_CR_rmW(daV_u, CVar, derivativeDict);
d2aV_uu = d2aV_uVec(4);
d2aV_uconju = d2aV_uVec(5);

bV = (1/norm_square_of_mu)*(mu1*lie_mu_vnormv(1)...
    +mu2*lie_mu_vnormv(2) +mu3*lie_mu_vnormv(3));
bVtest = u*y0*conj(daV_u) +y0^2/2*conj(d2aV_uu); 
% difference_bV = bV- bVtest;
% difference_bV = complex_simple3(difference_bV, u);

h11 = -i/2*(v1normv*lie_mu_conjmu(1) + v2normv*lie_mu_conjmu(2)...
    + v3normv*lie_mu_conjmu(3));
% h11test = i*y0^2*aV + i*y0^4/2*d2aV_uconju;
% difference_h11 = complex_simple3(h11-h11test, u);

T4 = 1/2*(v1normv*lie_mu_vnormv(1) +v2normv*lie_mu_vnormv(2)...
    +v3normv*lie_mu_vnormv(3));
% T4test = y0^2/2*daV_conju;
% difference_T4 = complex_simple3(T4-T4test, u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = G12_3 + G23_1 + G31_2;
d2aV_conjuVec = df_check1_NatTwoR_CR_rmW(daV_conju, CVar, derivativeDict);
d2aV_conjuconju = d2aV_conjuVec(5);
d3aV_uuVec = df_check1_NatTwoR_CR_rmW(d2aV_uu, CVar, derivativeDict);
d3aV_uuu = d3aV_uuVec(4);

d2aV_uconjuTest = 2/y0^2*(2*i*theta+conj(aV)-2*aV);
difference_d2aV_uconju = d2aV_uconju - d2aV_uconjuTest;
difference_d2aV_uconju = complex_simple3(difference_d2aV_uconju,u);

d2aV_conjuconjuTest = (-2*u/y0)*(daV_conju +conj(daV_u))-conj(d2aV_uu);
difference_d2aV_conjuconju = d2aV_conjuconju - d2aV_conjuconjuTest; 
difference_d2aV_conjuconju = complex_simple3(difference_d2aV_conjuconju, u);

d3aV_uuuTest = -(6*conj(u)^2/y0^2)*daV_u -(6*conj(u)/y0)*d2aV_uu;
difference_d3aV_uuu = d3aV_uuu - d3aV_uuuTest;
difference_d3aV_uuu = complex_simple3(difference_d3aV_uuu, u);

f = -i/2*(mu1*(G11_2*mu3 +G12_3*mu1 -G11_3*mu2)...
    +mu2*(-G22_1*mu3 + G22_3*mu1 +G23_1*mu2)...
    +mu3*(G31_2*mu3 -G33_2*mu1 +G33_1*mu2) );
fTest = -u*y0^3*daV_conju -y0^4/2*d2aV_conjuconju;
fTest2 = u*y0^3*conj(daV_u) +y0^4/2*conj(d2aV_uu);
difference_f = complex_simple3(f-fTest, u);
difference_f2 = complex_simple3(f-fTest2, u);

h11Test2 = y0^2*(-2*theta +i*conj(aV)-i*aV);
difference_h112 = complex_simple3(h11-h11Test2, u);

dfVec = df_check1_NatTwoR_CR_rmW(f, CVar, derivativeDict);
df_u = complex_simple3(dfVec(4), u);
d2f_uVec = df_check1_NatTwoR_CR_rmW(df_u, CVar, derivativeDict);
d2f_uu = complex_simple3(d2f_uVec(4), u);
phiF = d2f_uu -6*conj(u)/y0*df_u + 12*conj(u)^2/y0^2*f;
phiFtest = 8*i*theta -6*aV +6*conj(aV);
difference_phiF = complex_simple3(phiF - phiFtest, u);

rhoF = i*phiF - i*conj(phiF);
rhoFtest = -16*theta -12*i*(aV-conj(aV));
difference_rhoF = complex_simple3(rhoF-rhoFtest, u);






