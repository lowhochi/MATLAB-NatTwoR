% RicCurv_NatTwoR_CR_rmW.m
load('DataMain4_NatTwoR_CR_rmW.mat');
syms theta % theta = G12_3 + G23_1 + G31_2;

variable_NatTwoR_RicCurv_CR_rmW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RicCurv(m,n) = Ric(u_m,u_n)
% Analyze the Ricci curvature tensor.
% (1) Replace Ricci terms by aV and it derivatives.
symSetRic = symvar(RicCurv);
wSet = [u, w, dw_Mu, dw_conjMu, dw_u, dw_vnormv,...
    d2w_uMu, d2w_uconjMu, d2w_uu, d2w_uvnormv,...
    d3w_uuMu, d3w_uuconjMu, d3w_uuu, d3w_uuvnormv]; 
    
variableSet1 = [aM, bV, T4, h11, rho];

variableSet2 = [daM_Mu, daM_conjMu, daM_vnormv, daM_u, daM_conju,...
    dT4_Mu, dT4_conjMu, dT4_vnormv, dT4_u, dT4_conju,...
    dbV_conjMu, dbV_u, dbV_conju,...
    dh11_Mu, dh11_conjMu, dh11_vnormv, dh11_u, dh11_conju,...
    drho_Mu, drho_conjMu, drho_vnormv, drho_u, drho_conju];

variableSet3 = [d2aM_conjuMu, d2aM_conjuconjMu, d2aM_conjuconju,...
    d2aM_conjuvnormv, d2aM_uMu, d2aM_uconju, d2aM_uu,...
    d2h11_conjuMu, d2h11_conjuconjMu, d2h11_conjuconju,...
    d2h11_conjuvnormv, d2h11_uMu, d2h11_uconjMu, d2h11_uconju,...
    d2h11_uu, d2h11_uvnormv, d2rho_conjuconjMu, d2rho_conjuconju,...
    d2rho_uMu, d2rho_uconju, d2rho_uu];

variableSet4 = [d2T4_conjuconjMu, d2T4_conjuconju, d2T4_uMu, d2T4_uconjMu, ...
    d2T4_uconju, d2T4_uu, d2T4_uvnormv];

variableSet5 = [d3T4_uconjuconjMu, d3T4_uconjuconju, d3T4_uuMu,...
    d3T4_uuconju, d3T4_uuu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet1 = [aM, bV, T4, h11, rho];
aMSub = -(1+u*conj(u))^2*daV_u;
T4Sub = ((1+u*conj(u))^2/2)*daV_conju;
bVSub = u*(1+u*conj(u))*conj(daV_u) + (1+u*conj(u))^2/2*conj(d2aV_uu);
h11Sub = i*(1+u*conj(u))^2*aV + (i/2)*(1+u*conj(u))^4*d2aV_uconju;

phiW = d2w_uu - 6*conj(u)/(1+u*conj(u))*dw_u + 12*conj(u)^2/(1+u*conj(u))^2*w;
rhoSub = i*(phiW-conj(phiW)) + 16*theta + 12*i*(aV-conj(aV));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet2 = [daM_Mu, daM_conjMu, daM_vnormv, daM_u, daM_conju,...
%     dT4_Mu, dT4_conjMu, dT4_vnormv, dT4_u, dT4_conju,...
%     dbV_conjMu, dbV_u, dbV_conju,...
%     dh11_Mu, dh11_conjMu, dh11_vnormv, dh11_u, dh11_conju,...
%     drho_Mu, drho_conjMu, drho_vnormv, drho_u, drho_conju];

daMVec = df_NatTwo_MuSet_CR_rmW(aMSub, CVarRicCurv, dCVarRicCurv);
daM_MuSub = daMVec(1);
daM_conjMuSub = daMVec(2);
daM_vnormvSub = daMVec(3);
daM_uSub = daMVec(4);
daM_conjuSub = daMVec(5);

dT4Vec = df_NatTwo_MuSet_CR_rmW(T4Sub, CVarRicCurv, dCVarRicCurv);
dT4_MuSub = dT4Vec(1);
dT4_conjMuSub = dT4Vec(2);
dT4_vnormvSub = dT4Vec(3);
dT4_uSub = dT4Vec(4);
dT4_conjuSub = dT4Vec(5);

dbVVec = df_NatTwo_MuSet_CR_rmW(bVSub, CVarRicCurv, dCVarRicCurv);
dbV_conjMuSub = dbVVec(2);
dbV_uSub = dbVVec(4);
dbV_conjuSub = dbVVec(5);

dh11Vec = df_NatTwo_MuSet_CR_rmW(h11Sub, CVarRicCurv, dCVarRicCurv);
dh11_MuSub = dh11Vec(1);
dh11_conjMuSub = dh11Vec(2);
dh11_vnormvSub = dh11Vec(3);
dh11_uSub = dh11Vec(4);
dh11_conjuSub = dh11Vec(5);

drhoVec = df_NatTwo_MuSet_CR_rmW(rhoSub, CVarRicCurv, dCVarRicCurv);
drho_MuSub = drhoVec(1);
drho_conjMuSub = drhoVec(2);
drho_vnormvSub = drhoVec(3);
drho_uSub = drhoVec(4);
drho_conjuSub = drhoVec(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet3 = [d2aM_conjuMu, d2aM_conjuconjMu, d2aM_conjuconju,...
%     d2aM_conjuvnormv, d2aM_uMu, d2aM_uconju, d2aM_uu,...
%     d2h11_conjuMu, d2h11_conjuconjMu, d2h11_conjuconju,...
%     d2h11_conjuvnormv, d2h11_uMu, d2h11_uconjMu, d2h11_uconju,...
%     d2h11_uu, d2h11_uvnormv, d2rho_conjuconjMu, d2rho_conjuconju,...
%     d2rho_uMu, d2rho_uconju, d2rho_uu];

diffVariable2 = symvar([daM_uSub, daM_conjuSub, dT4_uSub, dT4_conjuSub,...
    dh11_uSub, dh11_conjuSub, drho_uSub, drho_conjuSub]);

d2aM_uVec = df_NatTwo_MuSet_CR_rmW(daM_uSub, CVarRicCurv, dCVarRicCurv);
d2aM_uMuSub = d2aM_uVec(1);
d2aM_uuSub = d2aM_uVec(4);
d2aM_uconjuSub = d2aM_uVec(5);

d2aM_conjuVec = df_NatTwo_MuSet_CR_rmW(daM_conjuSub, CVarRicCurv, dCVarRicCurv);
d2aM_conjuMuSub = d2aM_conjuVec(1);
d2aM_conjuconjMuSub = d2aM_conjuVec(2);
d2aM_conjuvnormvSub = d2aM_conjuVec(3);
d2aM_conjuconjuSub = d2aM_conjuVec(5);

d2h11_uVec = df_NatTwo_MuSet_CR_rmW(dh11_uSub, CVarRicCurv, dCVarRicCurv);
d2h11_uMuSub = d2h11_uVec(1);
d2h11_uconjMuSub = d2h11_uVec(2);
d2h11_uvnormvSub = d2h11_uVec(3);
d2h11_uuSub = d2h11_uVec(4);
d2h11_uconjuSub = d2h11_uVec(5);

d2h11_conjuVec = df_NatTwo_MuSet_CR_rmW(dh11_conjuSub, CVarRicCurv, dCVarRicCurv);
d2h11_conjuMuSub = d2h11_conjuVec(1);
d2h11_conjuconjMuSub = d2h11_conjuVec(2);
d2h11_conjuvnormvSub = d2h11_conjuVec(3);
d2h11_conjuconjuSub = d2h11_conjuVec(5);

d2rho_uVec = df_NatTwo_MuSet_CR_rmW(drho_uSub, CVarRicCurv, dCVarRicCurv);
d2rho_conjuVec = df_NatTwo_MuSet_CR_rmW(drho_conjuSub, CVarRicCurv, dCVarRicCurv);

d2rho_conjuconjMuSub = d2rho_conjuVec(2);
d2rho_conjuconjuSub = d2rho_conjuVec(5); % updated 
d2rho_uMuSub = d2rho_uVec(1);
d2rho_uuSub = d2rho_uVec(4);
d2rho_uconjuSub = d2rho_uVec(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet4 = [d2T4_conjuconjMu, d2T4_conjuconju, d2T4_uMu, d2T4_uconjMu, ...
%     d2T4_uconju, d2T4_uu, d2T4_uvnormv];

d2T4_uVec = df_NatTwo_MuSet_CR_rmW(dT4_uSub, CVarRicCurv, dCVarRicCurv);
d2T4_uMuSub = d2T4_uVec(1);
d2T4_uconjMuSub = d2T4_uVec(2);
d2T4_uvnormvSub = d2T4_uVec(3);
d2T4_uuSub = d2T4_uVec(4);
d2T4_uconjuSub = d2T4_uVec(5);

d2T4_conjuVec = df_NatTwo_MuSet_CR_rmW(dT4_conjuSub, CVarRicCurv, dCVarRicCurv);
d2T4_conjuconjMuSub = d2T4_conjuVec(2);
d2T4_conjuconjuSub = d2T4_conjuVec(5);

% variableSet5 = [d3T4_uconjuconjMu, d3T4_uconjuconju, d3T4_uuMu,...
%     d3T4_uuconju, d3T4_uuu];

diffVariable4 = symvar([d2T4_uuSub, d2T4_uconjuSub]);
d3T4_uuVec = df_NatTwo_MuSet_CR_rmW(d2T4_uuSub, CVarRicCurv, dCVarRicCurv);
d3T4_uconjuVec = df_NatTwo_MuSet_CR_rmW(d2T4_uconjuSub, CVarRicCurv, dCVarRicCurv);

d3T4_uuMuSub = d3T4_uuVec(1);
d3T4_uuuSub = d3T4_uuVec(4);
d3T4_uuconjuSub = d3T4_uuVec(5);
d3T4_uconjuconjMuSub = d3T4_uconjuVec(2);
d3T4_uconjuconjuSub = d3T4_uconjuVec(5);

subSet1 = [aMSub, bVSub, T4Sub, h11Sub, rhoSub];

subSet2 = [daM_MuSub, daM_conjMuSub, daM_vnormvSub, daM_uSub, daM_conjuSub,...
    dT4_MuSub, dT4_conjMuSub, dT4_vnormvSub, dT4_uSub, dT4_conjuSub,...
    dbV_conjMuSub, dbV_uSub, dbV_conjuSub,...
    dh11_MuSub, dh11_conjMuSub, dh11_vnormvSub, dh11_uSub, dh11_conjuSub,...
    drho_MuSub, drho_conjMuSub, drho_vnormvSub, drho_uSub, drho_conjuSub];

subSet3 = [d2aM_conjuMuSub, d2aM_conjuconjMuSub, d2aM_conjuconjuSub,...
    d2aM_conjuvnormvSub, d2aM_uMuSub, d2aM_uconjuSub, d2aM_uuSub,...
    d2h11_conjuMuSub, d2h11_conjuconjMuSub, d2h11_conjuconjuSub,...
    d2h11_conjuvnormvSub, d2h11_uMuSub, d2h11_uconjMuSub, d2h11_uconjuSub,...
    d2h11_uuSub, d2h11_uvnormvSub, d2rho_conjuconjMuSub, d2rho_conjuconjuSub,...
    d2rho_uMuSub, d2rho_uconjuSub, d2rho_uuSub];

subSet4 = [d2T4_conjuconjMuSub, d2T4_conjuconjuSub, d2T4_uMuSub, d2T4_uconjMuSub, ...
    d2T4_uconjuSub, d2T4_uuSub, d2T4_uvnormvSub];

subSet5 = [d3T4_uconjuconjMuSub, d3T4_uconjuconjuSub, d3T4_uuMuSub,...
    d3T4_uuconjuSub, d3T4_uuuSub];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Substitution
RicCurv_in_aV = sym('Ric',[6 6]);

MVarRicCurv1 = [u, w, dw_Mu, dw_conjMu, dw_u, dw_vnormv,...
    d2w_uMu, d2w_uconjMu, d2w_uu, d2w_uvnormv,...
    d3w_uuMu, d3w_uuconjMu, d3w_uuu, d3w_uuvnormv,...
    d4w_uuuMu, d4w_uuuconjMu, d4w_uuuvnormv, d4w_uuuu,...
    theta, dtheta_Mu dtheta_conjMu dtheta_vnormv];

MVarRicCurv2 = [aV, daV_Mu, daV_conjMu, daV_vnormv, daV_u, daV_conju,...
    d2aV_MuMu, d2aV_MuconjMu, d2aV_Muvnormv, d2aV_conjMuconjMu,...
    d2aV_conjMuvnormv, d2aV_vnormvvnormv, d2aV_uMu, d2aV_uconjMu,...
    d2aV_uvnormv, d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv,...
    d2aV_uu, d2aV_uconju, d2aV_conjuconju,...
    d3aV_uconjuconju, d3aV_uuconju, d3aV_uuu,...
    d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv, d3aV_uconjMuconjMu,...
    d3aV_uconjMuvnormv, d3aV_uvnormvvnormv, d3aV_uuMu, d3aV_uuconjMu,...
    d3aV_uuvnormv, d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv,...
    d3aV_conjuconjuMu, d3aV_conjuconjuconjMu, d3aV_conjuconjuvnormv,...
    d3aV_conjuconjuconju, d4aV_uconjuMuMu, d4aV_uconjuMuconjMu, d4aV_uconjuMuvnormv,...
    d4aV_uconjuconjMuconjMu, d4aV_uconjuconjMuvnormv, d4aV_uconjuvnormvvnormv,...
    d4aV_uuconjuMu, d4aV_uuconjuconjMu, d4aV_uuconjuvnormv,...
    d4aV_uconjuconjuMu, d4aV_uconjuconjuconjMu, d4aV_uconjuconjuvnormv,...
    d4aV_uuuconju, d4aV_uuconjuconju, d4aV_uconjuconjuconju];
 
MVarRicCurv = cat(2, MVarRicCurv1, MVarRicCurv2);

for m=1:6
    for n=1:6
        temp = RicCurv(m,n);
        temp = subs(temp, variableSet5, subSet5);
        temp = subs(temp, variableSet4, subSet4);
        temp = subs(temp, variableSet3, subSet3);
        temp = subs(temp, variableSet2, subSet2);
        temp = subs(temp, variableSet1, subSet1);
        temp = complex_simple3(temp, MVarRicCurv);
        RicCurv_in_aV(m,n) = temp;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars m n temp myNumber

save('DataRicCurv1_Oct7_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%