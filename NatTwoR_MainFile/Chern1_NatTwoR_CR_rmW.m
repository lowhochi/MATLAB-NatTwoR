% Important Update. Oct 9 2019.

load('DataMain1Ch_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Replace variables in symSetChern by aV and its derivatives.

% variable = [T4, aM, aV, bV, h11,... 
%     d2h11_uu, d2h11_uMu, d2h11_uconju, d2h11_conjuMu, d2h11_uconjMu,...
%     d2h11_MuconjMu, d2h11_conjuconju, d2h11_conjuconjMu,...
%     dT4_u, dT4_conju, dT4_conjMu, daM_u, daM_Mu, daM_conju, daM_conjMu,...
%     dh11_u, dh11_Mu, dh11_conju, dh11_conjMu ,dh11_vnormv];
% wSet =  [d2w_uu, d2w_uconjMu, d2w_conjMuconjMu, dw_u, dw_conjMu,...
%     dw_vnormv, w];

variable_NatTwoR_Chern1_CR_rmW

indexChern = [1,1,1,1; 1,1,1,2; 1,1,2,1; 1,1,2,2;
    1,2,1,1; 1,2,1,2; 1,2,2,1; 1,2,2,2;
    2,1,1,1; 2,1,1,2; 2,1,2,1; 2,1,2,2;
    2,2,1,1; 2,2,1,2; 2,2,2,1; 2,2,2,2];

aMSub = -(1+u*conj(u))^2*daV_u;
T4Sub = ((1+u*conj(u))^2/2)*daV_conju;
bVSub = u*(1+u*conj(u))*conj(daV_u) + (1+u*conj(u))^2/2*conj(d2aV_uu);
h11Sub = i*(1+u*conj(u))^2*aV + (i/2)*(1+u*conj(u))^4*d2aV_uconju;

% df = df_NatTwo_MuSet_CR_rmW(f, CVarChern1, dCVarChern1);
dh11Vec = df_NatTwo_MuSet_CR_rmW(h11Sub, CVarChern1, dCVarChern1);
dh11_MuSub = dh11Vec(1);
dh11_conjMuSub = dh11Vec(2);
dh11_vnormvSub = dh11Vec(3);
dh11_uSub = dh11Vec(4);
dh11_conjuSub = dh11Vec(5);

daMVec =  df_NatTwo_MuSet_CR_rmW(aMSub, CVarChern1, dCVarChern1);
daM_MuSub = daMVec(1);
daM_conjMuSub = daMVec(2);
daM_uSub = daMVec(4);
daM_conjuSub = daMVec(5);

dT4Vec = df_NatTwo_MuSet_CR_rmW(T4Sub, CVarChern1, dCVarChern1);
dT4_conjMuSub = dT4Vec(2);
dT4_uSub = dT4Vec(4);
dT4_conjuSub = dT4Vec(5);

d2h11_MuVec = df_NatTwo_MuSet_CR_rmW(dh11_MuSub, CVarChern1, dCVarChern1);
d2h11_MuconjMuSub = d2h11_MuVec(2);

d2h11_uVec = df_NatTwo_MuSet_CR_rmW(dh11_uSub, CVarChern1, dCVarChern1);
d2h11_uuSub = d2h11_uVec(4);
d2h11_uMuSub = d2h11_uVec(1);
d2h11_uconjMuSub = d2h11_uVec(2);
d2h11_uconjuSub = d2h11_uVec(5);
d2h11_conjuMuSub = conj(d2h11_uconjMuSub);
d2h11_conjuconjuSub = conj(d2h11_uuSub);
d2h11_conjuconjMuSub = conj(d2h11_uMuSub);

variableSet1 = [aM, bV, T4, h11]; %last
variableSet2 = [daM_Mu, daM_conjMu, daM_u, daM_conju,...
    dT4_conjMu, dT4_u, dT4_conju,...
    dh11_Mu, dh11_conjMu, dh11_vnormv, dh11_u, dh11_conju]; %second
variableSet3 = [d2h11_uu, d2h11_uconju, d2h11_conjuconju, d2h11_MuconjMu,...
    d2h11_uMu, d2h11_uconjMu, d2h11_conjuMu, d2h11_conjuconjMu]; %first

subSet1 = [aMSub, bVSub, T4Sub, h11Sub];
subSet2 = [daM_MuSub, daM_conjMuSub, daM_uSub, daM_conjuSub,...
    dT4_conjMuSub, dT4_uSub, dT4_conjuSub,...
    dh11_MuSub, dh11_conjMuSub, dh11_vnormvSub, dh11_uSub, dh11_conjuSub];
subSet3 = [d2h11_uuSub, d2h11_uconjuSub, d2h11_conjuconjuSub, d2h11_MuconjMuSub,...
    d2h11_uMuSub, d2h11_uconjMuSub, d2h11_conjuMuSub, d2h11_conjuconjMuSub];

Chern_in_aV = sym('Chern_in_aV',[2 2 2 2]);
for j=1:16
    m = indexChern(j,1);
    n = indexChern(j,2);
    k = indexChern(j,3);
    ll = indexChern(j,4);
    temp = Chern(m,n,k,ll);
    temp = subs(temp, variableSet3, subSet3);
    temp = subs(temp, variableSet2, subSet2);
    temp = subs(temp, variableSet1, subSet1);
    temp = complex_simple3(temp, MVarChern1);
    Chern_in_aV(m,n,k,ll) = temp;
end

save('DataChern1_Part1_NatTwoR_CR_rmW.mat');


