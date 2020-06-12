% Obsolete

% b/f WeylCurv_Part4_NatTwo_CR_rmW.m
% M is flat 
% and w = lambda_0+lambda_1*u+K*u^2-conj(lambda_1)*u^3+conj(lambda_0)*u^4. 
load('Data_WeylChern_Part4_May30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
clearvars countSet countSetTwo Chern ChernX0 dCVar3
clearvars WeylX0 WeylX WCdiffX0
clearvars symSetChern symSetWeylX indexChern indexWeyl
clearvars MVarX0 MVarS
% simplify set: symSetWeylS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countLinear = [1, 2, 1, 4;
    1, 2, 2, 3;
    1, 2, 3, 5;
    1, 2, 4, 5;
    1, 4, 1, 5;
    1, 4, 2, 5;
    1, 5, 2, 3;
    1, 5, 3, 5;
    1, 5, 4, 5;
    2, 3, 2, 5;
    2, 5, 3, 5;
    2, 5, 4, 5];
countLinearRef = [3, 6, 11, 13, 31, 35, 45, 50, 52, 68, 88, 90];

countNonLinear = [1, 2, 1, 2;
    1, 2, 1, 5;
    1, 2, 2, 5;
    1, 5, 1, 5;
    1, 5, 2, 5;
    2, 5, 2, 5];

WeylSLin = sym('WeylSlin',[12 5]);
symSetLin = []; % length = 13;
for j=1:12
    mm = countLinear(j,1);
    nn = countLinear(j,2);
    kk = countLinear(j,3);
    ll = countLinear(j,4);
    N = countLinearRef(j);
    temp = WeylSS(N,5);
    tempSet = symvar(temp);
    symSetLin = union(symSetLin, tempSet);
    WeylSLin(j,:) = [mm,nn,kk,ll,temp];
    clearvars mm nn kk ll temp tempSet N
end
v1normv = complex_simple3(v1/norm_of_v,[u]);
v2normv = complex_simple3(v2/norm_of_v,[u]);
v3normv = complex_simple3(v3/norm_of_v,[u]);
muVec = [mu1;mu2;mu3];
conjMuVec = [conj(mu1);conj(mu2);conj(mu3)];
vnormvVec = [v1normv;v2normv;v3normv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psi = Psi(p) holomorphic in p=x+i*z
% Sigma = Sigma(p,conj(p))
% H0 = H0(p,conj(p),y)
syms p
syms Psi dPsi_p d2Psi_pp
syms Sigma dSigma_p dSigma_conjp 
syms d2Sigma_pp d2Sigma_pconjp 
syms H0 dH0_y real
syms dH0_p 

d2Sigma_conjpconjp = conj(dH0_p) + conj(d2Sigma_pconjp);

MVarPsi2 = [u, Psi, Sigma, H0, dPsi_p, d2Psi_pp,...
    dSigma_p, dSigma_conjp, d2Sigma_pp, d2Sigma_pconjp, dH0_y, dH0_p];
lambda1Sub = Psi;
lambda0Sub = -y/2*dPsi_p + dSigma_p;
KSub = -dSigma_conjp - conj(dSigma_conjp) + H0;

dlambda1_MuSub = -2*dPsi_p;
dlambda1_conjMuSub = 2*conj(u)^2*dPsi_p;
dlambda1_vnormvSub = 2*conj(u)*dPsi_p/(1+u*conj(u));

dlambda0_MuSub = y*d2Psi_pp - u*dPsi_p - 2*d2Sigma_pp + 2*u^2*d2Sigma_pconjp;

dlambda0_conjMuSub = -y*conj(u)^2*d2Psi_pp - conj(u)*dPsi_p...
    + 2*conj(u)^2*d2Sigma_pp - 2*d2Sigma_pconjp;

dlambda0_vnormvSub = -conj(u)/(1+u*conj(u))*y*d2Psi_pp...
    + (u*conj(u)-1)/(2*(1+u*conj(u)))*dPsi_p...
    + 2*conj(u)/(1+u*conj(u))*d2Sigma_pp + 2*u/(1+u*conj(u))*d2Sigma_pconjp;

dK_MuSub = 2*d2Sigma_pconjp - 2*u^2*d2Sigma_conjpconjp...
    + 2*conj(d2Sigma_conjpconjp) - 2*u^2*conj(d2Sigma_pconjp)...
    - 2*dH0_p + 2*u^2*conj(dH0_p) + 2*u*dH0_y;

dK_conjMuSub = -2*conj(u)^2*d2Sigma_pconjp + 2*d2Sigma_conjpconjp...
    -2*conj(u)^2*conj(d2Sigma_conjpconjp) + 2*conj(d2Sigma_pconjp)...
    + 2*conj(u)^2*dH0_p - 2*conj(dH0_p) + 2*conj(u)*dH0_y;

dK_vnormv0 = -2*conj(u)*d2Sigma_pconjp - 2*u*d2Sigma_conjpconjp...
    -2*conj(u)*conj(d2Sigma_conjpconjp) - 2*u*conj(d2Sigma_pconjp)...
    + 2*conj(u)*dH0_p + 2*u*conj(dH0_p) + (1-u*conj(u))*dH0_y;
dK_vnormvSub = 1/(1+u*conj(u))*dK_vnormv0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subSetLin = sym('subSetLin', [13,1]);
for j=1:13
    variable = string(symSetLin(j));
    switch variable 
        case string(u)
            subSetLin(j) = u;
        case string(lambda0)
            subSetLin(j) = lambda0Sub;
        case string(lambda1)
            subSetLin(j) = lambda1Sub;
        case string(K)
            subSetLin(j) = KSub;
        case string(dlambda0_Mu)
            subSetLin(j) = dlambda0_MuSub;
        case string(dlambda0_conjMu)
            subSetLin(j) = dlambda0_conjMuSub;
        case string(dlambda0_vnormv)
            subSetLin(j) = dlambda0_vnormvSub;
        case string(dlambda1_Mu)
            subSetLin(j) = dlambda1_MuSub;
        case string(dlambda1_conjMu)
            subSetLin(j) = dlambda1_conjMuSub;
        case string(dlambda1_vnormv)
            subSetLin(j) = dlambda1_vnormvSub;
        case string(dK_Mu)
            subSetLin(j) = dK_MuSub;
        case string(dK_conjMu)
            subSetLin(j) = dK_conjMuSub;
        case string(dK_vnormv)    
            subSetLin(j) = dK_vnormvSub;
    end
    clear variable
end

for j=1:12
    temp = WeylSLin(j,5);
    temp = subs(temp, symSetLin, subSetLin);
    temp = complex_simple3(temp,MVarPsi2);
    WeylSLin(j,5) = temp;
    clear temp
end

clearvars j
save('DataTemp_WeylCurv4_Aug26.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
load('DataTemp_WeylCurv4_Aug26.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
assumeAlso([H0, dH0_y], 'real');
countNonLinearRef = [1, 4, 8, 43, 47, 85];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WeylSN = sym('WeylSN',[6 5]);
symSetNonLinear = [];
for j=1:6
    mm = countNonLinear(j,1);
    nn = countNonLinear(j,2);
    kk = countNonLinear(j,3);
    ll = countNonLinear(j,4);
    N = countNonLinearRef(j);
    temp =  WeylSS(N,5);
    tempSet = symvar(temp);
    symSetNonLinear = union(symSetNonLinear, tempSet);
    WeylSN(j,:) = [mm, nn, kk, ll, temp];
    clearvars mm nn kk ll N temp tempSet
end

syms d3Psi_ppp d3Sigma_ppp d3Sigma_ppconjp 
syms d3Sigma_pconjpconjp d2H0_yy real
syms d2H0_pp 

d3Sigma_conjpconjpconjp = conj(d2H0_pp) + conj(d3Sigma_ppconjp);

MVarTemp = [d3Psi_ppp, d3Sigma_ppp, d3Sigma_ppconjp,...
    d3Sigma_pconjpconjp, d2H0_yy, d2H0_pp];
MVarPsi3 = union(MVarPsi2, MVarTemp);
clear MVarTemp

d2lambda1_Mu1 = -2*d2Psi_pp;
d2lambda1_Mu2 = 0;
d2lambda1_Mu3 = -2*i*d2Psi_pp;
d2lambda1_conjMu1 = 2*conj(u)^2*d2Psi_pp;
d2lambda1_conjMu2 = 0;
d2lambda1_conjMu3 = 2*i*conj(u)^2*d2Psi_pp;
d2lambda1_vnormv1 = 2*conj(u)/(1+u*conj(u))*d2Psi_pp;
d2lambda1_vnormv2 = 0;
d2lambda1_vnormv3 = 2*i*conj(u)/(1+u*conj(u))*d2Psi_pp;
%
d2lambda1_MuMuSub = [d2lambda1_Mu1,d2lambda1_Mu2,d2lambda1_Mu3]*muVec;
d2lambda1_MuconjMuSub = [d2lambda1_Mu1,d2lambda1_Mu2,d2lambda1_Mu3]*conjMuVec;
d2lambda1_MuvnormvSub = [d2lambda1_Mu1,d2lambda1_Mu2,d2lambda1_Mu3]*vnormvVec;
d2lambda1_conjMuconjMuSub=[d2lambda1_conjMu1,d2lambda1_conjMu2,d2lambda1_conjMu3]*conjMuVec;
d2lambda1_conjMuvnormvSub=[d2lambda1_conjMu1,d2lambda1_conjMu2,d2lambda1_conjMu3]*vnormvVec;
d2lambda1_vnormvvnormvSub=[d2lambda1_vnormv1,d2lambda1_vnormv2,d2lambda1_vnormv3]*vnormvVec;
%
d2lambda0_Mu1 = y*d3Psi_ppp - u*d2Psi_pp - 2*d3Sigma_ppp...
    +(2*u^2-2)*d3Sigma_ppconjp + 2*u^2*d3Sigma_pconjpconjp;
d2lambda0_Mu2 = d2Psi_pp;
d2lambda0_Mu3 = i*y*d3Psi_ppp - i*u*d2Psi_pp - 2*i*d3Sigma_ppp...
    +(2*i+2*i*u^2)*d3Sigma_ppconjp - 2*i*u^2*d3Sigma_pconjpconjp;

d2lambda0_conjMu1 = -y*conj(u)^2*d3Psi_ppp - conj(u)*d2Psi_pp + 2*conj(u)^2*d3Sigma_ppp...
    + (2*conj(u)^2-2)*d3Sigma_ppconjp - 2*d3Sigma_pconjpconjp;
d2lambda0_conjMu2 = -conj(u)^2*d2Psi_pp;
d2lambda0_conjMu3=-i*conj(u)^2*y*d3Psi_ppp-i*conj(u)*d2Psi_pp+2*i*conj(u)^2*d3Sigma_ppp...
    +(-2*i*conj(u)^2-2*i)*d3Sigma_ppconjp + 2*i*d3Sigma_pconjpconjp;

d2lambda0_vnormv1 = 1/(1+u*conj(u))*(-conj(u)*y*d3Psi_ppp + (u*conj(u)-1)/2*d2Psi_pp...
    +2*conj(u)*d3Sigma_ppp + 2*(u+conj(u))*d3Sigma_ppconjp + 2*u*d3Sigma_pconjpconjp);
d2lambda0_vnormv2 = -conj(u)*d2Psi_pp/(1+u*conj(u));
d2lambda0_vnormv3 = 1/(1+u*conj(u))*(-i*conj(u)*y*d3Psi_ppp + i*(u*conj(u)-1)/2*d2Psi_pp...
    +2*i*conj(u)*d3Sigma_ppp + 2*i*(u-conj(u))*d3Sigma_ppconjp - 2*i*u*d3Sigma_pconjpconjp);
%
d2lambda0_MuMuSub = [d2lambda0_Mu1,d2lambda0_Mu2,d2lambda0_Mu3]*muVec;
d2lambda0_MuconjMuSub = [d2lambda0_Mu1,d2lambda0_Mu2,d2lambda0_Mu3]*conjMuVec;
d2lambda0_MuvnormvSub = [d2lambda0_Mu1,d2lambda0_Mu2,d2lambda0_Mu3]*vnormvVec;
d2lambda0_conjMuconjMuSub=[d2lambda0_conjMu1,d2lambda0_conjMu2,d2lambda0_conjMu3]*conjMuVec;
d2lambda0_conjMuvnormvSub=[d2lambda0_conjMu1,d2lambda0_conjMu2,d2lambda0_conjMu3]*vnormvVec;
d2lambda0_vnormvvnormvSub=[d2lambda0_vnormv1,d2lambda0_vnormv2,d2lambda0_vnormv3]*vnormvVec;
%
d2K_Mu1 = 2*d3Sigma_ppconjp + (2-2*u^2)*d3Sigma_pconjpconjp - 2*u^2*d3Sigma_conjpconjpconjp...
    + 2*conj(d3Sigma_conjpconjpconjp) + (2-2*u^2)*conj(d3Sigma_pconjpconjp)...
    - 2*u^2*conj(d3Sigma_ppconjp) - 2*d2H0_pp + 2*u^2*conj(d2H0_pp);
d2K_Mu2 = 2*u*d2H0_yy;
d2K_Mu3 = 2*i*d3Sigma_ppconjp + (-2*i-2*i*u^2)*d3Sigma_pconjpconjp...
    +2*i*u^2*d3Sigma_conjpconjpconjp + 2*i*conj(d3Sigma_conjpconjpconjp)...
    +(-2*i-2*i*u^2)*conj(d3Sigma_pconjpconjp)+ 2*i*u^2*conj(d3Sigma_ppconjp)...
    -2*i*d2H0_pp - 2*i*u^2*conj(d2H0_pp);

d2K_conjMu1 = conj(d2K_Mu1);
d2K_conjMu2 = conj(d2K_Mu2);
d2K_conjMu3 = conj(d2K_Mu3);

d2K_vnormv1Temp = -2*conj(u)*d3Sigma_ppconjp + (-2*conj(u)-2*u)*d3Sigma_pconjpconjp...
    -2*u*d3Sigma_conjpconjpconjp - 2*conj(u)*conj(d3Sigma_conjpconjpconjp)...
    + (-2*conj(u)-2*u)*conj(d3Sigma_pconjpconjp) - 2*u*conj(d3Sigma_ppconjp)...
    + 2*conj(u)*d2H0_pp + 2*u*conj(d2H0_pp);
d2K_vnormv1 = 1/(1+u*conj(u))*d2K_vnormv1Temp;
d2K_vnormv2 = (1-u*conj(u))/(1+u*conj(u))*d2H0_yy;
d2K_vnormv3Temp =-2*i*conj(u)*d3Sigma_ppconjp + (2*i*conj(u)-2*i*u)*d3Sigma_pconjpconjp...
    +2*i*u*d3Sigma_conjpconjpconjp - 2*i*conj(u)*conj(d3Sigma_conjpconjpconjp)...
    + (2*i*conj(u)-2*i*u)*conj(d3Sigma_pconjpconjp) + 2*i*u*conj(d3Sigma_ppconjp)...
    + 2*i*conj(u)*d2H0_pp - 2*i*u*conj(d2H0_pp);
d2K_vnormv3 = 1/(1+u*conj(u))*d2K_vnormv3Temp;
%
d2K_MuMuSub = [d2K_Mu1,d2K_Mu2,d2K_Mu3]*muVec;
d2K_MuconjMuSub = [d2K_Mu1,d2K_Mu2,d2K_Mu3]*conjMuVec;
d2K_MuvnormvSub = [d2K_Mu1,d2K_Mu2,d2K_Mu3]*vnormvVec;
d2K_conjMuconjMuSub=[d2K_conjMu1,d2K_conjMu2,d2K_conjMu3]*conjMuVec;
d2K_conjMuvnormvSub=[d2K_conjMu1,d2K_conjMu2,d2K_conjMu3]*vnormvVec;
d2K_vnormvvnormvSub=[d2K_vnormv1,d2K_vnormv2,d2K_vnormv3]*vnormvVec;
clearvars d2K_vnormv1Temp d2K_vnormv3Temp

% length(symSetNonLinear) = 32
symSet1 = [u, lambda0, lambda1, K, aMu, h11];
subSet1 = [u, lambda0Sub, lambda1Sub, KSub, 0, 0];

symSet2 = [dlambda0_Mu, dlambda0_conjMu, dlambda0_vnormv,...
    dlambda1_Mu, dlambda1_conjMu, dlambda1_vnormv,...
    dK_Mu, dK_conjMu, dK_vnormv];
subSet2 = [dlambda0_MuSub, dlambda0_conjMuSub, dlambda0_vnormvSub,...
    dlambda1_MuSub, dlambda1_conjMuSub, dlambda1_vnormvSub,...
    dK_MuSub, dK_conjMuSub, dK_vnormvSub];

symSet3 = [d2lambda0_MuMu, d2lambda0_MuconjMu, d2lambda0_Muvnormv,...
    d2lambda0_conjMuconjMu, d2lambda0_conjMuvnormv, d2lambda0_vnormvvnormv];
subSet3 = [d2lambda0_MuMuSub, d2lambda0_MuconjMuSub, d2lambda0_MuvnormvSub,...
    d2lambda0_conjMuconjMuSub, d2lambda0_conjMuvnormvSub, d2lambda0_vnormvvnormvSub];

symSet4 = [d2lambda1_MuMu, d2lambda1_MuconjMu, d2lambda1_Muvnormv,...
    d2lambda1_conjMuconjMu, d2lambda1_conjMuvnormv, d2lambda1_vnormvvnormv];
subSet4 = [d2lambda1_MuMuSub, d2lambda1_MuconjMuSub, d2lambda1_MuvnormvSub,...
    d2lambda1_conjMuconjMuSub, d2lambda1_conjMuvnormvSub, d2lambda1_vnormvvnormvSub];

symSet5=[d2K_MuconjMu, d2K_Muvnormv, d2K_conjMuconjMu, d2K_conjMuvnormv, d2K_vnormvvnormv];
subSet5=[d2K_MuconjMuSub, d2K_MuvnormvSub, d2K_conjMuconjMuSub,...
    d2K_conjMuvnormvSub, d2K_vnormvvnormvSub];

for j=1:6
    temp = WeylSN(j,5);
    temp = subs(temp, symSet1, subSet1);
    temp = subs(temp, symSet2, subSet2);
    temp = subs(temp, symSet3, subSet3);
    temp = subs(temp, symSet4, subSet4);
    temp = subs(temp, symSet5, subSet5);
    temp = complex_simple3(temp, MVarPsi3);
    WeylSN(j,5) = temp;
    clear temp
end
clear j

save('DataTemp_WeylCurv4Two_Aug28.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
load('DataTemp_WeylCurv4Two_Aug28.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
assumeAlso([H0, dH0_y, d3Sigma_pconjpconjp, d2H0_yy], 'real');

WN1212 = WeylSN(1,5)/(1+u*conj(u))^2;
WN1215 = WeylSN(2,5)/(1+u*conj(u));
WN1225 = WeylSN(3,5)/(1+u*conj(u));
WN1515 = WeylSN(4,5);
WN1525 = WeylSN(5,5);
WN2525 = WeylSN(6,5);
WNvector = [WN1212, WN1215, WN1225, WN1515, WN1525, WN2525];
syms conju

for j=1:6
    temp = WNvector(j); 
    temp = subs(temp, conj(u), conju);
    WNvector(j) = temp;
    clear temp
end

indexNonLinear = {'u1212', 'u1215', 'u1225', 'u1515', 'u1525', 'u2525';
    'term1212', 'term1215', 'term1225', 'term1515', 'term1525',...
    'term2525'};

[CoeffSN.term1212, CoeffSN.u1212] = coeffs(WNvector(1), [u,conju,y]);
[CoeffSN.term1215, CoeffSN.u1215] = coeffs(WNvector(2), [u,conju,y]);
[CoeffSN.term1225, CoeffSN.u1225] = coeffs(WNvector(3), [u,conju,y]);
[CoeffSN.term1515, CoeffSN.u1515] = coeffs(WNvector(4), [u,conju,y]);
[CoeffSN.term1525, CoeffSN.u1525] = coeffs(WNvector(5), [u,conju,y]);
[CoeffSN.term2525, CoeffSN.u2525] = coeffs(WNvector(6), [u,conju,y]);

latex_WeylSN = fopen('latex_WeylCurv_Part4SN_Aug23.txt','w');

fprintf(latex_WeylSN, 'Non-linear terms in WeylSN\n');
for j=1:6
    mm = WeylSN(j,1);
    nn = WeylSN(j,2);
    kk = WeylSN(j,3);
    ll = WeylSN(j,4);
    if (j==1)
        fprintf(latex_WeylSN, 'WeylSN(%d,%d,%d,%d) is (1+u*conj(u))^2 *\n',...
            [mm,nn,kk,ll]);
    elseif (j==2)
        fprintf(latex_WeylSN, 'WeylSN(%d,%d,%d,%d) is (1+u*conj(u))*\n',...
            [mm,nn,kk,ll]);
    elseif (j==3)
        fprintf(latex_WeylSN, 'WeylSN(%d,%d,%d,%d) is (1+u*conj(u))*\n',...
            [mm,nn,kk,ll]);
    else
        fprintf(latex_WeylSN, 'WeylSN(%d,%d,%d,%d) is\n',[mm,nn,kk,ll]);
    end   
    uVecTemp = CoeffSN.(indexNonLinear{1,j});
    termVecTemp = CoeffSN.(indexNonLinear{2,j});
    N = length(uVecTemp);
    for jj=1:N
        latexU = latex(uVecTemp(jj));
        latexTerm = latex(termVecTemp(jj));
        fprintf(latex_WeylSN, '$%s $: $%s $ \n', latexU, latexTerm);
        fprintf(latex_WeylSN, '%s\n', ' ');
        clearvars latexU latexTerm
    end
    fprintf(latex_WeylSN, '%s\n', ' ');
    clearvars mm nn kk ll jj  uVecTemp termVecTemp N
end
fclose(latex_WeylSN);