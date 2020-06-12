% Update from temp_WeylCurv_Part4_Aug23_2019.m

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
countNonLinearRef = [1, 4, 8, 43, 47, 85];

v1normv = complex_simple3(v1/norm_of_v,[u]);
v2normv = complex_simple3(v2/norm_of_v,[u]);
v3normv = complex_simple3(v3/norm_of_v,[u]);
muVec = [mu1;mu2;mu3];
conjMuVec = [conj(mu1);conj(mu2);conj(mu3)];
vnormvVec = [v1normv;v2normv;v3normv];

WeylSLin = sym('WeylSlin',[12 5]);
WeylSN = sym('WeylSN',[6 5]);

symSetLin = [u, lambda0,lambda1,K,dlambda0_Mu,dlambda0_conjMu,dlambda0_vnormv...
    dlambda1_Mu,dlambda1_conjMu,dlambda1_vnormv,...
    dK_Mu, dK_conjMu, dK_vnormv]; % length = 13;

for j=1:12
    mm = countLinear(j,1);
    nn = countLinear(j,2);
    kk = countLinear(j,3);
    ll = countLinear(j,4);
    N = countLinearRef(j);
    temp = WeylSS(N,5);
    WeylSLin(j,:) = [mm,nn,kk,ll,temp];
    clearvars mm nn kk ll temp tempSet N
end

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms p
syms Psi dPsi_p d2Psi_pp
syms Zeta1 dZeta1_p d2Zeta1_pp
syms Zeta2 dZeta2_p d2Zeta2_pp
syms K0 dK0_y real

% Linear terms: WeylSLin = sym('WeylSlin',[12 5]).
MVarPsi2 = [u, p, Psi, dPsi_p, d2Psi_pp, Zeta1, dZeta1_p, d2Zeta1_pp,...
    Zeta2, dZeta2_p, d2Zeta2_pp];
lambda1Sub = Psi;
lambda0Sub = -1/2*y*dPsi_p -i/2*conj(p)*dZeta1_p -i/2*dZeta2_p;
KSub = K0 + i*(Zeta1 - conj(Zeta1));

dlambda1_x = dPsi_p;
dlambda1_y = 0;
dlambda1_z = i*dPsi_p;
dlambda0_x = -y/2*d2Psi_pp-i/2*(conj(p)*d2Zeta1_pp + dZeta1_p)...
    -i/2*d2Zeta2_pp;
dlambda0_y = -1/2*dPsi_p;
dlambda0_z = -i/2*y*d2Psi_pp + 1/2*conj(p)*d2Zeta1_pp - 1/2*dZeta1_p...
    + 1/2*d2Zeta2_pp;
dK_x = i*(dZeta1_p - conj(dZeta1_p));
dK_y = dK0_y;
dK_z = -dZeta1_p - conj(dZeta1_p);

dlambda1_MuSub = [dlambda1_x, dlambda1_y, dlambda1_z]*muVec;
dlambda1_conjMuSub = [dlambda1_x, dlambda1_y, dlambda1_z]*conjMuVec;
dlambda1_vnormvSub = [dlambda1_x, dlambda1_y, dlambda1_z]*vnormvVec;

dlambda0_MuSub = [dlambda0_x, dlambda0_y, dlambda0_z]*muVec;
dlambda0_conjMuSub = [dlambda0_x, dlambda0_y, dlambda0_z]*conjMuVec;
dlambda0_vnormvSub = [dlambda0_x, dlambda0_y, dlambda0_z]*vnormvVec;

dK_MuSub = [dK_x, dK_y, dK_z]*muVec;
dK_conjMuSub = [dK_x, dK_y, dK_z]*conjMuVec;
dK_vnormvSub = [dK_x, dK_y, dK_z]*vnormvVec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subSetLin = [u, lambda0Sub, lambda1Sub, KSub,...
    dlambda0_MuSub, dlambda0_conjMuSub, dlambda0_vnormvSub,...
    dlambda1_MuSub, dlambda1_conjMuSub, dlambda1_vnormvSub,...
    dK_MuSub, dK_conjMuSub, dK_vnormvSub];

for j=1:12
    temp = WeylSLin(j,5);
    temp = subs(temp, symSetLin, subSetLin);
    temp = complex_simple3(temp,MVarPsi2);
    WeylSLin(j,5) = temp;
    clear temp
end
clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Non-linear terms: WeylSN = sym('WeylSN',[6 5]).
syms p d3Psi_ppp d3Zeta1_ppp d3Zeta2_ppp
syms d2K0_yy real

CVar45N = [p, Psi, dPsi_p, d2Psi_pp, Zeta1, dZeta1_p, d2Zeta1_pp,...
    Zeta2, dZeta2_p, d2Zeta2_pp];
dCVar45N = sym('dCVar45',[3,10]);
dCVar45N(:,1) = [1; 0; 0]; %p
dCVar45N(:,2) = [dPsi_p; 0; 0]; %Psi
dCVar45N(:,3) = [d2Psi_pp; 0; 0]; %dPsi_p
dCVar45N(:,4) = [d3Psi_ppp; 0; 0]; %d2Psi_pp
dCVar45N(:,5) = [dZeta1_p; 0; 0]; %Zeta1
dCVar45N(:,6) = [d2Zeta1_pp; 0; 0]; %dZeta1_p
dCVar45N(:,7) = [d3Zeta1_ppp; 0; 0]; %d2Zeta1_pp
dCVar45N(:,8) = [dZeta2_p; 0; 0]; %Zeta2
dCVar45N(:,9) = [d2Zeta2_pp; 0; 0]; %dZeta2_p
dCVar45N(:,10) = [d3Zeta2_ppp; 0; 0]; %d2Zeta2_pp

RVar45N = [y, K0, dK0_y];
dRVar45N = sym('dRVar45N',[3,3]);
dRVar45N(:,1) = [0; 0; 1]; %y
dRVar45N(:,2) = [0; 0; dK0_y]; %K0
dRVar45N(:,3) = [0; 0; d2K0_yy]; %dK0_y


MVar45N = union(CVar45N, [u, d3Psi_ppp, d3Zeta1_ppp, d3Zeta2_ppp]);
% df = df_NatTwo_WeylCurv_Part45(f,CVar45N,dCVar45N,...
%     RVar45N, dRVar45N, MVar45N);
% df = [df/dp, df/dconj(p), df/dy];

dlambda1_MuRow = df_NatTwo_WeylCurv_Part45(dlambda1_MuSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2lambda1_Mux = dlambda1_MuRow(1) + dlambda1_MuRow(2);
d2lambda1_Muy = dlambda1_MuRow(3);
d2lambda1_Muz = i*dlambda1_MuRow(1) - i*dlambda1_MuRow(2);

dlambda1_conjMuRow = df_NatTwo_WeylCurv_Part45(dlambda1_conjMuSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2lambda1_conjMux = dlambda1_conjMuRow(1) + dlambda1_conjMuRow(2);
d2lambda1_conjMuy = dlambda1_conjMuRow(3);
d2lambda1_conjMuz = i*dlambda1_conjMuRow(1) - i*dlambda1_conjMuRow(2);

dlambda1_vnormvRow = df_NatTwo_WeylCurv_Part45(dlambda1_vnormvSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2lambda1_vnormvx = dlambda1_vnormvRow(1) + dlambda1_vnormvRow(2);
d2lambda1_vnormvy = dlambda1_vnormvRow(3);
d2lambda1_vnormvz = i*dlambda1_vnormvRow(1) - i*dlambda1_vnormvRow(2);

d2lambda1_MuMuSub=[d2lambda1_Mux,d2lambda1_Muy,d2lambda1_Muz]*muVec;
d2lambda1_MuconjMuSub=[d2lambda1_Mux,d2lambda1_Muy,d2lambda1_Muz]*conjMuVec;
d2lambda1_MuvnormvSub=[d2lambda1_Mux,d2lambda1_Muy,d2lambda1_Muz]*vnormvVec;
d2lambda1_conjMuconjMuSub = [d2lambda1_conjMux, d2lambda1_conjMuy,...
    d2lambda1_conjMuz]*conjMuVec;
d2lambda1_conjMuvnormvSub = [d2lambda1_conjMux, d2lambda1_conjMuy,...
    d2lambda1_conjMuz]*vnormvVec;
d2lambda1_vnormvvnormvSub = [d2lambda1_vnormvx, d2lambda1_vnormvy,...
    d2lambda1_vnormvz]*vnormvVec;
% %
dlambda0_MuRow = df_NatTwo_WeylCurv_Part45(dlambda0_MuSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2lambda0_Mux = dlambda0_MuRow(1) + dlambda0_MuRow(2);
d2lambda0_Muy = dlambda0_MuRow(3);
d2lambda0_Muz = i*dlambda0_MuRow(1) - i*dlambda0_MuRow(2);

dlambda0_conjMuRow = df_NatTwo_WeylCurv_Part45(dlambda0_conjMuSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2lambda0_conjMux = dlambda0_conjMuRow(1) + dlambda0_conjMuRow(2);
d2lambda0_conjMuy = dlambda0_conjMuRow(3);
d2lambda0_conjMuz = i*dlambda0_conjMuRow(1) - i*dlambda0_conjMuRow(2);

dlambda0_vnormvRow = df_NatTwo_WeylCurv_Part45(dlambda0_vnormvSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2lambda0_vnormvx = dlambda0_vnormvRow(1) + dlambda0_vnormvRow(2);
d2lambda0_vnormvy = dlambda0_vnormvRow(3);
d2lambda0_vnormvz = i*dlambda0_vnormvRow(1) - i*dlambda0_vnormvRow(2);

d2lambda0_MuMuSub=[d2lambda0_Mux,d2lambda0_Muy,d2lambda0_Muz]*muVec;
d2lambda0_MuconjMuSub=[d2lambda0_Mux,d2lambda0_Muy,d2lambda0_Muz]*conjMuVec;
d2lambda0_MuvnormvSub=[d2lambda0_Mux,d2lambda0_Muy,d2lambda0_Muz]*vnormvVec;
d2lambda0_conjMuconjMuSub = [d2lambda0_conjMux, d2lambda0_conjMuy,...
    d2lambda0_conjMuz]*conjMuVec;
d2lambda0_conjMuvnormvSub = [d2lambda0_conjMux, d2lambda0_conjMuy,...
    d2lambda0_conjMuz]*vnormvVec;
d2lambda0_vnormvvnormvSub = [d2lambda0_vnormvx, d2lambda0_vnormvy,...
    d2lambda0_vnormvz]*vnormvVec;
% %
dK_MuRow = df_NatTwo_WeylCurv_Part45(dK_MuSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2K_Mux = dK_MuRow(1) + dK_MuRow(2);
d2K_Muy = dK_MuRow(3);
d2K_Muz = i*dK_MuRow(1) - i*dK_MuRow(2);

dK_conjMuRow = df_NatTwo_WeylCurv_Part45(dK_conjMuSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2K_conjMux = dK_conjMuRow(1) + dK_conjMuRow(2);
d2K_conjMuy = dK_conjMuRow(3);
d2K_conjMuz = i*dK_conjMuRow(1) - i*dK_conjMuRow(2);

dK_vnormvRow = df_NatTwo_WeylCurv_Part45(dK_vnormvSub,...
    CVar45N,dCVar45N,RVar45N, dRVar45N, MVar45N);
d2K_vnormvx = dK_vnormvRow(1) + dK_vnormvRow(2);
d2K_vnormvy = dK_vnormvRow(3);
d2K_vnormvz = i*dK_vnormvRow(1) - i*dK_vnormvRow(2);

d2K_MuMuSub=[d2K_Mux,d2K_Muy,d2K_Muz]*muVec;
d2K_MuconjMuSub=[d2K_Mux,d2K_Muy,d2K_Muz]*conjMuVec;
d2K_MuvnormvSub=[d2K_Mux,d2K_Muy,d2K_Muz]*vnormvVec;
d2K_conjMuconjMuSub = [d2K_conjMux, d2K_conjMuy,d2K_conjMuz]*conjMuVec;
d2K_conjMuvnormvSub = [d2K_conjMux, d2K_conjMuy,d2K_conjMuz]*vnormvVec;
d2K_vnormvvnormvSub = [d2K_vnormvx, d2K_vnormvy,d2K_vnormvz]*vnormvVec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symSet1 = [u, lambda0, lambda1, K, aMu, h11, dlambda0_Mu, dlambda0_conjMu,...
    dlambda0_vnormv, dlambda1_Mu, dlambda1_conjMu, dlambda1_vnormv,...
    dK_Mu, dK_conjMu, dK_vnormv];
subSet1 = [u, lambda0Sub, lambda1Sub, KSub, 0, 0, dlambda0_MuSub,...
    dlambda0_conjMuSub, dlambda0_vnormvSub, dlambda1_MuSub, ...
    dlambda1_conjMuSub, dlambda1_vnormvSub, dK_MuSub, dK_conjMuSub,...
    dK_vnormvSub];

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
    temp = subs(temp, symSet3, subSet3);
    temp = subs(temp, symSet4, subSet4);
    temp = subs(temp, symSet5, subSet5);
    temp = complex_simple3(temp, MVar45N);
    WeylSN(j,5) = temp;
    clear temp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('DataTemp_WeylCurv_Part45_Aug31.mat');

%%
load('DataTemp_WeylCurv_Part45_Aug31.mat');
assumeAlso([x y z u1 u2 gamma K0 dK0_y d2K0_yy], 'real');
syms conju
MVarPsi3 = [u, p, Psi, dPsi_p, d2Psi_pp, d3Psi_ppp,...
    Zeta1, dZeta1_p, d2Zeta1_pp, d3Zeta1_ppp,...
    Zeta2, dZeta2_p, d2Zeta2_pp, d3Zeta2_ppp];
% countNonLinear = [1, 2, 1, 2;
%     1, 2, 1, 5;
%     1, 2, 2, 5;
%     1, 5, 1, 5;
%     1, 5, 2, 5;
%     2, 5, 2, 5];
WNfactor.f1212 = factor(WeylSN(1,5));
WNfactor.f1215 = factor(WeylSN(2,5));
WNfactor.f1225 = factor(WeylSN(3,5));
WNfactor.f1515 = factor(WeylSN(4,5));
WNfactor.f1525 = factor(WeylSN(5,5));
WNfactor.f2525 = factor(WeylSN(6,5));

WN1212 = WeylSN(1,5)/(1+u*conj(u))^2;
WN1215 = WeylSN(2,5)/(1+u*conj(u));
WN1225 = WeylSN(3,5)/(1+u*conj(u));
WN1515 = WeylSN(4,5);
WN1525 = WeylSN(5,5);
WN2525 = WeylSN(6,5);
WNvector = [WN1212, WN1215, WN1225, WN1515, WN1525, WN2525];
for j=1:6
    temp = WNvector(j);
    temp = complex_simple3(temp, MVarPsi3);
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

latex_WeylSN = fopen('latex_WeylCurv_Part4SN_Aug31.txt','w');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
