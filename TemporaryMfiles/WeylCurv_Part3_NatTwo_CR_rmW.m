load('Data_WeylChern_Part1_May30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WeylCurv5_NatTwo_CR_rmW.m
% WeylS = Weyl when M is a flat manifold. 
% Put w(x,u) = lambda0 + lambda1*u +k*u^2 
%   - conj(lambda1)*u^3 + conj(lambda0)*u^4.
WeylS = sym('WeylS5',[6 6 6 6]);
subSetWeyl5 = sym('subSetWeyl5' ,[135,1]);
for j=1:135
    variable = string(symSetWeyl(j));
    switch variable
        case string(u)
            subSetWeyl5(j)=u;
        case string(w)
            subSetWeyl5(j)=w;
        case string(dw_Mu)
            subSetWeyl5(j)=dw_Mu;
        case string(dw_conjMu)
            subSetWeyl5(j)=dw_conjMu;
        case string(dw_vnormv)
            subSetWeyl5(j)=dw_vnormv;
        case string(dw_u)
            subSetWeyl5(j)=dw_u;
        case string(d2w_uu)
            subSetWeyl5(j)=d2w_uu;
        case string(d2w_uMu)
            subSetWeyl5(j)=d2w_uMu;
        case string(d2w_uconjMu)
            subSetWeyl5(j)=d2w_uconjMu;
        case string(d2w_uvnormv)
            subSetWeyl5(j)=d2w_uvnormv;
        case string(d2w_MuconjMu)
            subSetWeyl5(j)=d2w_MuconjMu;
        case string(d2w_Muvnormv)
            subSetWeyl5(j)=d2w_Muvnormv;
        case string(d2w_conjMuconjMu)
            subSetWeyl5(j)=d2w_conjMuconjMu;
        case string(d2w_conjMuvnormv)
            subSetWeyl5(j)=d2w_conjMuvnormv;
        case string(d2w_vnormvvnormv)
            subSetWeyl5(j)=d2w_vnormvvnormv;
        case string(d3w_uuu)
            subSetWeyl5(j)=d3w_uuu;
        case string(d3w_uuMu)
            subSetWeyl5(j)=d3w_uuMu;
        case string(d3w_uuconjMu)
            subSetWeyl5(j)=d3w_uuconjMu;
        case string(d3w_uuvnormv)
            subSetWeyl5(j)=d3w_uuvnormv;
        case string(d3w_uMuconjMu)
            subSetWeyl5(j)=d3w_uMuconjMu;
        case string(d3w_uMuvnormv)
            subSetWeyl5(j)=d3w_uMuvnormv;
        case string(d3w_uconjMuconjMu)
            subSetWeyl5(j)=d3w_uconjMuconjMu;
        case string(d3w_uconjMuvnormv)
            subSetWeyl5(j)=d3w_uconjMuvnormv;           
        otherwise 
            subSetWeyl5(j) = 0;
    end
end
for j=1:1296
    mm = countSet(j,1);
    nn = countSet(j,2);
    pp = countSet(j,3);
    qq = countSet(j,4);
    temp = Weyl(mm,nn,pp,qq);
    temp = subs(temp,symSetWeyl,subSetWeyl5);
    WeylS(mm,nn,pp,qq) = temp;
end
MVarS = setdiff(subSetWeyl5, [0]);
clearvars mm nn pp qq j temp subSetWeyl5 variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms lambda0 lambda1 
syms dlambda0_Mu dlambda0_conjMu dlambda0_vnormv
syms d2lambda0_MuMu d2lambda0_MuconjMu d2lambda0_Muvnormv
syms d2lambda0_conjMuconjMu d2lambda0_conjMuvnormv d2lambda0_vnormvvnormv
syms dlambda1_Mu dlambda1_conjMu dlambda1_vnormv
syms d2lambda1_MuMu d2lambda1_MuconjMu d2lambda1_Muvnormv
syms d2lambda1_conjMuconjMu d2lambda1_conjMuvnormv d2lambda1_vnormvvnormv
syms K % K=K(x) real
syms dK_Mu dK_conjMu dK_vnormv
syms d2K_MuMu d2K_MuconjMu d2K_Muvnormv
syms d2K_conjMuconjMu d2K_conjMuvnormv d2K_vnormvvnormv

dlambda0Row = [dlambda0_Mu, dlambda0_conjMu, dlambda0_vnormv];
d2lambda0_conjMuMu = d2lambda0_MuconjMu + dlambda0Row*lieMu{1,2};
d2lambda0_vnormvMu = d2lambda0_Muvnormv + dlambda0Row*lieMu{1,3};
d2lambda0_vnormvconjMu = d2lambda0_conjMuvnormv + dlambda0Row*lieMu{2,3};
d2lambda0_Muu = dlambda0Row*lieMu{4,1};
d2lambda0_conjMuconju = dlambda0Row*lieMu{4,2};
d2lambda0_vnormvu = dlambda0Row*lieMu{4,3};
d2lambda0_vnormvconju = dlambda0Row*lieMu{5,3};

dlambda1Row = [dlambda1_Mu, dlambda1_conjMu, dlambda1_vnormv];
d2lambda1_conjMuMu = d2lambda1_MuconjMu + dlambda1Row*lieMu{1,2};
d2lambda1_vnormvMu = d2lambda1_Muvnormv + dlambda1Row*lieMu{1,3};
d2lambda1_vnormvconjMu = d2lambda1_conjMuvnormv + dlambda1Row*lieMu{2,3};
d2lambda1_Muu = dlambda1Row*lieMu{4,1};
d2lambda1_conjMuconju = dlambda1Row*lieMu{4,2};
d2lambda1_vnormvu = dlambda1Row*lieMu{4,3};
d2lambda1_vnormvconju = dlambda1Row*lieMu{5,3};

dKRow = [dK_Mu, dK_conjMu, dK_vnormv];
d2K_conjMuMu = d2K_MuconjMu + dKRow*lieMu{1,2};
d2K_vnormvMu = d2K_Muvnormv + dKRow*lieMu{1,3};
d2K_vnormvconjMu = d2K_conjMuvnormv + dKRow*lieMu{2,3};
d2K_Muu = dKRow*lieMu{4,1};
d2K_conjMuconju = dKRow*lieMu{4,2};
d2K_vnormvu = dKRow*lieMu{4,3};
d2K_vnormvconju = dKRow*lieMu{5,3};

wS5 = lambda0 + lambda1*u + K*u^2 - conj(lambda1)*u^3 + conj(lambda0)*u^4;
dwS5_Mu = dlambda0_Mu + dlambda1_Mu*u + dK_Mu*u^2 ...
    - conj(dlambda1_conjMu)*u^3 + conj(dlambda0_conjMu)*u^4;
dwS5_conjMu = dlambda0_conjMu + dlambda1_conjMu*u + dK_conjMu*u^2 ...
    - conj(dlambda1_Mu)*u^3 + conj(dlambda0_Mu)*u^4;
dwS5_vnormv = dlambda0_vnormv + dlambda1_vnormv*u + dK_vnormv*u^2 ...
    - conj(dlambda1_vnormv)*u^3 + conj(dlambda0_vnormv)*u^4;
dwS5_u = lambda1 + 2*K*u - 3*conj(lambda1)*u^2 + 4*conj(lambda0)*u^3;
%
d2wS5_uu = 2*K - 6*conj(lambda1)*u + 12*conj(lambda0)*u^2;
d2wS5_uMu = dlambda1_Mu + 2*dK_Mu*u - 3*conj(dlambda1_conjMu)*u^2 ...
    + 4*conj(dlambda0_conjMu)*u^3;
d2wS5_uconjMu = dlambda1_conjMu + 2*dK_conjMu*u - 3*conj(dlambda1_Mu)*u^2 ...
    + 4*conj(dlambda0_Mu)*u^3;
d2wS5_uvnormv = dlambda1_vnormv+2*dK_vnormv*u-3*conj(dlambda1_vnormv)*u^2 ...
    + 4*conj(dlambda0_vnormv)*u^3;
d2wS5_MuconjMu = d2lambda0_MuconjMu + d2lambda1_MuconjMu*u ...
    + d2K_MuconjMu*u^2 - conj(d2lambda1_conjMuMu)*u^3 ...
    + conj(d2lambda0_conjMuMu)*u^4;
d2wS5_Muvnormv = d2lambda0_Muvnormv + d2lambda1_Muvnormv*u ...
    + d2K_Muvnormv*u^2 - conj(d2lambda1_conjMuvnormv)*u^3 ...
    + conj(d2lambda0_conjMuvnormv)*u^4;
d2wS5_conjMuconjMu = d2lambda0_conjMuconjMu + d2lambda1_conjMuconjMu*u ...
    + d2K_conjMuconjMu*u^2 - conj(d2lambda1_MuMu)*u^3 ...
    + conj(d2lambda0_MuMu)*u^4;
d2wS5_conjMuvnormv = d2lambda0_conjMuvnormv + d2lambda1_conjMuvnormv*u ...
    + d2K_conjMuvnormv*u^2 - conj(d2lambda1_Muvnormv)*u^3 ...
    + conj(d2lambda0_Muvnormv)*u^4;
d2wS5_vnormvvnormv = d2lambda0_vnormvvnormv + d2lambda1_vnormvvnormv*u ...
    + d2K_vnormvvnormv*u^2 - conj(d2lambda1_vnormvvnormv)*u^3 ...
    + conj(d2lambda0_vnormvvnormv)*u^4;
%
d3wS5_uuu = -6*conj(lambda1) + 24*conj(lambda0)*u;
d3wS5_uuMu = 2*dK_Mu-6*conj(dlambda1_conjMu)*u+12*conj(dlambda0_conjMu)*u^2;
d3wS5_uuconjMu = 2*dK_conjMu-6*conj(dlambda1_Mu)*u+12*conj(dlambda0_Mu)*u^2;
d3wS5_uuvnormv = 2*dK_vnormv-6*conj(dlambda1_vnormv)*u ...
    +12*conj(dlambda0_vnormv)*u^2;
d3wS5_uMuconjMu = d2lambda1_MuconjMu + 2*d2K_MuconjMu*u ...
    - 3*conj(d2lambda1_conjMuMu)*u^2 + 4*conj(d2lambda0_conjMuMu)*u^3;
d3wS5_uMuvnormv = d2lambda1_Muvnormv + 2*d2K_Muvnormv*u ...
    -3*conj(d2lambda1_conjMuvnormv)*u^2 + 4*conj(d2lambda0_conjMuvnormv)*u^3;
d3wS5_uconjMuconjMu = d2lambda1_conjMuconjMu + 2*d2K_conjMuconjMu*u ...
    - 3*conj(d2lambda1_MuMu)*u^2 + 4*conj(d2lambda0_MuMu)*u^3;
d3wS5_uconjMuvnormv = d2lambda1_conjMuvnormv + 2*d2K_conjMuvnormv*u ...
    - 3*conj(d2lambda1_Muvnormv)*u^2 + 4*conj(d2lambda0_Muvnormv)*u^3; 

subSetwS5 = sym('subsSetwS5',[23,1]);
% length(MVarS) = 23
for j=1:23
    variable = string(MVarS(j));
    switch variable
        case string(u)
            subSetwS5(j)=u;
        case string(w)
            subSetwS5(j) = wS5;
        case string(dw_Mu)
            subSetwS5(j) = dwS5_Mu; 
        case string(dw_conjMu)
            subSetwS5(j) = dwS5_conjMu;
        case string(dw_vnormv)
            subSetwS5(j) = dwS5_vnormv;
        case string(dw_u)
            subSetwS5(j) = dwS5_u;
        case string(d2w_uu)
            subSetwS5(j) = d2wS5_uu;
        case string(d2w_uMu)
            subSetwS5(j) = d2wS5_uMu;
        case string(d2w_uconjMu)
            subSetwS5(j) = d2wS5_uconjMu;
        case string(d2w_uvnormv)
            subSetwS5(j) = d2wS5_uvnormv;
        case string(d2w_MuconjMu)
            subSetwS5(j) = d2wS5_MuconjMu;
        case string(d2w_Muvnormv)
            subSetwS5(j) = d2wS5_Muvnormv;
        case string(d2w_conjMuconjMu)
            subSetwS5(j) = d2wS5_conjMuconjMu;
        case string(d2w_conjMuvnormv)
            subSetwS5(j) = d2wS5_conjMuvnormv;
        case string(d2w_vnormvvnormv)
            subSetwS5(j) = d2wS5_vnormvvnormv;
        case string(d3w_uuu)
            subSetwS5(j) = d3wS5_uuu;
        case string(d3w_uuMu)
            subSetwS5(j) = d3wS5_uuMu;
        case string(d3w_uuconjMu)
            subSetwS5(j) = d3wS5_uuconjMu;
        case string(d3w_uuvnormv)
            subSetwS5(j) = d3wS5_uuvnormv;
        case string(d3w_uMuconjMu)
            subSetwS5(j) = d3wS5_uMuconjMu;
        case string(d3w_uMuvnormv)
            subSetwS5(j) = d3wS5_uMuvnormv;
        case string(d3w_uconjMuconjMu)
            subSetwS5(j) = d3wS5_uconjMuconjMu;
        case string(d3w_uconjMuvnormv)        
            subSetwS5(j) = d3wS5_uconjMuvnormv;
        otherwise % This part is added later.
            subSetwS5(j)=0;
    end
end

for j=1:1296
    mm = countSet(j,1);
    nn = countSet(j,2);
    pp = countSet(j,3);
    qq = countSet(j,4);
    tempS5 = WeylS(mm,nn,pp,qq);
    tempS5 = subs(tempS5, MVarS, subSetwS5); %Key!
    WeylS(mm,nn,pp,qq) = tempS5;
end
clearvars mm nn pp qq j tempS5 subSetwS5 variable
clearvars wS5 dwS5_Mu dwS5_conjMu dwS5_vnormv dwS5_u
clearvars d2wS5_MuconjMu d2wS5_Muvnormv d2wS5_conjMuconjMu
clearvars d2wS5_conjMuvnormv d2wS5_uMu d2wS5_uconjMu 
clearvars d2wS5_uu d2wS5_uvnormv d2wS5_vnormvvnormv 
clearvars d3wS5_uMuconjMu d3wS5_uMuvnormv d3wS5_uconjMuconjMu
clearvars d3wS5_uconjMuvnormv d3wS5_uuMu d3wS5_uuconjMu d3wS5_uuu 
clearvars d3wS5_uuvnormv
clearvars d2aMu d2aV d2h11 d2rho d2T4 d2w d2wRow d3T4_u d3w_u
clearvars d2T4Row daMuRow daVRow dh11Row dKRow drhoRow dlambda1Row dlambda0Row
clearvars dwRow dT4Row
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Data_WeylChern_Part3_May30.mat');