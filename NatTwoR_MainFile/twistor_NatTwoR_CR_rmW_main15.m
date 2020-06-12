% twistor_NatTwoR_CR_rmW_main15.m
% b/f twistor_NatTwoR_CR_rmW_main1.m
% syms u w
% syms h11 T4 aM aV bV
% syms dw_Mu dw_conjMu dw_vnormv dw_u
% syms dT4_Mu dT4_conjMu dT4_vnormv dT4_u dT4_conju
% syms dh11_Mu dh11_conjMu dh11_vnormv dh11_u dh11_conju

syms daM_Mu daM_conjMu daM_vnormv daM_u daM_conju % Chang aMu to aM
syms daV_Mu daV_conjMu daV_vnormv daV_u daV_conju
syms dbV_Mu dbV_conjMu dbV_vnormv dbV_u dbV_conju

load('DataMain1_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lieMu: Lie brackets between vectors mu, conjMu, vnormv, u and conju;
% lieMu is used to construct 2nd derivatives of variables aM, aV, bV, etc;
% lie_Mu_conjMu = aM*mu - conj(aM)*conjMu + 2*i*h11*vnormv;
lieMu = cell(5,5);
lie_Mu_conjMu = [aM; -conj(aM); 2*i*h11];
lie_Mu_vnormv = [aV; bV; 2*T4];
lie_Mu_ddu = [-2*conj(u)/(1+u*conj(u)); 0; -2];
lie_conjMu_vnormv = [conj(bV); conj(aV); 2*conj(T4)];
lie_conjMu_ddconju = [0; -2*u/(1+u*conj(u)); -2];
lie_vnormv_ddu = [0; 1/(1+u*conj(u))^2; 0];
lie_vnormv_ddconju = [1/(1+u*conj(u))^2; 0; 0];

lieMu{1,1} = zeros(3,1);
lieMu{1,2} = lie_Mu_conjMu;
lieMu{1,3} = lie_Mu_vnormv;
lieMu{1,4} = lie_Mu_ddu;
lieMu{1,5} = zeros(3,1);
lieMu{2,1} = -lie_Mu_conjMu;
lieMu{2,2} = zeros(3,1);
lieMu{2,3} = lie_conjMu_vnormv;
lieMu{2,4} = zeros(3,1);
lieMu{2,5} = lie_conjMu_ddconju;
lieMu{3,1} = -lie_Mu_vnormv;
lieMu{3,2} = -lie_conjMu_vnormv;
lieMu{3,3} = zeros(3,1);
lieMu{3,4} = lie_vnormv_ddu;
lieMu{3,5} = lie_vnormv_ddconju;
lieMu{4,1} = -lie_Mu_ddu;
lieMu{4,2} = zeros(3,1);
lieMu{4,3} = -lie_vnormv_ddu;
lieMu{4,4} = zeros(3,1);
lieMu{4,5} = zeros(3,1);
lieMu{5,1} = zeros(3,1);
lieMu{5,2} = -lie_conjMu_ddconju;
lieMu{5,3} = -lie_vnormv_ddconju;
lieMu{5,4} = zeros(3,1);
lieMu{5,5} = zeros(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d2w_uu d2w_uMu d2w_uconjMu d2w_uvnormv
syms d2w_MuMu d2w_MuconjMu d2w_Muvnormv d2w_conjMuconjMu 
syms d2w_conjMuvnormv d2w_vnormvvnormv
dwRow = [dw_Mu, dw_conjMu, dw_vnormv];
d2w_conjMuMu = d2w_MuconjMu + dwRow*lieMu{1,2};
d2w_vnormvMu = d2w_Muvnormv + dwRow*lieMu{1,3}; 
d2w_vnormvconjMu = d2w_conjMuvnormv + dwRow*lieMu{2,3};
d2w_Muu = d2w_uMu + dwRow*lieMu{4,1}; 
d2w_conjMuu = d2w_uconjMu + dwRow*lieMu{4,2};
d2w_vnormvu = d2w_uvnormv + dwRow*lieMu{4,3};
d2w_conjMuconju = dwRow*lieMu{5,2};
d2w_vnormvconju = dwRow*lieMu{5,3};
%
d2w = [d2w_MuMu, d2w_MuconjMu, d2w_Muvnormv, d2w_Muu, 0;
    d2w_conjMuMu, d2w_conjMuconjMu, d2w_conjMuvnormv, d2w_conjMuu, d2w_conjMuconju;
    d2w_vnormvMu, d2w_vnormvconjMu, d2w_vnormvvnormv, d2w_vnormvu, d2w_vnormvconju;
    d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_uu, 0;
    0, 0, 0, 0, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d2h11_conjuconju d2h11_uconju d2h11_uu
syms d2h11_MuMu d2h11_MuconjMu d2h11_Muvnormv d2h11_conjMuconjMu 
syms d2h11_conjMuvnormv d2h11_vnormvvnormv
syms d2h11_uMu d2h11_uconjMu d2h11_uvnormv 
syms d2h11_conjuMu d2h11_conjuconjMu d2h11_conjuvnormv
dh11Row = [dh11_Mu, dh11_conjMu, dh11_vnormv];
d2h11_conjMuMu = d2h11_MuconjMu + dh11Row*lieMu{1,2};
d2h11_vnormvMu = d2h11_Muvnormv + dh11Row*lieMu{1,3};
d2h11_vnormvconjMu = d2h11_conjMuvnormv + dh11Row*lieMu{2,3};
d2h11_Muu = d2h11_uMu + dh11Row*lieMu{4,1};
d2h11_conjMuu = d2h11_uconjMu + dh11Row*lieMu{4,2};
d2h11_vnormvu = d2h11_uvnormv + dh11Row*lieMu{4,3};
d2h11_Muconju = d2h11_conjuMu + dh11Row*lieMu{5,1};
d2h11_conjMuconju = d2h11_conjuconjMu + dh11Row*lieMu{5,2};
d2h11_vnormvconju = d2h11_conjuvnormv + dh11Row*lieMu{5,3}; 

d2h11 = [d2h11_MuMu, d2h11_MuconjMu, d2h11_Muvnormv, d2h11_Muu, d2h11_Muconju;
    d2h11_conjMuMu, d2h11_conjMuconjMu, d2h11_conjMuvnormv, ...
        d2h11_conjMuu, d2h11_conjMuconju;
    d2h11_vnormvMu, d2h11_vnormvconjMu, d2h11_vnormvvnormv, ...
        d2h11_vnormvu, d2h11_vnormvconju;
    d2h11_uMu, d2h11_uconjMu, d2h11_uvnormv, d2h11_uu, d2h11_uconju;
    d2h11_conjuMu, d2h11_conjuconjMu, d2h11_conjuvnormv, d2h11_uconju, d2h11_conjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d2T4_uMu d2T4_uconjMu d2T4_uvnormv d2T4_uu d2T4_uconju d2T4_conjuconju
syms d2T4_MuMu d2T4_MuconjMu d2T4_Muvnormv d2T4_conjMuconjMu
syms d2T4_conjMuvnormv d2T4_vnormvvnormv
syms d2T4_conjuMu d2T4_conjuconjMu d2T4_conjuvnormv 
dT4Row = [dT4_Mu, dT4_conjMu, dT4_vnormv];
d2T4_conjMuMu = d2T4_MuconjMu + dT4Row*lieMu{1,2};
d2T4_vnormvMu = d2T4_Muvnormv + dT4Row*lieMu{1,3};
d2T4_vnormvconjMu = d2T4_conjMuvnormv + dT4Row*lieMu{2,3};
d2T4_Muu = d2T4_uMu + dT4Row*lieMu{4,1};
d2T4_conjMuu = d2T4_uconjMu + dT4Row*lieMu{4,2};
d2T4_vnormvu = d2T4_uvnormv + dT4Row*lieMu{4,3};
d2T4_Muconju = d2T4_conjuMu + dT4Row*lieMu{5,1};
d2T4_conjMuconju = d2T4_conjuconjMu + dT4Row*lieMu{5,2};
d2T4_vnormvconju = d2T4_conjuvnormv + dT4Row*lieMu{5,3}; 

d2T4 = [d2T4_MuMu, d2T4_MuconjMu, d2T4_Muvnormv, d2T4_Muu, d2T4_Muconju;
    d2T4_conjMuMu, d2T4_conjMuconjMu, d2T4_conjMuvnormv, ...
        d2T4_conjMuu, d2T4_conjMuconju;
    d2T4_vnormvMu, d2T4_vnormvconjMu, d2T4_vnormvvnormv, ...
        d2T4_vnormvu, d2T4_vnormvconju;
    d2T4_uMu, d2T4_uconjMu, d2T4_uvnormv, d2T4_uu, d2T4_uconju;
    d2T4_conjuMu, d2T4_conjuconjMu, d2T4_conjuvnormv, d2T4_uconju, d2T4_conjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVarMain15 = [u, w, h11, T4, aM, aV, bV,...
    dw_u, dw_Mu, dw_conjMu, dw_vnormv,...
    dh11_u, dh11_conju, dh11_Mu, dh11_conjMu, dh11_vnormv,...
    dT4_u, dT4_conju, dT4_Mu, dT4_conjMu, dT4_vnormv];
dCVarMain15 = sym('dCVarMain15', [5,21]);

dCVarMain15(:,1) = [0; 0; 0; 1; 0]; %u
dCVarMain15(:,2) = [dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0]; %w
dCVarMain15(:,3) = [dh11_Mu; dh11_conjMu; dh11_vnormv; dh11_u; dh11_conju]; %h11
dCVarMain15(:,4) = [dT4_Mu; dT4_conjMu; dT4_vnormv; dT4_u; dT4_conju]; %T4
dCVarMain15(:,5) = [daM_Mu; daM_conjMu; daM_vnormv; daM_u; daM_conju]; %aM
dCVarMain15(:,6) = [daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVarMain15(:,7) = [dbV_Mu; dbV_conjMu; dbV_vnormv; dbV_u; dbV_conju]; %bV
dCVarMain15(:,8) = transpose(d2w(4,:)); %dw_u
dCVarMain15(:,9) = transpose(d2w(1,:)); %dw_Mu
dCVarMain15(:,10) = transpose(d2w(2,:)); %dw_conjMu
dCVarMain15(:,11) = transpose(d2w(3,:)); %dw_vnormv
dCVarMain15(:,12) = transpose(d2h11(4,:)); %dh11_u
dCVarMain15(:,13) = transpose(d2h11(5,:)); %dh11_conju
dCVarMain15(:,14) = transpose(d2h11(1,:)); %dh11_Mu
dCVarMain15(:,15) = transpose(d2h11(2,:)); %dh11_conjMu
dCVarMain15(:,16) = transpose(d2h11(3,:)); %dh11_vnormv
dCVarMain15(:,17) = transpose(d2T4(4,:)); %dT4_u
dCVarMain15(:,18) = transpose(d2T4(5,:)); %dT4_conju
dCVarMain15(:,19) = transpose(d2T4(1,:)); %dT4_Mu
dCVarMain15(:,20) = transpose(d2T4(2,:)); %dT4_conjMu
dCVarMain15(:,21) = transpose(d2T4(3,:)); %dT4_vnormv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% df = df_NatTwo_MuSet_CR_rmW(f,CVarMain15,dCVarMain15)
% dGamma.holo{m,n,k} = d(\Gamma_{mn}^k) in 
% dGamma.antiholo{m,n,k} = d(\Gamma_{\overline{m}n}^k)
% dGamma.T{n,k} = d(\Gamma_{0n}^k) 
% in (Mu, conjMu, vnormv, d/du, d/dconju). 
dGamma.holo = cell(2,2,2);
dGamma.antiholo = cell(2,2,2);
dGamma.T = cell(2,2);
for n=1:2
    for k=1:2
        for m=1:2
            temp1 = Gamma.holo(m,n,k);
            temp2 = Gamma.antiholo(m,n,k);
            dGamma.holo{m,n,k}=df_NatTwo_MuSet_CR_rmW(temp1,CVarMain15,dCVarMain15);
            dGamma.antiholo{m,n,k}= df_NatTwo_MuSet_CR_rmW(temp2,CVarMain15,dCVarMain15);
        end
        temp3 = Gamma.T(n,k);
        dGamma.T{n,k}=df_NatTwo_MuSet_CR_rmW(temp3,CVarMain15,dCVarMain15);
    end
end
% The curvature tensor R is in 2 x 2 x 2 x 2.
% R(m,n,k,l) = R_m^n_k_ol where R(X_k, conj_X_l)X_m = R_m^n_k_ol X_n. 
R = sym('R_%d_%d_%d_%d', [2,2,2,2]);
for m=1:2
    for n=1:2
        for k=1:2
            for ll=1:2
                tempPart1 = dGamma.antiholo{ll,m,n}*CRvector(:,2*k-1)...
                    - dGamma.holo{k,m,n}*CRvector(:,2*ll);
                tempPart2 = 0;
                for p=1:2
                    tempPart2=tempPart2-Gamma.holo(k,m,p)*Gamma.antiholo(ll,p,n)...
                        + Gamma.antiholo(ll,m,p)*Gamma.holo(k,p,n)...
                        + Gamma.antiholo(ll,k,p)*Gamma.holo(p,m,n)...
                        - conj(Gamma.antiholo(k,ll,p))*Gamma.antiholo(p,m,n);
                end
                tempPart3 = 2*i*Gamma.T(m,n)*h(k,ll);
                R(m,n,k,ll) = tempPart1 + tempPart2 + tempPart3;
            end
        end
    end
end

% ric(m,n) = ric(X_m, conj(X_n)) = R_{m\bar n}.
ric = sym('ric_%d_%d', [2,2]);
rho = 0;
for m=1:2
    for n=1:2
        temp = 0;
        for k=1:2
            temp = temp + R(m,k,k,n);
        end        
        ric(m,n) = temp;
    end
end
% rho = scalar curvature of the Tanaka Webster connection.
for m = 1:2
    for n = 1:2
        rho = rho + hInv(n,m)*ric(m,n);
    end
end
ricci_curv = {ric, rho};
rho = complex_simple3(rho,symvar(rho));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars d2w d2T4 d2h11 dwRow dh11Row dT4Row m n k ll
clearvars temp temp1 temp2 temp3 tempPart1 tempPart2 tempPart3
list_of_variables_Main15 = who;

save('DataMain15_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%