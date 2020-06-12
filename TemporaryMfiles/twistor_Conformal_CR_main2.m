% twistor_Conformal_CR_main2.m
load('Data_Conformal_CR_Main1.mat');
% variableSet = union(symvar(GammaC.holo), symvar(GammaC.antiholo));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVar2 = [u, w, f, dw_u, dw_Mu, dw_conjMu, df_Mu, df_u];
derivativeDict.dw_u = [d2w_uMu; d2w_uconjMu; d2w_uvnormv; d2w_uu; 0];
derivativeDict.dw_Mu = [d2w_MuMu; d2w_MuconjMu; d2w_Muvnormv; d2w_Muu; 0];
derivativeDict.dw_conjMu = [d2w_conjMuMu; d2w_conjMuconjMu; 
    d2w_conjMuvnormv; d2w_conjMuu; d2w_conjMuconju];
% dGamma.holo{m,n,k} = d(\Gamma_{mn}^k)
% dGamma.antiholo{m,n,k} = d(\Gamma_{\overline{m}n}^k)
% dGamma.T{n,k} = d(\Gamma_{0n}^k) 
% in (Mu, conjMu, vnormv, d/du, d/dconju). 
dGammaC.holo = cell(2,2,2);
dGammaC.antiholo = cell(2,2,2);
for n=1:2
    for k=1:2
        for m=1:2
            temp1 = GammaC.holo(m,n,k);
            temp2 = GammaC.antiholo(m,n,k);
            dGammaC.holo{m,n,k}= df_Conformal_main1(temp1,CVar2,derivativeDict);
            dGammaC.antiholo{m,n,k}= df_Conformal_main1(temp2,CVar2,derivativeDict);
        end
    end
end
clearvars n k m temp1 temp2
% The curvature tensor R is in 2 x 2 x 2 x 2.
% R(m,n,k,l) = R_m^n_k_ol where R(X_k, conj_X_l)X_m = R_m^n_k_ol X_n. 
Rc = sym('Rc_%d_%d_%d_%d', [2,2,2,2]);
for m=1:2
    for n=1:2
        for k=1:2
            for ll=1:2
                tempPart1 = dGammaC.antiholo{ll,m,n}*Uvector(:,2*k-1)...
                    - dGammaC.holo{k,m,n}*Uvector(:,2*ll);
                tempPart2 = 0;
                for p=1:2
                    tempPart2=tempPart2-GammaC.holo(k,m,p)*GammaC.antiholo(ll,p,n)...
                        + GammaC.antiholo(ll,m,p)*GammaC.holo(k,p,n)...
                        + GammaC.antiholo(ll,k,p)*GammaC.holo(p,m,n)...
                        - conj(GammaC.antiholo(k,ll,p))*GammaC.antiholo(p,m,n);
                end
                tempPart3 = 2*i*GammaC.T(m,n)*hC(k,ll);
                Rc(m,n,k,ll) = tempPart1 + tempPart2 + tempPart3;
            end
        end
    end
end
MVarRc = symvar(Rc);
for m=1:2
    for n=1:2
        for k=1:2
            for ll=1:2
                temp = Rc(m,n,k,ll);
                temp = subs(temp, conjRealVariable, realVariable);
                Rc(m,n,k,ll) = complex_simple3(temp, MVarRc); 
            end
        end
    end
end
clearvars m n k ll p temp tempPart1 tempPart2 tempPart3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ric(m,n) = ric(X_m, conj(X_n)) = R_{m\bar n}.
ricC = sym('ricC_%d_%d', [2,2]);
rhoC = 0;
for m=1:2
    for n=1:2
        ricC(m,n) = Rc(m,1,1,n)+ Rc(m,2,2,n);
        rhoC = rhoC + hCinv(n,m)*ricC(m,n);
    end
end
% rho = scalar curvature of the Tanaka Webster connection.
rhoC = complex_simple3(rhoC, MVarRc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars m n 
save('Data_Conformal_CR_Main2.mat');