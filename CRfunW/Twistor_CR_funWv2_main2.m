% Twistor_CR_funWv2_main2.m
load('DataMain1_CR_funWv2.mat');
realVariables = [x,y,z];
assumeAlso(realVariables,'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \nabla_{X_m} X_n = \Gamma_{mn}^k X_k
% \nabla_{X_{\bar m}} X_n = Gamma_{\bar m n}^k X_k
% \nabla_T X_n = \Gamma_{0n}^k X_k
% Gamma.holo(m,n,k) = \Gamma_{mn}^k
% Gamma.antiholo(m,n,k) = Gamma_{\bar m n}^k 
% Gamma.T(n,k) = Gamma_{0n}^k
Gamma.holo = sym('Gamma_holo_%d_%d_%d',[2 2 2]);
Gamma.antiholo = sym('Gamma_antiholo_%d_%d_%d',[2 2 2]);
Gamma.T = sym('Gamma0_%d_%d', [2 2]);
% CVar1 = [x,y,z,u,w];
hInv = inv(h);
dh = cell(2,2);
dh{1,1} = df_main1_CR_funWv2(h(1,1),CVar1,derivativeDict);
dh{1,2} = df_main1_CR_funWv2(h(1,2),CVar1,derivativeDict);
dh{2,1} = df_main1_CR_funWv2(h(2,1),CVar1,derivativeDict);
dh{2,2} = df_main1_CR_funWv2(h(2,2),CVar1,derivativeDict);

% Define Gamma.holo(m,n,k) and Gamma.antiholo(m,n,k).
for m=1:2
    for n=1:2
        for k=1:2
            temp1 = 0; %Gamma.holo: temp1, temp2
            temp3 = 0; %Gamma.antiholo: temp3
            Xm = CRVector(:,2*m-1);
            Xn = CRVector(:,2*n-1); 
            lie_Xn_conjXm = lieMain1{2*n-1,2*m};
            for p=1:2
                lie_Xm_conjXp = lieMain1{2*m-1,2*p};
                conjXp = CRVector(:,2*p);
                temp2 = dh{n,p}*Xm - transpose(Xn)*g*lie_Xm_conjXp;
                temp1 = temp1 + hInv(p,k)*temp2;
                temp3 = temp3 - hInv(p,k)*transpose(conjXp)*g*lie_Xn_conjXm;
            end
            Gamma.holo(m,n,k) = temp1;
            Gamma.antiholo(m,n,k) = temp3;
        end
    end
end

% Define Gamma.T(n,k)
for n=1:2
    for k=1:2
        temp4 = 0;
        lie_Xn_T = lieMain1{2*n-1,5};
        for p=1:2
            conjXp = CRVector(:,2*p);
            temp4 = temp4 - hInv(p,k)*transpose(conjXp)*g*lie_Xn_T;
        end
        Gamma.T(n,k) = temp4;
    end
end

% Simplify Gamma
MVarMain2 = [u, w, dw_x, dw_y, dw_z, dw_u];
for n=1:2
    for k=1:2
        Gamma.T(n,k) = complex_simple3(Gamma.T(n,k), MVarMain2);
        for m=1:2
            Gamma.holo(m,n,k) = complex_simple3(Gamma.holo(m,n,k), MVarMain2);
            Gamma.antiholo(m,n,k) = complex_simple3(Gamma.antiholo(m,n,k), MVarMain2);
        end
    end
end

clearvars m n k p temp1 temp2 temp3 temp4 
clearvars lie_Xn_T lie_Xn_conjXm lie_Xm_conjXp
clearvars conjXp Xm Xn
save('DataMain2_CR_funWv2.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%