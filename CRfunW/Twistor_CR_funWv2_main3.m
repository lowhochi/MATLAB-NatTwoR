% Twistor_CR_funWv2_main3.m
load('DataMain2_CR_funWv2.mat');
assumeAlso(realVariables,'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d2w_xx d2w_xy d2w_xz d2w_yy d2w_yz d2w_zz
syms d2w_xu d2w_yu d2w_zu d2w_uu
CVar2 = [x, y, z, u, w, dw_x, dw_y, dw_z, dw_u];
derivativeDict.dw_x = [d2w_xx; d2w_xy; d2w_xz; d2w_xu; 0];
derivativeDict.dw_y = [d2w_xy; d2w_yy; d2w_yz; d2w_yu; 0];
derivativeDict.dw_z = [d2w_xz; d2w_yz; d2w_zz; d2w_zu; 0];
derivativeDict.dw_u = [d2w_xu; d2w_yu; d2w_zu; d2w_uu; 0];

% dGamma.holo{m,n,k} = d(Gamma.holo(m,n,k)) in [1,5]-row vector.
% dGamma.antiholo{m,n,k} = d(Gamma.antiholo(m,n,k)) in [1,5]-row vector.
% dGamma.T{n,k} = d(Gamma.T(n,k)) in [1,5]-row vector.
dGamma.holo = cell(2,2,2); 
dGamma.antiholo = cell(2,2,2);
dGamma.T = cell(2,2);
for n=1:2
    for k=1:2
        temp0 = Gamma.T(n,k);
        dGamma.T{n,k} = df_main1_CR_funWv2(temp0,CVar2,derivativeDict);
        for m=1:2
            temp1 = Gamma.holo(m,n,k);
            temp2 = Gamma.antiholo(m,n,k);
            dGamma.holo{m,n,k} = df_main1_CR_funWv2(temp1,CVar2,derivativeDict);
            dGamma.antiholo{m,n,k} = df_main1_CR_funWv2(temp2,CVar2,derivativeDict);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = sym('R', [2,2,2,2]);
% The curvature tensor R is in 2 x 2 x 2 x 2.
% R(m,n,k,ll) = R_m^n_k_oll where R(X_k, conjXll)X_m = R_m^n_k_oll \cdot X_n. 

indexChern = [1,1,1,1; 1,1,1,2; 1,1,2,1; 1,1,2,2;
    1,2,1,1; 1,2,1,2; 1,2,2,1; 1,2,2,2;
    2,1,1,1; 2,1,1,2; 2,1,2,1; 2,1,2,2;
    2,2,1,1; 2,2,1,2; 2,2,2,1; 2,2,2,2];

MVarMain3 = [];
for jj=1:16
    m = indexChern(jj,1);
    n = indexChern(jj,2);
    k = indexChern(jj,3);
    ll = indexChern(jj,4);
    temp1 = dGamma.antiholo{ll,m,n}*CRVector(:,2*k-1)...
        - dGamma.holo(k,m,n)*CRVector(:,2*ll);
    
    temp2 = 0;
    for p=1:2
        temp2=temp2-Gamma.holo(k,m,p)*Gamma.antiholo(ll,p,n)...
            + Gamma.antiholo(ll,m,p)*Gamma.holo(k,p,n)...
            + Gamma.antiholo(ll,k,p)*Gamma.holo(p,m,n)...
            - conj(Gamma.antiholo(k,ll,p))*Gamma.antiholo(p,m,n);
    end
    temp3 = 2*i*Gamma.T(m,n)*h(k,ll);
    R(m,n,k,ll) = temp1 + temp2 + temp3;
    MVarMain3 = union( MVarMain3, symvar(R(m,n,k,ll)) );
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
clearvars m n k ll jj p temp1 temp2 temp3 temp
clearvars temp0 
save('DataMain3_CR_funWv2.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%