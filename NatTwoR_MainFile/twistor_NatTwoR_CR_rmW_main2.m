% twistor_NatTwoR_CR_rmW_main2.m
% syms x y z u1 u2 gamma real % u=u1+i*u2
syms gamma real % gamma is real by default
syms u w
syms h11 T4 aM aV bV 
syms dw_Mu dw_conjMu dw_vnormv dw_u 
syms dT4_Mu dT4_conjMu dT4_vnormv dT4_u dT4_conju
syms dh11_Mu dh11_conjMu dh11_vnormv dh11_u dh11_conju
syms daM_Mu daM_conjMu daM_vnormv daM_u daM_conju % change aMu to aM
syms daV_Mu daV_conjMu daV_vnormv daV_u daV_conju
syms dbV_Mu dbV_conjMu dbV_vnormv dbV_u dbV_conju
syms d2h11_uu d2h11_uconju d2h11_conjuconju
syms rho % scalar curvature
syms drho_Mu drho_conjMu drho_vnormv drho_u drho_conju
syms d2w_uu d2w_uMu d2w_uconjMu d2w_uvnormv
syms d2T4_uMu d2T4_uconjMu d2T4_uvnormv d2T4_uu d2T4_uconju
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1 =u^2-1; mu2= 2*u; mu3= i*(u^2+1);
v1 = i*(mu2*conj(mu3)-mu3*conj(mu2));
v2 = i*(mu3*conj(mu1)-mu1*conj(mu3));
v3 = i*(mu1*conj(mu2)-mu2*conj(mu1));
v1 = complex_simple3(v1,[u]);
v2 = complex_simple3(v2,[u]);
v3 = complex_simple3(v3,[u]);
norm_of_v = sqrt(v1*v1+v2*v2+v3*v3);
norm_of_v = complex_simple3(norm_of_v,[u]);
v1normv = v1/norm_of_v; % change T1 to v1normv etc
v2normv = v2/norm_of_v;
v3normv = v3/norm_of_v;
v1normv = complex_simple3(v1normv, [u,w]);
v2normv = complex_simple3(v2normv, [u,w]);
v3normv = complex_simple3(v3normv, [u,w]);

conj_X1 = [u^2-1; 2*u; i*(u^2+1); w; 0];
X1 = [conj(u)^2-1; 2*conj(u); -i*(conj(u)^2+1); 0; conj(w)];
X2 = [0;0;0;1;0];
conj_X2 = [0;0;0;0;1];
T = [v1normv; v2normv; v3normv; T4; conj(T4)];

h = [h11, -i;
    i, 0];
hInv = inv(h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fefferman metric F0 in [X1, conj_X1, X2, conj_X2, T, d/d\gamma];
% Change on Apr 20 2019.
F0 = sym('F0',[6 6]);

F_T_T = -1/24*rho + i/2*(-conj(aV)+2*u*conj(T4)/(1+u*conj(u))-dT4_u);
F_X1_T =i/4*(conj(dw_u)- 2*u/(1+u*conj(u))*conj(w)-2*conj(T4)+aM);
F_X2_T = -i*conj(u)/(2*(1+u*conj(u)));

F0 = [0, h11, 0, -i, F_X1_T, 0;
    h11, 0, i, 0, conj(F_X1_T), 0;
    0, i, 0, 0, F_X2_T, 0;
    -i, 0, 0, 0, conj(F_X2_T), 0;
    F_X1_T, conj(F_X1_T), F_X2_T, conj(F_X2_T), F_T_T, 1/4;
    0, 0, 0, 0, 1/4, 0];
F0Inv = myInverse(F0); % myInverse.m and so on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVarMain2= [u, w, h11, T4, aM, aV, rho, dw_u, dT4_u]; % = symvar(F0)
dCVarMain2= sym('dCVarMain2', [5,9]);
dCVarMain2(:,1) = [0; 0; 0; 1; 0];
dCVarMain2(:,2) = [dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0];
dCVarMain2(:,3) = [dh11_Mu; dh11_conjMu; dh11_vnormv; dh11_u; dh11_conju];
dCVarMain2(:,4) = [dT4_Mu; dT4_conjMu; dT4_vnormv; dT4_u; dT4_conju];
dCVarMain2(:,5) = [daM_Mu; daM_conjMu; daM_vnormv; daM_u; daM_conju];
dCVarMain2(:,6) = [daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju];
dCVarMain2(:,7) = [drho_Mu; drho_conjMu; drho_vnormv; drho_u; drho_conju];
dCVarMain2(:,8) = [d2w_uMu; d2w_uconjMu; d2w_uvnormv; d2w_uu; 0];
dCVarMain2(:,9) = [d2T4_uMu; d2T4_uconjMu; d2T4_uvnormv; d2T4_uu; d2T4_uconju];

% Lie bracket: lieTwo in [u1, u2, u3, u4, u5, u6];
% u1 = X1, u2 = conj_X1, u3 = X2, u4 = conj_X2, u5 = T, u6 = d/dgamma.
lieTwo = cell(6,6);
% [X1,conj_X2] = (-2*u/(1+u*conj(u)))*X1 + 2*T4*X2
%   + (2*u*conj(w)/(1+u*conj(u)) + 2*conj(T4) - conj(dw_u))*conj_X2 - 2T.
lieTwo_X1_conj_X2 = [-2*u/(1+u*conj(u)); 0; 2*T4; ...
    2*u*conj(w)/(1+u*conj(u)) + 2*conj(T4) - conj(dw_u); -2; 0];
lieTwo_X1_conj_X1 = [conj(aM); -aM; aM*w + 2*i*h11*T4 + dw_conjMu;...
    -conj(aM)*conj(w) + 2*i*h11*conj(T4) - conj(dw_conjMu); -2*i*h11; 0];
lieTwo_conj_X1_T = [bV - w/(1+u*conj(u))^2; aV - 2*conj(u)/(1+u*conj(u))*T4; ...
    -dw_vnormv - T4*dw_u + (-aV + dT4_u + 2*conj(u)/(1+u*conj(u))*T4)*w + dT4_Mu;...
    w*conj(w)/(1+u*conj(u))^2 - bV*conj(w) + conj(dT4_conjMu) + w*conj(dT4_conju);...
    0; 0];
lieTwo_X2_T = [-1/(1+u*conj(u))^2; 0; dT4_u; conj(w)/(1+u*conj(u))^2 + conj(dT4_conju); ...
    0; 0];
lieTwo_X1_T = [conj(lieTwo_conj_X1_T(2)); conj(lieTwo_conj_X1_T(1));...
    conj(lieTwo_conj_X1_T(4)); conj(lieTwo_conj_X1_T(3));...
    conj(lieTwo_conj_X1_T(5)); conj(lieTwo_conj_X1_T(6))];
lieTwo_conj_X1_X2 = [conj(lieTwo_X1_conj_X2(2)); conj(lieTwo_X1_conj_X2(1));...
    conj(lieTwo_X1_conj_X2(4)); conj(lieTwo_X1_conj_X2(3));
    conj(lieTwo_X1_conj_X2(5)); conj(lieTwo_X1_conj_X2(6))];
lieTwo_conj_X2_T = [conj(lieTwo_X2_T(2)); conj(lieTwo_X2_T(1));
    conj(lieTwo_X2_T(4)); conj(lieTwo_X2_T(3));
    conj(lieTwo_X2_T(5)); conj(lieTwo_X2_T(6))];

lieTwo{1,1} = zeros(6,1);
lieTwo{1,2} = lieTwo_X1_conj_X1;
lieTwo{1,3} = zeros(6,1);
lieTwo{1,4} = lieTwo_X1_conj_X2;
lieTwo{1,5} = lieTwo_X1_T;
lieTwo{1,6} = zeros(6,1);
lieTwo{2,1} = -lieTwo_X1_conj_X1;
lieTwo{2,2} = zeros(6,1);
lieTwo{2,3} = lieTwo_conj_X1_X2;
lieTwo{2,4} = zeros(6,1);
lieTwo{2,5} = lieTwo_conj_X1_T;
lieTwo{2,6} = zeros(6,1);
lieTwo{3,1} = zeros(6,1);
lieTwo{3,2} = -lieTwo_conj_X1_X2;
lieTwo{3,3} = zeros(6,1);
lieTwo{3,4} = zeros(6,1);
lieTwo{3,5} = lieTwo_X2_T;
lieTwo{3,6} = zeros(6,1);
lieTwo{4,1} = -lieTwo_X1_conj_X2;
lieTwo{4,2} = zeros(6,1);
lieTwo{4,3} = zeros(6,1);
lieTwo{4,4} = zeros(6,1);
lieTwo{4,5} = lieTwo_conj_X2_T;
lieTwo{4,6} = zeros(6,1);
lieTwo{5,1} = -lieTwo_X1_T;
lieTwo{5,2} = -lieTwo_conj_X1_T;
lieTwo{5,3} = -lieTwo_X2_T;
lieTwo{5,4} = -lieTwo_conj_X2_T;
lieTwo{5,5} = zeros(6,1);
lieTwo{5,6} = zeros(6,1);
lieTwo{6,1} = zeros(6,1);
lieTwo{6,2} = zeros(6,1);
lieTwo{6,3} = zeros(6,1);
lieTwo{6,4} = zeros(6,1);
lieTwo{6,5} = zeros(6,1);
lieTwo{6,6} = zeros(6,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task: Find christoffel symbols in basis [X1, conj_X1, X2, conj_X2, T, d/d\gamma].
% Uvector = [X1, conj(X1), X2, conj(X2), T, d/dgamma] 
% in (Mu, conjMu, vnormv, u, conj(u), gamma). 
Uvector = [0, 1, 0, 0, 0, 0;
           1, 0, 0, 0, 0, 0;
           0, 0, 0, 0, 1, 0;
           0, w, 1, 0, T4, 0;
           conj(w), 0, 0, 1, conj(T4), 0;
           0, 0, 0, 0, 0, 1];
       
RGamma0 = sym('RGamma0',[6 6 6]);
dF0 = cell(6,6);
for j=1:6
    for k=1:6
        dF0{j,k} = df_NatTwo_MuGamma_CR_rmW(F0(j,k),CVarMain2, dCVarMain2, gamma);
    end
end
% RGamma0(m,n,k) = \Phi_{mn,k} on P.43
% RGamma(m,n,k) = \Phi_{mn}^k on P.43
for m=1:6
    for n=1:6
        for k=1:6
            temp1 = dF0{n,k}*Uvector(:,m) - dF0{m,n}*Uvector(:,k)...
                + dF0{m,k}*Uvector(:,n);
            temp2 = 0;
            for ll=1:6
                % Change on Apr 20 2019.
                temp2 = temp2 -F0(m,ll)*lieTwo{n,k}(ll)...
                    + F0(k,ll)*lieTwo{m,n}(ll) + F0(n,ll)*lieTwo{k,m}(ll);
            end
            RGamma0(m,n,k) = 1/2*(temp1+ temp2);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RGamma = sym('RGamma', [6 6 6]);
symSetRGamma = [];
for m=1:6
    for n=1:6
        for k=1:6
            temp3 = 0;
            for ll=1:6
                temp3 = temp3 + RGamma0(m,n,ll)*F0Inv(ll,k);
            end
            RGamma(m,n,k) = temp3;
            tempSet = symvar(temp3);
            symSetRGamma = union(symSetRGamma, tempSet);
        end
    end
end

for m=1:6
    for n=1:6
        for k=1:6
            RGamma(m,n,k) = complex_simple3(RGamma(m,n,k), symSetRGamma);
        end
    end
end

% save('DataTempMain2_NatTwoR_CR_rmW.mat');

clearvars m n k ll j temp1 temp2 temp3 tempSet
clearvars RGamma0

list_of_variables_Main2 = who;
save('DataMain2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

