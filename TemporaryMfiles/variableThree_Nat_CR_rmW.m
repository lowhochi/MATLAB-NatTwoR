% variableThree_Nat_CR_rmW.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Mu, conjMu] and other Lie brackets in (Mu, conjMu, vnormv)
lie_Mu_conjMu = [aMu; -conj(aMu); 2*i*h11];
lie_Mu_vnormv = [aV; bV; 2*T4];
lie_Mu_ddu = [-2*conj(u)/(1+u*conj(u)); 0; -2];
lie_conjMu_vnormv = [conj(bV); conj(aV); 2*conj(T4)];
lie_conjMu_ddconju = [0; -2*u/(1+u*conj(u)); -2];
lie_vnormv_ddu = [0; 1/(1+u*conj(u))^2; 0];
lie_vnormv_ddconju = [1/(1+u*conj(u))^2; 0; 0];

lieMu = cell(5,5);
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
% Part(1) w-terms
% syms d2w_uu d2w_uMu d2w_uconjMu d2w_uvnormv (Already defined)
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
% Part(2) h11-terms 
% syms d2h11_conjuconju d2h11_uconju d2h11_uu (Already defined)
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

% d2w_vnormvu = d2w_uvnormv - dwRow*lieMu{3,4};
d2h11 = [d2h11_MuMu, d2h11_MuconjMu, d2h11_Muvnormv, d2h11_Muu, d2h11_Muconju;
    d2h11_conjMuMu, d2h11_conjMuconjMu, d2h11_conjMuvnormv, ...
        d2h11_conjMuu, d2h11_conjMuconju;
    d2h11_vnormvMu, d2h11_vnormvconjMu, d2h11_vnormvvnormv, ...
        d2h11_vnormvu, d2h11_vnormvconju;
    d2h11_uMu, d2h11_uconjMu, d2h11_uvnormv, d2h11_uu, d2h11_uconju;
    d2h11_conjuMu, d2h11_conjuconjMu, d2h11_conjuvnormv, d2h11_uconju, d2h11_conjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(3) T4-terms
% syms d2T4_uMu d2T4_uconjMu d2T4_uvnormv d2T4_uu d2T4_uconju
syms d2T4_MuMu d2T4_MuconjMu d2T4_Muvnormv d2T4_conjMuconjMu
syms d2T4_conjMuvnormv d2T4_vnormvvnormv
syms d2T4_conjuMu d2T4_conjuconjMu d2T4_conjuvnormv 
syms d2T4_conjuconju
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
% Part(4) aMu-terms
syms d2aMu_MuMu d2aMu_MuconjMu d2aMu_Muvnormv d2aMu_conjMuconjMu
syms d2aMu_conjMuvnormv d2aMu_vnormvvnormv
syms d2aMu_uMu d2aMu_uconjMu d2aMu_uvnormv 
syms d2aMu_conjuMu d2aMu_conjuconjMu d2aMu_conjuvnormv 
syms d2aMu_uu d2aMu_uconju d2aMu_conjuconju
daMuRow = [daMu_Mu, daMu_conjMu, daMu_vnormv];
d2aMu_conjMuMu = d2aMu_MuconjMu + daMuRow*lieMu{1,2};
d2aMu_vnormvMu = d2aMu_Muvnormv + daMuRow*lieMu{1,3};
d2aMu_vnormvconjMu = d2aMu_conjMuvnormv + daMuRow*lieMu{2,3};
d2aMu_Muu = d2aMu_uMu + daMuRow*lieMu{4,1};
d2aMu_conjMuu = d2aMu_uconjMu + daMuRow*lieMu{4,2};
d2aMu_vnormvu = d2aMu_uvnormv + daMuRow*lieMu{4,3};
d2aMu_Muconju = d2aMu_conjuMu + daMuRow*lieMu{5,1};
d2aMu_conjMuconju = d2aMu_conjuconjMu + daMuRow*lieMu{5,2};
d2aMu_vnormvconju = d2aMu_conjuvnormv + daMuRow*lieMu{5,3}; 

d2aMu = [d2aMu_MuMu, d2aMu_MuconjMu, d2aMu_Muvnormv, d2aMu_Muu, d2aMu_Muconju;
    d2aMu_conjMuMu, d2aMu_conjMuconjMu, d2aMu_conjMuvnormv, ...
        d2aMu_conjMuu, d2aMu_conjMuconju;
    d2aMu_vnormvMu, d2aMu_vnormvconjMu, d2aMu_vnormvvnormv, ...
        d2aMu_vnormvu, d2aMu_vnormvconju;
    d2aMu_uMu, d2aMu_uconjMu, d2aMu_uvnormv, d2aMu_uu, d2aMu_uconju;
    d2aMu_conjuMu, d2aMu_conjuconjMu, d2aMu_conjuvnormv, d2aMu_uconju, d2aMu_conjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(5) aV-terms 
syms d2aV_MuMu d2aV_MuconjMu d2aV_Muvnormv d2aV_conjMuconjMu
syms d2aV_conjMuvnormv d2aV_vnormvvnormv
syms d2aV_uMu d2aV_uconjMu d2aV_uvnormv 
syms d2aV_conjuMu d2aV_conjuconjMu d2aV_conjuvnormv 
syms d2aV_uu d2aV_uconju d2aV_conjuconju
daVRow = [daV_Mu, daV_conjMu, daV_vnormv];
d2aV_conjMuMu = d2aV_MuconjMu + daVRow*lieMu{1,2};
d2aV_vnormvMu = d2aV_Muvnormv + daVRow*lieMu{1,3};
d2aV_vnormvconjMu = d2aV_conjMuvnormv + daVRow*lieMu{2,3};
d2aV_Muu = d2aV_uMu + daVRow*lieMu{4,1};
d2aV_conjMuu = d2aV_uconjMu + daVRow*lieMu{4,2};
d2aV_vnormvu = d2aV_uvnormv + daVRow*lieMu{4,3};
d2aV_Muconju = d2aV_conjuMu + daVRow*lieMu{5,1};
d2aV_conjMuconju = d2aV_conjuconjMu + daVRow*lieMu{5,2};
d2aV_vnormvconju = d2aV_conjuvnormv + daVRow*lieMu{5,3}; 
d2aV = [d2aV_MuMu, d2aV_MuconjMu, d2aV_Muvnormv, d2aV_Muu, d2aV_Muconju;
    d2aV_conjMuMu, d2aV_conjMuconjMu, d2aV_conjMuvnormv, ...
        d2aV_conjMuu, d2aV_conjMuconju;
    d2aV_vnormvMu, d2aV_vnormvconjMu, d2aV_vnormvvnormv, ...
        d2aV_vnormvu, d2aV_vnormvconju;
    d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv, d2aV_uu, d2aV_uconju;
    d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv, d2aV_uconju, d2aV_conjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(6) rho-terms
syms d2rho_MuMu d2rho_MuconjMu d2rho_Muvnormv d2rho_conjMuconjMu
syms d2rho_conjMuvnormv d2rho_vnormvvnormv
syms d2rho_uMu d2rho_uconjMu d2rho_uvnormv 
syms d2rho_conjuMu d2rho_conjuconjMu d2rho_conjuvnormv 
syms d2rho_uu d2rho_uconju d2rho_conjuconju

drhoRow = [drho_Mu, drho_conjMu, drho_vnormv];
d2rho_conjMuMu = d2rho_MuconjMu + drhoRow*lieMu{1,2};
d2rho_vnormvMu = d2rho_Muvnormv + drhoRow*lieMu{1,3};
d2rho_vnormvconjMu = d2rho_conjMuvnormv + drhoRow*lieMu{2,3};
d2rho_Muu = d2rho_uMu + drhoRow*lieMu{4,1};
d2rho_conjMuu = d2rho_uconjMu + drhoRow*lieMu{4,2};
d2rho_vnormvu = d2rho_uvnormv + drhoRow*lieMu{4,3};
d2rho_Muconju = d2rho_conjuMu + drhoRow*lieMu{5,1};
d2rho_conjMuconju = d2rho_conjuconjMu + drhoRow*lieMu{5,2};
d2rho_vnormvconju = d2rho_conjuvnormv + drhoRow*lieMu{5,3}; 

d2rho = [d2rho_MuMu, d2rho_MuconjMu, d2rho_Muvnormv, d2rho_Muu, d2rho_Muconju;
    d2rho_conjMuMu, d2rho_conjMuconjMu, d2rho_conjMuvnormv, ...
        d2rho_conjMuu, d2rho_conjMuconju;
    d2rho_vnormvMu, d2rho_vnormvconjMu, d2rho_vnormvvnormv, ...
        d2rho_vnormvu, d2rho_vnormvconju;
    d2rho_uMu, d2rho_uconjMu, d2rho_uvnormv, d2rho_uu, d2rho_uconju;
    d2rho_conjuMu, d2rho_conjuconjMu, d2rho_conjuvnormv, d2rho_uconju, d2rho_conjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(7) dw_u terms 
d2wRow = [d2w_uMu, d2w_uconjMu, d2w_uvnormv];
syms d3w_uMuMu d3w_uMuconjMu d3w_uMuvnormv 
syms d3w_uconjMuconjMu d3w_uconjMuvnormv d3w_uvnormvvnormv
syms d3w_uuMu  d3w_uuconjMu d3w_uuvnormv d3w_uuu

d3w_uconjMuMu = d3w_uMuconjMu + d2wRow*lieMu{1,2};
d3w_uvnormvMu = d3w_uMuvnormv + d2wRow*lieMu{1,3};
d3w_uvnormvconjMu = d3w_uconjMuvnormv + d2wRow*lieMu{2,3};
d3w_uMuu = d3w_uuMu + d2wRow*lieMu{4,1};
d3w_uconjMuu = d3w_uuconjMu + d2wRow*lieMu{4,2};
d3w_uvnormvu = d3w_uuvnormv + d2wRow*lieMu{4,3};
d3w_uconjMuconju = d2wRow*lieMu{5,2};
d3w_uvnormvconju = d2wRow*lieMu{5,3};

d3w_u = [d3w_uMuMu, d3w_uMuconjMu, d3w_uMuvnormv, d3w_uMuu,  0;
    d3w_uconjMuMu, d3w_uconjMuconjMu, d3w_uconjMuvnormv, d3w_uconjMuu, d3w_uconjMuconju;
    d3w_uvnormvMu, d3w_uvnormvconjMu, d3w_uvnormvvnormv, d3w_uvnormvu, d3w_uvnormvconju;
    d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv, d3w_uuu, 0;
    0, 0, 0, 0, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(8) dT4_u terms 
d2T4Row = [d2T4_uMu, d2T4_uconjMu, d2T4_uvnormv];
syms d3T4_uuu d3T4_uuconju d3T4_uconjuconju
syms d3T4_uMuMu d3T4_uMuconjMu d3T4_uMuvnormv
syms d3T4_uconjMuconjMu d3T4_uconjMuvnormv d3T4_uvnormvvnormv
syms d3T4_uuMu d3T4_uuconjMu d3T4_uuvnormv 
syms d3T4_uconjuMu d3T4_uconjuconjMu d3T4_uconjuvnormv

d3T4_uconjMuMu = d3T4_uMuconjMu + d2T4Row*lieMu{1,2};
d3T4_uvnormvMu = d3T4_uMuvnormv + d2T4Row*lieMu{1,3};
d3T4_uvnormvconjMu = d3T4_uconjMuvnormv + d2T4Row*lieMu{2,3};
d3T4_uMuu = d3T4_uuMu + d2T4Row*lieMu{4,1};
d3T4_uconjMuu = d3T4_uuconjMu + d2T4Row*lieMu{4,2};
d3T4_uvnormvu = d3T4_uuvnormv + d2T4Row*lieMu{4,3};
d3T4_uMuconju = d3T4_uconjuMu + d2T4Row*lieMu{5,1};
d3T4_uconjMuconju = d3T4_uconjuconjMu + d2T4Row*lieMu{5,2};
d3T4_uvnormvconju = d3T4_uconjuvnormv + d2T4Row*lieMu{5,3};

d3T4_u = [d3T4_uMuMu, d3T4_uMuconjMu, d3T4_uMuvnormv, d3T4_uMuu, d3T4_uMuconju;
    d3T4_uconjMuMu, d3T4_uconjMuconjMu, d3T4_uconjMuvnormv,...
        d3T4_uconjMuu, d3T4_uconjMuconju;
    d3T4_uvnormvMu, d3T4_uvnormvconjMu, d3T4_uvnormvvnormv,...
        d3T4_uvnormvu, d3T4_uvnormvconju;
    d3T4_uuMu, d3T4_uuconjMu, d3T4_uuvnormv, d3T4_uuu, d3T4_uuconju;
    d3T4_uconjuMu, d3T4_uconjuconjMu, d3T4_uconjuvnormv,...
        d3T4_uuconju, d3T4_uconjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CVar4 and dCVar4 
CVar4 = [u, w, h11, T4, aMu, aV, bV, rho,...
    dw_u, dw_Mu, dw_conjMu, dw_vnormv,...
    dh11_u, dh11_conju, dh11_Mu, dh11_conjMu, dh11_vnormv,...
    dT4_u, dT4_conju, dT4_Mu, dT4_conjMu, dT4_vnormv,...
    daMu_u, daMu_conju, daMu_Mu, daMu_conjMu, daMu_vnormv,...
    daV_u, daV_conju, daV_Mu, daV_conjMu, daV_vnormv,...
    drho_u, drho_conju, drho_Mu, drho_conjMu, drho_vnormv,...
    d2w_uu, d2w_uMu, d2w_uconjMu, d2w_uvnormv,...
    d2T4_uu, d2T4_uconju, d2T4_uMu, d2T4_uconjMu, d2T4_uvnormv];

dCVar4 = sym('dCVar4', [5, 46]);
dCVar4(:,1) = [0;0;0;1;0]; %u
dCVar4(:,2) = [dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0];%w
dCVar4(:,3) = [dh11_Mu; dh11_conjMu; dh11_vnormv; dh11_u; dh11_conju]; %h11
dCVar4(:,4) = [dT4_Mu; dT4_conjMu; dT4_vnormv; dT4_u; dT4_conju]; %T4
dCVar4(:,5) = [daMu_Mu; daMu_conjMu; daMu_vnormv; daMu_u; daMu_conju]; %aMu
dCVar4(:,6) = [daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVar4(:,7) = [dbV_Mu; dbV_conjMu; dbV_vnormv; dbV_u; dbV_conju]; %bV
dCVar4(:,8) = [drho_Mu; drho_conjMu; drho_vnormv; drho_u; drho_conju]; %rho
%
dCVar4(:,9) = transpose(d2w(4,:)); %dw_u 
dCVar4(:,10) = transpose(d2w(1,:)); %dw_Mu
dCVar4(:,11) = transpose(d2w(2,:)); %dw_conjMu
dCVar4(:,12) = transpose(d2w(3,:)); %dw_vnormv
dCVar4(:,13) = transpose(d2h11(4,:)); %dh11_u
dCVar4(:,14) = transpose(d2h11(5,:)); %dh11_conju
dCVar4(:,15) = transpose(d2h11(1,:)); %dh11_Mu
dCVar4(:,16) = transpose(d2h11(2,:)); %dh11_conjMu
dCVar4(:,17) = transpose(d2h11(3,:)); %dh11_vnormv 
dCVar4(:,18) = transpose(d2T4(4,:)); %dT4_u
dCVar4(:,19) = transpose(d2T4(5,:)); %dT4_conju
dCVar4(:,20) = transpose(d2T4(1,:)); %dT4_Mu
dCVar4(:,21) = transpose(d2T4(2,:)); %dT4_conjMu
dCVar4(:,22) = transpose(d2T4(3,:)); %dT4_vnormv
dCVar4(:,23) = transpose(d2aMu(4,:)); %daMu_u
dCVar4(:,24) = transpose(d2aMu(5,:)); %daMu_conju
dCVar4(:,25) = transpose(d2aMu(1,:)); %daMu_Mu
dCVar4(:,26) = transpose(d2aMu(2,:)); %daMu_conjMu
dCVar4(:,27) = transpose(d2aMu(3,:)); %daMu_vnormv
dCVar4(:,28) = transpose(d2aV(4,:)); %daV_u
dCVar4(:,29) = transpose(d2aV(5,:)); %daV_conju
dCVar4(:,30) = transpose(d2aV(1,:)); %daV_Mu
dCVar4(:,31) = transpose(d2aV(2,:)); %daV_conjMu
dCVar4(:,32) = transpose(d2aV(3,:)); %daV_vnormv
dCVar4(:,33) = transpose(d2rho(4,:)); %drho_u
dCVar4(:,34) = transpose(d2rho(5,:)); %drho_conju
dCVar4(:,35) = transpose(d2rho(1,:)); %drho_Mu
dCVar4(:,36) = transpose(d2rho(2,:)); %drho_conjMu
dCVar4(:,37) = transpose(d2rho(3,:)); %drho_vnormv
%
dCVar4(:,38) = transpose(d3w_u(4,:)); %d2w_uu
dCVar4(:,39) = transpose(d3w_u(1,:)); %d2w_uMu
dCVar4(:,40) = transpose(d3w_u(2,:)); %d2w_uconjMu
dCVar4(:,41) = transpose(d3w_u(3,:)); %d2w_uvnormv
dCVar4(:,42) = transpose(d3T4_u(4,:)); %d2T4_uu
dCVar4(:,43) = transpose(d3T4_u(5,:)); %d2T4_uconju
dCVar4(:,44) = transpose(d3T4_u(1,:)); %d2T4_uMu
dCVar4(:,45) = transpose(d3T4_u(2,:)); %d2T4_uconjMu
dCVar4(:,46) = transpose(d3T4_u(3,:)); %d2T4_uvnormv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%