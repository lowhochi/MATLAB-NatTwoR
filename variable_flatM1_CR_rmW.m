% variable_flatM1_CR_rmW.m
% update lieMu for flat-case
clear lieMu
clearvars lie_Mu_conjMu lie_Mu_vnormv lie_Mu_ddu lie_conjMu_vnormv
clearvars lie_conjMu_ddconju lie_vnormv_ddu lie_vnormv_ddconju
lie_Mu_conjMu = [0; 0; 0];
lie_Mu_vnormv = [0; 0; 0];
lie_Mu_ddu = [-2*conj(u)/(1+u*conj(u)); 0; -2];
lie_conjMu_vnormv = [0; 0; 0];
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

syms Y
syms GH11_1 GH11_2 GT01_2 phiW %check that phiW is not previously defined.

syms dGH11_1_Mu dGH11_1_conjMu dGH11_1_vnormv dGH11_1_conju
dGH11_1_u = -2*conj(w)/y^2;
syms d2GH11_1_MuMu d2GH11_1_MuconjMu d2GH11_1_Muvnormv d2GH11_1_conjMuconjMu
syms d2GH11_1_conjMuvnormv d2GH11_1_vnormvvnormv
syms d2GH11_1_conjuMu d2GH11_1_conjuconjMu d2GH11_1_conjuvnormv
dGH11_1Row = [dGH11_1_Mu, dGH11_1_conjMu, dGH11_1_vnormv];
d2GH11_1_conjMuMu = d2GH11_1_MuconjMu + dGH11_1Row*lieMu{1,2};
d2GH11_1_vnormvMu = d2GH11_1_Muvnormv + dGH11_1Row*lieMu{1,3};
d2GH11_1_vnormvconjMu = d2GH11_1_conjMuvnormv + dGH11_1Row*lieMu{2,3};
d2GH11_1_Muconju = d2GH11_1_conjuMu + dGH11_1Row*lieMu{5,1};
d2GH11_1_conjMuconju = d2GH11_1_conjuconjMu + dGH11_1Row*lieMu{5,2};
d2GH11_1_vnormvconju = d2GH11_1_conjuvnormv + dGH11_1Row*lieMu{5,3}; 

MVarFlat12 = [u, w, dw_u, dw_Mu, dw_vnormv, d2rho_conjuconjMu,...
    d2rho_conjuconju, d2rho_uMu, d2rho_uu, drho_Mu, drho_conjMu,...
    drho_conju, drho_u, GH11_1, GH11_2, phiW,...
    dGH11_1_Mu, dGH11_1_conjMu, dGH11_1_vnormv, dGH11_1_conju,...
    d2GH11_1_MuMu, d2GH11_1_MuconjMu, d2GH11_1_Muvnormv, d2GH11_1_conjMuconjMu,...
    d2GH11_1_conjMuvnormv, d2GH11_1_vnormvvnormv,...
    d2GH11_1_conjuMu, d2GH11_1_conjuconjMu, d2GH11_1_conjuvnormv, Y];
