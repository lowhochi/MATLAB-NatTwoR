% variable_NatTwoR_RicCurv_CR_rmW.m

% d3aV_uXX
syms d3aV_uconjuconju d3aV_uuconju d3aV_uuu
syms d3aV_uMuMu d3aV_uMuconjMu d3aV_uMuvnormv d3aV_uconjMuconjMu 
syms d3aV_uconjMuvnormv d3aV_uvnormvvnormv
syms d3aV_uuMu d3aV_uuconjMu d3aV_uuvnormv 
syms d3aV_uconjuMu d3aV_uconjuconjMu d3aV_uconjuvnormv
d2aV_uRow = [d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv];
d3aV_uconjMuMu = d3aV_uMuconjMu + d2aV_uRow*lieMu{1,2};
d3aV_uvnormvMu = d3aV_uMuvnormv + d2aV_uRow*lieMu{1,3};
d3aV_uvnormvconjMu = d3aV_uconjMuvnormv + d2aV_uRow*lieMu{2,3};
d3aV_uMuu = d3aV_uuMu + d2aV_uRow*lieMu{4,1};
d3aV_uconjMuu = d3aV_uuconjMu + d2aV_uRow*lieMu{4,2};
d3aV_uvnormvu = d3aV_uuvnormv + d2aV_uRow*lieMu{4,3};
d3aV_uMuconju = d3aV_uconjuMu + d2aV_uRow*lieMu{5,1};
d3aV_uconjMuconju = d3aV_uconjuconjMu + d2aV_uRow*lieMu{5,2};
d3aV_uvnormvconju = d3aV_uconjuvnormv + d2aV_uRow*lieMu{5,3}; 

% d3aV_conjuconjuX
syms d3aV_conjuconjuMu d3aV_conjuconjuconjMu d3aV_conjuconjuvnormv
% exist d3aV_uconjuconju
syms d3aV_conjuconjuconju

% d4aV_uconjuXX
syms d4aV_uconjuMuMu d4aV_uconjuMuconjMu d4aV_uconjuMuvnormv
syms d4aV_uconjuconjMuconjMu d4aV_uconjuconjMuvnormv d4aV_uconjuvnormvvnormv
syms d4aV_uuconjuMu d4aV_uuconjuconjMu d4aV_uuconjuvnormv
syms d4aV_uconjuconjuMu d4aV_uconjuconjuconjMu d4aV_uconjuconjuvnormv
syms d4aV_uuuconju d4aV_uuconjuconju d4aV_uconjuconjuconju

d3aV_uconjuRow = [d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv];
d4aV_uconjuconjMuMu = d4aV_uconjuMuconjMu + d3aV_uconjuRow*lieMu{1,2};
d4aV_uconjuvnormvMu = d4aV_uconjuMuvnormv + d3aV_uconjuRow*lieMu{1,3};
d4aV_uconjuvnormvconjMu = d4aV_uconjuconjMuvnormv + d3aV_uconjuRow*lieMu{2,3};
d4aV_uconjuMuu = d4aV_uuconjuMu + d3aV_uconjuRow*lieMu{4,1};
d4aV_uconjuconjMuu = d4aV_uuconjuconjMu + d3aV_uconjuRow*lieMu{4,2};
d4aV_uconjuvnormvu = d4aV_uuconjuvnormv + d3aV_uconjuRow*lieMu{4,3};
d4aV_uconjuMuconju = d4aV_uconjuconjuMu + d3aV_uconjuRow*lieMu{5,1};
d4aV_uconjuconjMuconju = d4aV_uconjuconjuconjMu + d3aV_uconjuRow*lieMu{5,2};
d4aV_uconjuvnormvconju = d4aV_uconjuconjuvnormv + d3aV_uconjuRow*lieMu{5,3};

% d4w_uuuX
syms d4w_uuuMu d4w_uuuconjMu d4w_uuuvnormv d4w_uuuu

% dtheta
syms dtheta_Mu dtheta_conjMu dtheta_vnormv

CVarRicCurv = [u, aV, daV_u, daV_conju, d2aV_uu, d2aV_uconju,...
    w, dw_u, d2w_uu, theta,...
    d2aV_conjuconju, d3aV_uuconju, d3aV_uconjuconju, d3w_uuu];

myNumber = length(CVarRicCurv);
dCVarRicCurv = sym('dCVar',[5,myNumber]);

dCVarRicCurv(:,1)=[0; 0; 0; 1; 0]; %u
dCVarRicCurv(:,2)=[daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVarRicCurv(:,3)=[d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv; d2aV_uu; d2aV_uconju]; %daV_u
dCVarRicCurv(:,4)=[d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;...
    d2aV_uconju; d2aV_conjuconju]; %daV_conju
dCVarRicCurv(:,5)=[d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv; d3aV_uuu;...
    d3aV_uuconju]; %d2aV_uu
dCVarRicCurv(:,6)=[d3aV_uconjuMu; d3aV_uconjuconjMu; d3aV_uconjuvnormv;...
    d3aV_uuconju; d3aV_uconjuconju];%d2aV_uconju
dCVarRicCurv(:,7)=[dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0]; %w
dCVarRicCurv(:,8)=[d2w_uMu; d2w_uconjMu; d2w_uvnormv; d2w_uu; 0]; %dw_u
dCVarRicCurv(:,9)=[d3w_uuMu; d3w_uuconjMu; d3w_uuvnormv; d3w_uuu; 0]; %d2w_uu
dCVarRicCurv(:,10)=[dtheta_Mu; dtheta_conjMu; dtheta_vnormv;
    0; 0]; %theta
% %
dCVarRicCurv(:,11)=[d3aV_conjuconjuMu; d3aV_conjuconjuconjMu;
    d3aV_conjuconjuvnormv; d3aV_uconjuconju; d3aV_conjuconjuconju]; %d2aV_conjuconju
dCVarRicCurv(:,12)=[d4aV_uuconjuMu; d4aV_uuconjuconjMu; d4aV_uuconjuvnormv;
    d4aV_uuuconju; d4aV_uuconjuconju]; %d3aV_uuconju
dCVarRicCurv(:,13)=[d4aV_uconjuconjuMu; d4aV_uconjuconjuconjMu;
    d4aV_uconjuconjuvnormv; d4aV_uuconjuconju; d4aV_uconjuconjuconju]; %d3aV_uconjuconju
dCVarRicCurv(:,14)=[d4w_uuuMu; d4w_uuuconjMu; d4w_uuuvnormv; d4w_uuuu; 0]; %d3w_uuu
% %

