% variable_NatTwoR_Chern1_CR_rmW.m
syms d2aV_conjuconju d2aV_uconju d2aV_uu
syms d2aV_MuMu d2aV_MuconjMu d2aV_Muvnormv d2aV_conjMuconjMu 
syms d2aV_conjMuvnormv d2aV_vnormvvnormv
syms d2aV_uMu d2aV_uconjMu d2aV_uvnormv 
syms d2aV_conjuMu d2aV_conjuconjMu d2aV_conjuvnormv

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

CVarChern1 = sym('CVarChern1',[1,10]);
CVarChern1(1) = u;
CVarChern1(2) = aV;
CVarChern1(3) = daV_u;
CVarChern1(4) = daV_conju;
CVarChern1(5) = d2aV_uu;
CVarChern1(6) = d2aV_uconju;
CVarChern1(7) = d3aV_uuconju;
CVarChern1(8) = d3aV_uconjuMu;
CVarChern1(9) = d2aV_uMu;
CVarChern1(10) = daV_Mu;

dCVarChern1 = sym('dCVarChern1',[5,10]);

dCVarChern1(:,1) = [0;0;0;1;0];%u
dCVarChern1(:,2) = [daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVarChern1(:,3) = [d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv; d2aV_uu; d2aV_uconju]; %daV_u
dCVarChern1(:,4) = [d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;...
    d2aV_uconju; d2aV_conjuconju]; %daV_conju
dCVarChern1(:,5) = [d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv; d3aV_uuu;...
    d3aV_uuconju]; %d2aV_uu
dCVarChern1(:,6) = [d3aV_uconjuMu; d3aV_uconjuconjMu; d3aV_uconjuvnormv;...
    d3aV_uuconju; d3aV_uconjuconju];%d2aV_uconju
dCVarChern1(:,7) = [d4aV_uuconjuMu; d4aV_uuconjuconjMu; d4aV_uuconjuvnormv;...
    d4aV_uuuconju; d4aV_uuconjuconju];%d3aV_uuconju;
dCVarChern1(:,8) = [d4aV_uconjuMuMu; d4aV_uconjuMuconjMu; d4aV_uconjuMuvnormv;...
    d4aV_uconjuMuu; d4aV_uconjuMuconju]; %d3aV_uconjuMu;
dCVarChern1(:,9) = [d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv,...
    d3aV_uMuu, d3aV_uMuconju]; %d2aV_uMu;
dCVarChern1(:,10) = [d2aV_MuMu; d2aV_MuconjMu; d2aV_Muvnormv; d2aV_Muu;
    d2aV_Muconju]; %daV_Mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MVarChern1 = [u, aV, daV_Mu, daV_conjMu, daV_vnormv, daV_u, daV_conju,...
    d2aV_conjuconju, d2aV_uconju, d2aV_uu,...
    d2aV_MuMu, d2aV_MuconjMu, d2aV_Muvnormv, d2aV_conjMuconjMu,... 
    d2aV_conjMuvnormv, d2aV_vnormvvnormv, d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv,... 
	d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv,...
    d3aV_uconjuconju, d3aV_uuconju, d3aV_uuu,...
    d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv, d3aV_uconjMuconjMu,... 
    d3aV_uconjMuvnormv, d3aV_uvnormvvnormv,...
    d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv,... 
    d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv,...
    d4aV_uconjuMuMu, d4aV_uconjuMuconjMu, d4aV_uconjuMuvnormv,...
    d4aV_uconjuconjMuconjMu, d4aV_uconjuconjMuvnormv, d4aV_uconjuvnormvvnormv,...
    d4aV_uuconjuMu, d4aV_uuconjuconjMu, d4aV_uuconjuvnormv,...
    d4aV_uconjuconjuMu, d4aV_uconjuconjuconjMu, d4aV_uconjuconjuvnormv,...
    d4aV_uuuconju, d4aV_uuconjuconju, d4aV_uconjuconjuconju,...
    d2w_uu, d2w_uconjMu, d2w_conjMuconjMu, dw_u, dw_conjMu,...
    dw_vnormv, w];
