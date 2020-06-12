syms G12_3 G23_1 G31_2 G11_2 G11_3 G22_1 G22_3 G33_1 G33_2 % real
syms dG12_3_Mu dG12_3_conjMu dG12_3_vnormv
syms dG23_1_Mu dG23_1_conjMu dG23_1_vnormv
syms dG31_2_Mu dG31_2_conjMu dG31_2_vnormv 
syms dG11_2_Mu dG11_2_conjMu dG11_2_vnormv 
syms dG11_3_Mu dG11_3_conjMu dG11_3_vnormv 
syms dG22_1_Mu dG22_1_conjMu dG22_1_vnormv 
syms dG22_3_Mu dG22_3_conjMu dG22_3_vnormv
syms dG33_1_Mu dG33_1_conjMu dG33_1_vnormv
syms dG33_2_Mu dG33_2_conjMu dG33_2_vnormv
% The second derivatives
% d2G12_3
syms d2G12_3_MuMu d2G12_3_MuconjMu d2G12_3_Muvnormv
syms d2G12_3_conjMuconjMu d2G12_3_conjMuvnormv d2G12_3_vnormvvnormv
dG12_3Row = [dG12_3_Mu, dG12_3_conjMu, dG12_3_vnormv];
d2G12_3_conjMuMu = d2G12_3_MuconjMu + dG12_3Row*lieMu{1,2};
d2G12_3_vnormvMu = d2G12_3_Muvnormv + dG12_3Row*lieMu{1,3};
d2G12_3_vnormvconjMu = d2G12_3_conjMuvnormv + dG12_3Row*lieMu{2,3}; 
%
d2G12_3_Muu = dG12_3Row*lieMu{4,1};
d2G12_3_conjMuconju = conj(d2G12_3_Muu);
d2G12_3_vnormvu = dG12_3Row*lieMu{4,3};
d2G12_3_vnormvconju = dG12_3Row*lieMu{5,3};
d2G12_3 =[d2G12_3_MuMu, d2G12_3_MuconjMu, d2G12_3_Muvnormv,... 
    d2G12_3_Muu, 0;
    d2G12_3_conjMuMu, d2G12_3_conjMuconjMu, d2G12_3_conjMuvnormv,...
    0, d2G12_3_conjMuconju;
    d2G12_3_vnormvMu, d2G12_3_vnormvconjMu, d2G12_3_vnormvvnormv,...
    d2G12_3_vnormvu, d2G12_3_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G23_1
syms d2G23_1_MuMu d2G23_1_MuconjMu d2G23_1_Muvnormv
syms d2G23_1_conjMuconjMu d2G23_1_conjMuvnormv d2G23_1_vnormvvnormv
dG23_1Row = [dG23_1_Mu, dG23_1_conjMu, dG23_1_vnormv];
d2G23_1_conjMuMu = d2G23_1_MuconjMu + dG23_1Row*lieMu{1,2};
d2G23_1_vnormvMu = d2G23_1_Muvnormv + dG23_1Row*lieMu{1,3};
d2G23_1_vnormvconjMu = d2G23_1_conjMuvnormv + dG23_1Row*lieMu{2,3}; 
d2G23_1_Muu = dG23_1Row*lieMu{4,1};
d2G23_1_conjMuconju = conj(d2G23_1_Muu);
d2G23_1_vnormvu = dG23_1Row*lieMu{4,3};
d2G23_1_vnormvconju = dG23_1Row*lieMu{5,3};
d2G23_1 =[d2G23_1_MuMu, d2G23_1_MuconjMu, d2G23_1_Muvnormv,...
    d2G23_1_Muu, 0;
    d2G23_1_conjMuMu, d2G23_1_conjMuconjMu, d2G23_1_conjMuvnormv,...
    0, d2G23_1_conjMuconju;
    d2G23_1_vnormvMu, d2G23_1_vnormvconjMu, d2G23_1_vnormvvnormv,...
    d2G23_1_vnormvu, d2G23_1_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G31_2
syms d2G31_2_MuMu d2G31_2_MuconjMu d2G31_2_Muvnormv
syms d2G31_2_conjMuconjMu d2G31_2_conjMuvnormv d2G31_2_vnormvvnormv 
dG31_2Row = [dG31_2_Mu, dG31_2_conjMu, dG31_2_vnormv];
d2G31_2_conjMuMu = d2G31_2_MuconjMu + dG31_2Row*lieMu{1,2};
d2G31_2_vnormvMu = d2G31_2_Muvnormv + dG31_2Row*lieMu{1,3};
d2G31_2_vnormvconjMu = d2G31_2_conjMuvnormv + dG31_2Row*lieMu{2,3}; 
d2G31_2_Muu = dG31_2Row*lieMu{4,1};
d2G31_2_conjMuconju = conj(d2G31_2_Muu);
d2G31_2_vnormvu = dG31_2Row*lieMu{4,3};
d2G31_2_vnormvconju = dG31_2Row*lieMu{5,3};
d2G31_2 =[d2G31_2_MuMu, d2G31_2_MuconjMu, d2G31_2_Muvnormv,...
    d2G31_2_Muu, 0;
    d2G31_2_conjMuMu, d2G31_2_conjMuconjMu, d2G31_2_conjMuvnormv,...
    0, d2G31_2_conjMuconju;
    d2G31_2_vnormvMu, d2G31_2_vnormvconjMu, d2G31_2_vnormvvnormv,...
    d2G31_2_vnormvu, d2G31_2_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G11_2
syms d2G11_2_MuMu d2G11_2_MuconjMu d2G11_2_Muvnormv
syms d2G11_2_conjMuconjMu d2G11_2_conjMuvnormv d2G11_2_vnormvvnormv 
dG11_2Row = [dG11_2_Mu, dG11_2_conjMu, dG11_2_vnormv];
d2G11_2_conjMuMu = d2G11_2_MuconjMu + dG11_2Row*lieMu{1,2};
d2G11_2_vnormvMu = d2G11_2_Muvnormv + dG11_2Row*lieMu{1,3};
d2G11_2_vnormvconjMu = d2G11_2_conjMuvnormv + dG11_2Row*lieMu{2,3};
d2G11_2_Muu = dG11_2Row*lieMu{4,1};
d2G11_2_conjMuconju = conj(d2G11_2_Muu);
d2G11_2_vnormvu = dG11_2Row*lieMu{4,3};
d2G11_2_vnormvconju = dG11_2Row*lieMu{5,3};
d2G11_2 =[d2G11_2_MuMu, d2G11_2_MuconjMu, d2G11_2_Muvnormv,... 
    d2G11_2_Muu, 0;
    d2G11_2_conjMuMu, d2G11_2_conjMuconjMu, d2G11_2_conjMuvnormv,...
    0, d2G11_2_conjMuconju;
    d2G11_2_vnormvMu, d2G11_2_vnormvconjMu, d2G11_2_vnormvvnormv,...
    d2G11_2_vnormvu, d2G11_2_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G11_3
syms d2G11_3_MuMu d2G11_3_MuconjMu d2G11_3_Muvnormv
syms d2G11_3_conjMuconjMu d2G11_3_conjMuvnormv d2G11_3_vnormvvnormv 
dG11_3Row = [dG11_3_Mu, dG11_3_conjMu, dG11_3_vnormv];
d2G11_3_conjMuMu = d2G11_3_MuconjMu + dG11_3Row*lieMu{1,2};
d2G11_3_vnormvMu = d2G11_3_Muvnormv + dG11_3Row*lieMu{1,3};
d2G11_3_vnormvconjMu = d2G11_3_conjMuvnormv + dG11_3Row*lieMu{2,3};
d2G11_3_Muu = dG11_3Row*lieMu{4,1};
d2G11_3_conjMuconju = conj(d2G11_3_Muu);
d2G11_3_vnormvu = dG11_3Row*lieMu{4,3};
d2G11_3_vnormvconju = dG11_3Row*lieMu{5,3};
d2G11_3 =[d2G11_3_MuMu, d2G11_3_MuconjMu, d2G11_3_Muvnormv,...
    d2G11_3_Muu, 0;
    d2G11_3_conjMuMu, d2G11_3_conjMuconjMu, d2G11_3_conjMuvnormv,...
    0, d2G11_3_conjMuconju;
    d2G11_3_vnormvMu, d2G11_3_vnormvconjMu, d2G11_3_vnormvvnormv,...
    d2G11_3_vnormvu, d2G11_3_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G22_1
syms d2G22_1_MuMu d2G22_1_MuconjMu d2G22_1_Muvnormv
syms d2G22_1_conjMuconjMu d2G22_1_conjMuvnormv d2G22_1_vnormvvnormv
dG22_1Row = [dG22_1_Mu, dG22_1_conjMu, dG22_1_vnormv];
d2G22_1_conjMuMu = d2G22_1_MuconjMu + dG22_1Row*lieMu{1,2};
d2G22_1_vnormvMu = d2G22_1_Muvnormv + dG22_1Row*lieMu{1,3};
d2G22_1_vnormvconjMu = d2G22_1_conjMuvnormv + dG22_1Row*lieMu{2,3};
d2G22_1_Muu = dG22_1Row*lieMu{4,1};
d2G22_1_conjMuconju = conj(d2G22_1_Muu);
d2G22_1_vnormvu = dG22_1Row*lieMu{4,3};
d2G22_1_vnormvconju = dG22_1Row*lieMu{5,3};
d2G22_1 =[d2G22_1_MuMu, d2G22_1_MuconjMu, d2G22_1_Muvnormv,...
    d2G22_1_Muu, 0;
    d2G22_1_conjMuMu, d2G22_1_conjMuconjMu, d2G22_1_conjMuvnormv,...
    0, d2G22_1_conjMuconju;
    d2G22_1_vnormvMu, d2G22_1_vnormvconjMu, d2G22_1_vnormvvnormv,...
    d2G22_1_vnormvu, d2G22_1_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G22_3
syms d2G22_3_MuMu d2G22_3_MuconjMu d2G22_3_Muvnormv
syms d2G22_3_conjMuconjMu d2G22_3_conjMuvnormv d2G22_3_vnormvvnormv
dG22_3Row = [dG22_3_Mu, dG22_3_conjMu, dG22_3_vnormv];
d2G22_3_conjMuMu = d2G22_3_MuconjMu + dG22_3Row*lieMu{1,2};
d2G22_3_vnormvMu = d2G22_3_Muvnormv + dG22_3Row*lieMu{1,3};
d2G22_3_vnormvconjMu = d2G22_3_conjMuvnormv + dG22_3Row*lieMu{2,3}; 
d2G22_3_Muu = dG22_3Row*lieMu{4,1};
d2G22_3_conjMuconju = conj(d2G22_3_Muu);
d2G22_3_vnormvu = dG22_3Row*lieMu{4,3};
d2G22_3_vnormvconju = dG22_3Row*lieMu{5,3};
d2G22_3 =[d2G22_3_MuMu, d2G22_3_MuconjMu, d2G22_3_Muvnormv,...
    d2G22_3_Muu, 0;
    d2G22_3_conjMuMu, d2G22_3_conjMuconjMu, d2G22_3_conjMuvnormv,...
    0, d2G22_3_conjMuconju;
    d2G22_3_vnormvMu, d2G22_3_vnormvconjMu, d2G22_3_vnormvvnormv,...
    d2G22_3_vnormvu, d2G22_3_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G33_1
syms d2G33_1_MuMu d2G33_1_MuconjMu d2G33_1_Muvnormv
syms d2G33_1_conjMuconjMu d2G33_1_conjMuvnormv d2G33_1_vnormvvnormv
dG33_1Row = [dG33_1_Mu, dG33_1_conjMu, dG33_1_vnormv];
d2G33_1_conjMuMu = d2G33_1_MuconjMu + dG33_1Row*lieMu{1,2};
d2G33_1_vnormvMu = d2G33_1_Muvnormv + dG33_1Row*lieMu{1,3};
d2G33_1_vnormvconjMu = d2G33_1_conjMuvnormv + dG33_1Row*lieMu{2,3}; 
d2G33_1_Muu = dG33_1Row*lieMu{4,1};
d2G33_1_conjMuconju = conj(d2G33_1_Muu);
d2G33_1_vnormvu = dG33_1Row*lieMu{4,3};
d2G33_1_vnormvconju = dG33_1Row*lieMu{5,3};
d2G33_1 =[d2G33_1_MuMu, d2G33_1_MuconjMu, d2G33_1_Muvnormv,...
    d2G33_1_Muu, 0;
    d2G33_1_conjMuMu, d2G33_1_conjMuconjMu, d2G33_1_conjMuvnormv,...
    0, d2G33_1_conjMuconju;
    d2G33_1_vnormvMu, d2G33_1_vnormvconjMu, d2G33_1_vnormvvnormv,...
    d2G33_1_vnormvu, d2G33_1_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d2G33_2
syms d2G33_2_MuMu d2G33_2_MuconjMu d2G33_2_Muvnormv
syms d2G33_2_conjMuconjMu d2G33_2_conjMuvnormv d2G33_2_vnormvvnormv
dG33_2Row = [dG33_2_Mu, dG33_2_conjMu, dG33_2_vnormv];
d2G33_2_conjMuMu = d2G33_2_MuconjMu + dG33_2Row*lieMu{1,2};
d2G33_2_vnormvMu = d2G33_2_Muvnormv + dG33_2Row*lieMu{1,3};
d2G33_2_vnormvconjMu = d2G33_2_conjMuvnormv + dG33_2Row*lieMu{2,3}; 
d2G33_2_Muu = dG33_2Row*lieMu{4,1};
d2G33_2_conjMuconju = conj(d2G33_2_Muu);
d2G33_2_vnormvu = dG33_2Row*lieMu{4,3};
d2G33_2_vnormvconju = dG33_2Row*lieMu{5,3};
d2G33_2 = [d2G33_2_MuMu, d2G33_2_MuconjMu, d2G33_2_Muvnormv,... 
    d2G33_2_Muu, 0;
    d2G33_2_conjMuMu, d2G33_2_conjMuconjMu, d2G33_2_conjMuvnormv,...
    0, d2G33_2_conjMuconju;
    d2G33_2_vnormvMu, d2G33_2_vnormvconjMu, d2G33_2_vnormvvnormv,...
    d2G33_2_vnormvu, d2G33_2_vnormvconju;
    0, 0, 0, 0, 0;
    0, 0, 0, 0, 0];

% d3w_u is set up already.
% d4w
syms d4w_uuMuMu d4w_uuMuconjMu d4w_uuMuvnormv d4w_uuconjMuconjMu
syms d4w_uuconjMuconjMu d4w_uuconjMuvnormv d4w_uuvnormvvnormv 
syms d4w_uuuMu d4w_uuuconjMu d4w_uuuvnormv d4w_uuuu

d3w_uuRow = [d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv];
d4w_uuconjMuMu = d4w_uuMuconjMu + d3w_uuRow*lieMu{1,2};
d4w_uuvnormvMu = d4w_uuMuvnormv + d3w_uuRow*lieMu{1,3};
d4w_uuvnormvconjMu = d4w_uuconjMuvnormv + d3w_uuRow*lieMu{2,3};

d4w_uuMuu = d4w_uuuMu + d3w_uuRow*lieMu{4,1};
d4w_uuconjMuu = d4w_uuuconjMu + d3w_uuRow*lieMu{4,2};
d4w_uuvnormvu = d4w_uuuvnormv + d3w_uuRow*lieMu{4,3};
d4w_uuconjMuconju = d3w_uuRow*lieMu{5,2};
d4w_uuvnormvconju = d3w_uuRow*lieMu{5,3};
d4w_uu = [d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuMuvnormv, d4w_uuMuu, 0;
    d4w_uuconjMuMu, d4w_uuconjMuconjMu, d4w_uuconjMuvnormv, d4w_uuconjMuu, ...
    d4w_uuconjMuconju;
    d4w_uuvnormvMu, d4w_uuvnormvconjMu, d4w_uuvnormvvnormv, ...
    d4w_uuvnormvu, d4w_uuvnormvconju;
    d4w_uuuMu, d4w_uuuconjMu, d4w_uuuvnormv, d4w_uuuu, 0; 
    0, 0, 0, 0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVar5 = [G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3, G33_1, G33_2,...
    dG12_3_Mu, dG12_3_conjMu, dG12_3_vnormv,...
    dG23_1_Mu, dG23_1_conjMu, dG23_1_vnormv,...
    dG31_2_Mu, dG31_2_conjMu, dG31_2_vnormv,...
    dG11_2_Mu, dG11_2_conjMu, dG11_2_vnormv,...
    dG11_3_Mu, dG11_3_conjMu, dG11_3_vnormv,...
    dG22_1_Mu, dG22_1_conjMu, dG22_1_vnormv,...
    dG22_3_Mu, dG22_3_conjMu, dG22_3_vnormv,...
    dG33_1_Mu, dG33_1_conjMu, dG33_1_vnormv,...
    dG33_2_Mu, dG33_2_conjMu, dG33_2_vnormv,...
    u, w, dw_Mu, dw_conjMu, dw_vnormv, dw_u,...
    d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_uu,...
    d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv, d3w_uuu];

dCVar5 = sym('dCVar5', [5,50]);
% dCVar5 is in the format [Mu; conjMu; vnormv; du; dconju];
dCVar5(:,1) = [dG12_3_Mu; dG12_3_conjMu; dG12_3_vnormv; 0; 0]; %G12_3
dCVar5(:,2) = [dG23_1_Mu; dG23_1_conjMu; dG23_1_vnormv; 0; 0]; %G23_1
dCVar5(:,3) = [dG31_2_Mu; dG31_2_conjMu; dG31_2_vnormv; 0; 0]; %G31_2
dCVar5(:,4) = [dG11_2_Mu; dG11_2_conjMu; dG11_2_vnormv; 0; 0]; %G11_2
dCVar5(:,5) = [dG11_3_Mu; dG11_3_conjMu; dG11_3_vnormv; 0; 0]; %G11_3
dCVar5(:,6) = [dG22_1_Mu; dG22_1_conjMu; dG22_1_vnormv; 0; 0]; %G22_1
dCVar5(:,7) = [dG22_3_Mu; dG22_3_conjMu; dG22_3_vnormv; 0; 0]; %G22_3
dCVar5(:,8) = [dG33_1_Mu; dG33_1_conjMu; dG33_1_vnormv; 0; 0]; %G33_1
dCVar5(:,9) = [dG33_2_Mu; dG33_2_conjMu; dG33_2_vnormv; 0; 0]; %G33_2
%
dCVar5(:,10)= transpose(d2G12_3(1,:)); %dG12_3_Mu
dCVar5(:,11)= transpose(d2G12_3(2,:)); %dG12_3_conjMu
dCVar5(:,12)= transpose(d2G12_3(3,:)); %dG12_3_vnormv
dCVar5(:,13)= transpose(d2G23_1(1,:)); %dG23_1_Mu
dCVar5(:,14)= transpose(d2G23_1(2,:)); %dG23_1_conjMu
dCVar5(:,15)= transpose(d2G23_1(3,:)); %dG23_1_vnormv
dCVar5(:,16)= transpose(d2G31_2(1,:)); %dG31_2_Mu
dCVar5(:,17)= transpose(d2G31_2(2,:)); %dG31_2_conjMu
dCVar5(:,18)= transpose(d2G31_2(3,:)); %dG31_2_vnormv
dCVar5(:,19)= transpose(d2G11_2(1,:)); %dG11_2_Mu
dCVar5(:,20)= transpose(d2G11_2(2,:)); %dG11_2_conjMu
dCVar5(:,21)= transpose(d2G11_2(3,:)); %dG11_2_vnormv
dCVar5(:,22)= transpose(d2G11_3(1,:)); %dG11_3_Mu
dCVar5(:,23)= transpose(d2G11_3(2,:)); %dG11_3_conjMu
dCVar5(:,24)= transpose(d2G11_3(3,:)); %dG11_3_vnormv
dCVar5(:,25)= transpose(d2G22_1(1,:)); %dG22_1_Mu
dCVar5(:,26)= transpose(d2G22_1(2,:)); %dG22_1_conjMu
dCVar5(:,27)= transpose(d2G22_1(3,:)); %dG22_1_vnormv
dCVar5(:,28)= transpose(d2G22_3(1,:)); %dG22_3_Mu
dCVar5(:,29)= transpose(d2G22_3(2,:)); %dG22_3_conjMu
dCVar5(:,30)= transpose(d2G22_3(3,:)); %dG22_3_vnormv
dCVar5(:,31)= transpose(d2G33_1(1,:)); %dG33_1_Mu
dCVar5(:,32)= transpose(d2G33_1(2,:)); %dG33_1_conjMu
dCVar5(:,33)= transpose(d2G33_1(3,:)); %dG33_1_vnormv
dCVar5(:,34)= transpose(d2G33_2(1,:)); %dG33_2_Mu
dCVar5(:,35)= transpose(d2G33_2(2,:)); %dG33_2_conjMu
dCVar5(:,36)= transpose(d2G33_2(3,:)); %dG33_2_vnormv

dCVar5(:,37) = [0; 0; 0; 1; 0]; %u
dCVar5(:,38) = [dw_Mu; dw_conjMu; dw_vnormv; dw_u; 0]; %w
dCVar5(:,39) = transpose(d2w(1,:)); %dw_Mu
dCVar5(:,40) = transpose(d2w(2,:)); %dw_conjMu
dCVar5(:,41) = transpose(d2w(3,:)); %dw_vnormv
dCVar5(:,42) = transpose(d2w(4,:)); %dw_u
dCVar5(:,43) = transpose(d3w_u(1,:)); %d2w_uMu
dCVar5(:,44) = transpose(d3w_u(2,:)); %d2w_uconjMu
dCVar5(:,45) = transpose(d3w_u(3,:)); %d2w_uvnormv
dCVar5(:,46) = transpose(d3w_u(4,:)); %d2w_uu
dCVar5(:,47) = transpose(d4w_uu(1,:)); %d3w_uuMu
dCVar5(:,48) = transpose(d4w_uu(2,:)); %d3w_uuconjMu
dCVar5(:,49) = transpose(d4w_uu(3,:)); %d3w_uuvnormv
dCVar5(:,50) = transpose(d4w_uu(4,:)); %d3w_uuu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%