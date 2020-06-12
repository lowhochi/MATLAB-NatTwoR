% variable_NatTwoR_Weyl7Part2_CR_rmW.m
syms G12_3 G23_1 G31_2 G11_2 G11_3 G22_1 G22_3 G33_1 G33_2 real
syms dG12_3_Mu dG12_3_conjMu dG12_3_vnormv
syms dG23_1_Mu dG23_1_conjMu dG23_1_vnormv
syms dG31_2_Mu dG31_2_conjMu dG31_2_vnormv 
syms dG11_2_Mu dG11_2_conjMu dG11_2_vnormv 
syms dG11_3_Mu dG11_3_conjMu dG11_3_vnormv 
syms dG22_1_Mu dG22_1_conjMu dG22_1_vnormv 
syms dG22_3_Mu dG22_3_conjMu dG22_3_vnormv
syms dG33_1_Mu dG33_1_conjMu dG33_1_vnormv
syms dG33_2_Mu dG33_2_conjMu dG33_2_vnormv
% update lieMu;
lieMu{1,2} = [aMSub; -conj(aMSub); 2*i*h11SubTwo];
lieMu{1,3} = [aV; bVSub; 2*T4Sub]; 
lieMu{2,1} = [-aMSub; conj(aMSub); -2*i*h11SubTwo];
lieMu{2,3} = [conj(bVSub); conj(aV); 2*conj(T4Sub)];
lieMu{3,1} = [-aV; -bVSub; -2*T4Sub];
lieMu{3,2} = [-conj(bVSub); -conj(aV); -2*conj(T4Sub)];

% The second derivatives d2G12_3
syms d2G12_3_MuMu d2G12_3_MuconjMu d2G12_3_Muvnormv
syms d2G12_3_conjMuconjMu d2G12_3_conjMuvnormv d2G12_3_vnormvvnormv
dG12_3Row = [dG12_3_Mu, dG12_3_conjMu, dG12_3_vnormv];
d2G12_3_conjMuMu = d2G12_3_MuconjMu + dG12_3Row*lieMu{1,2};
d2G12_3_vnormvMu = d2G12_3_Muvnormv + dG12_3Row*lieMu{1,3};
d2G12_3_vnormvconjMu = d2G12_3_conjMuvnormv + dG12_3Row*lieMu{2,3}; 
d2G12_3_Muu = dG12_3Row*lieMu{4,1};
d2G12_3_conjMuconju = conj(d2G12_3_Muu);
d2G12_3_vnormvu = dG12_3Row*lieMu{4,3};
d2G12_3_vnormvconju = dG12_3Row*lieMu{5,3};
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

syms dG12_3_by1 dG12_3_by2 dG12_3_by3 real
syms dG23_1_by1 dG23_1_by2 dG23_1_by3 real 
syms dG31_2_by1 dG31_2_by2 dG31_2_by3 real 
syms dG11_2_by1 dG11_2_by2 dG11_2_by3 real 
syms dG11_3_by1 dG11_3_by2 dG11_3_by3 real 
syms dG22_1_by1 dG22_1_by2 dG22_1_by3 real 
syms dG22_3_by1 dG22_3_by2 dG22_3_by3 real 
syms dG33_1_by1 dG33_1_by2 dG33_1_by3 real 
syms dG33_2_by1 dG33_2_by2 dG33_2_by3 real 

muVec = [mu1;mu2;mu3];
conjMuVec = [conj(mu1);conj(mu2);conj(mu3)];
vnormvVec = [v1normv;v2normv;v3normv];

Gijk = [G11_2, G11_3, G12_3, G22_1, G22_3, G23_1, G31_2, G33_1, G33_2];
dGijkE = [dG12_3_by1, dG12_3_by2, dG12_3_by3,...
    dG23_1_by1, dG23_1_by2, dG23_1_by3,...
    dG31_2_by1, dG31_2_by2, dG31_2_by3,...
    dG11_2_by1, dG11_2_by2, dG11_2_by3,...
    dG11_3_by1, dG11_3_by2, dG11_3_by3,...
    dG22_1_by1, dG22_1_by2, dG22_1_by3,...
    dG22_3_by1, dG22_3_by2, dG22_3_by3,...
    dG33_1_by1, dG33_1_by2, dG33_1_by3,...
    dG33_2_by1, dG33_2_by2, dG33_2_by3]; 

% construct lieE{i,j} = [e_i, e_j] under the basis {e1,e2,e3}
lieE = cell(3,3);
lieE{1,1} = zeros(3,1);
lieE{1,2} = [-G11_2; G22_1; G12_3+G23_1];
lieE{1,3} = [-G11_3; -G12_3-G31_2; G33_1];
lieE{2,1} = -lieE{1,2};
lieE{2,2} = zeros(3,1);
lieE{2,3} = [G23_1+G31_2; -G22_3; G33_2];
lieE{3,1} = -lieE{1,3};
lieE{3,2} = -lieE{2,3};
lieE{3,3} = zeros(3,1);

syms d2G12_3_by11 d2G12_3_by12 d2G12_3_by13 real  
syms d2G12_3_by22 d2G12_3_by23 d2G12_3_by33 real 
dG12_3ERow = [dG12_3_by1, dG12_3_by2, dG12_3_by3];
d2G12_3_by21 = d2G12_3_by12 + dG12_3ERow*lieE{1,2};
d2G12_3_by31 = d2G12_3_by13 + dG12_3ERow*lieE{1,3};
d2G12_3_by32 = d2G12_3_by23 + dG12_3ERow*lieE{2,3};

d2G12_3 = [d2G12_3_by11, d2G12_3_by12, d2G12_3_by13;
    d2G12_3_by21, d2G12_3_by22, d2G12_3_by23; 
    d2G12_3_by31, d2G12_3_by32, d2G12_3_by33];

syms d2G23_1_by11 d2G23_1_by12 d2G23_1_by13 real 
syms d2G23_1_by22 d2G23_1_by23 d2G23_1_by33 real 
dG23_1ERow = [dG23_1_by1, dG23_1_by2, dG23_1_by3];
d2G23_1_by21 = d2G23_1_by12 + dG23_1ERow*lieE{1,2};
d2G23_1_by31 = d2G23_1_by13 + dG23_1ERow*lieE{1,3};
d2G23_1_by32 = d2G23_1_by23 + dG23_1ERow*lieE{2,3};

d2G23_1 = [d2G23_1_by11, d2G23_1_by12, d2G23_1_by13;
    d2G23_1_by21, d2G23_1_by22, d2G23_1_by23; 
    d2G23_1_by31, d2G23_1_by32, d2G23_1_by33];

syms d2G31_2_by11 d2G31_2_by12 d2G31_2_by13 real  
syms d2G31_2_by22 d2G31_2_by23 d2G31_2_by33 real  
dG31_2ERow = [dG31_2_by1, dG31_2_by2, dG31_2_by3];
d2G31_2_by21 = d2G31_2_by12 + dG31_2ERow*lieE{1,2};
d2G31_2_by31 = d2G31_2_by13 + dG31_2ERow*lieE{1,3};
d2G31_2_by32 = d2G31_2_by23 + dG31_2ERow*lieE{2,3};

d2G31_2 = [d2G31_2_by11, d2G31_2_by12, d2G31_2_by13;
    d2G31_2_by21, d2G31_2_by22, d2G31_2_by23; 
    d2G31_2_by31, d2G31_2_by32, d2G31_2_by33];

syms d2G11_2_by11 d2G11_2_by12 d2G11_2_by13 real 
syms d2G11_2_by22 d2G11_2_by23 d2G11_2_by33 real
dG11_2ERow = [dG11_2_by1, dG11_2_by2, dG11_2_by3];
d2G11_2_by21 = d2G11_2_by12 + dG11_2ERow*lieE{1,2};
d2G11_2_by31 = d2G11_2_by13 + dG11_2ERow*lieE{1,3};
d2G11_2_by32 = d2G11_2_by23 + dG11_2ERow*lieE{2,3};

d2G11_2 = [d2G11_2_by11, d2G11_2_by12, d2G11_2_by13;
    d2G11_2_by21, d2G11_2_by22, d2G11_2_by23; 
    d2G11_2_by31, d2G11_2_by32, d2G11_2_by33];

syms d2G11_3_by11 d2G11_3_by12 d2G11_3_by13 real 
syms d2G11_3_by22 d2G11_3_by23 d2G11_3_by33 real  
dG11_3ERow = [dG11_3_by1, dG11_3_by2, dG11_3_by3];
d2G11_3_by21 = d2G11_3_by12 + dG11_3ERow*lieE{1,2};
d2G11_3_by31 = d2G11_3_by13 + dG11_3ERow*lieE{1,3};
d2G11_3_by32 = d2G11_3_by23 + dG11_3ERow*lieE{2,3};

d2G11_3 = [d2G11_3_by11, d2G11_3_by12, d2G11_3_by13;
    d2G11_3_by21, d2G11_3_by22, d2G11_3_by23; 
    d2G11_3_by31, d2G11_3_by32, d2G11_3_by33];

syms d2G22_1_by11 d2G22_1_by12 d2G22_1_by13 real 
syms d2G22_1_by22 d2G22_1_by23 d2G22_1_by33 real  
dG22_1ERow = [dG22_1_by1, dG22_1_by2, dG22_1_by3];
d2G22_1_by21 = d2G22_1_by12 + dG22_1ERow*lieE{1,2};
d2G22_1_by31 = d2G22_1_by13 + dG22_1ERow*lieE{1,3};
d2G22_1_by32 = d2G22_1_by23 + dG22_1ERow*lieE{2,3};

d2G22_1 = [d2G22_1_by11, d2G22_1_by12, d2G22_1_by13;
    d2G22_1_by21, d2G22_1_by22, d2G22_1_by23; 
    d2G22_1_by31, d2G22_1_by32, d2G22_1_by33];

syms d2G22_3_by11 d2G22_3_by12 d2G22_3_by13 real 
syms d2G22_3_by22 d2G22_3_by23 d2G22_3_by33 real 
dG22_3ERow = [dG22_3_by1, dG22_3_by2, dG22_3_by3];
d2G22_3_by21 = d2G22_3_by12 + dG22_3ERow*lieE{1,2};
d2G22_3_by31 = d2G22_3_by13 + dG22_3ERow*lieE{1,3};
d2G22_3_by32 = d2G22_3_by23 + dG22_3ERow*lieE{2,3};

d2G22_3 = [d2G22_3_by11, d2G22_3_by12, d2G22_3_by13;
    d2G22_3_by21, d2G22_3_by22, d2G22_3_by23; 
    d2G22_3_by31, d2G22_3_by32, d2G22_3_by33];

syms d2G33_1_by11 d2G33_1_by12 d2G33_1_by13 real 
syms d2G33_1_by22 d2G33_1_by23 d2G33_1_by33 real  
dG33_1ERow = [dG33_1_by1, dG33_1_by2, dG33_1_by3];
d2G33_1_by21 = d2G33_1_by12 + dG33_1ERow*lieE{1,2};
d2G33_1_by31 = d2G33_1_by13 + dG33_1ERow*lieE{1,3};
d2G33_1_by32 = d2G33_1_by23 + dG33_1ERow*lieE{2,3};

d2G33_1 = [d2G33_1_by11, d2G33_1_by12, d2G33_1_by13;
    d2G33_1_by21, d2G33_1_by22, d2G33_1_by23; 
    d2G33_1_by31, d2G33_1_by32, d2G33_1_by33];

syms d2G33_2_by11 d2G33_2_by12 d2G33_2_by13 real 
syms d2G33_2_by22 d2G33_2_by23 d2G33_2_by33 real  
dG33_2ERow = [dG33_2_by1, dG33_2_by2, dG33_2_by3];
d2G33_2_by21 = d2G33_2_by12 + dG33_2ERow*lieE{1,2};
d2G33_2_by31 = d2G33_2_by13 + dG33_2ERow*lieE{1,3};
d2G33_2_by32 = d2G33_2_by23 + dG33_2ERow*lieE{2,3};

d2G33_2 = [d2G33_2_by11, d2G33_2_by12, d2G33_2_by13;
    d2G33_2_by21, d2G33_2_by22, d2G33_2_by23; 
    d2G33_2_by31, d2G33_2_by32, d2G33_2_by33];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dG12_3_MuSub = transpose(muVec)*[dG12_3_by1; dG12_3_by2; dG12_3_by3];
dG12_3_conjMuSub = transpose(conjMuVec)*[dG12_3_by1; dG12_3_by2; dG12_3_by3];
dG12_3_vnormvSub = transpose(vnormvVec)*[dG12_3_by1; dG12_3_by2; dG12_3_by3];
dG23_1_MuSub = transpose(muVec)*[dG23_1_by1; dG23_1_by2; dG23_1_by3];
dG23_1_conjMuSub = transpose(conjMuVec)*[dG23_1_by1; dG23_1_by2; dG23_1_by3];
dG23_1_vnormvSub = transpose(vnormvVec)*[dG23_1_by1; dG23_1_by2; dG23_1_by3];
dG31_2_MuSub = transpose(muVec)*[dG31_2_by1; dG31_2_by2; dG31_2_by3];
dG31_2_conjMuSub = transpose(conjMuVec)*[dG31_2_by1; dG31_2_by2; dG31_2_by3];
dG31_2_vnormvSub = transpose(vnormvVec)*[dG31_2_by1; dG31_2_by2; dG31_2_by3];
dG11_2_MuSub = transpose(muVec)*[dG11_2_by1; dG11_2_by2; dG11_2_by3];
dG11_2_conjMuSub = transpose(conjMuVec)*[dG11_2_by1; dG11_2_by2; dG11_2_by3];
dG11_2_vnormvSub = transpose(vnormvVec)*[dG11_2_by1; dG11_2_by2; dG11_2_by3];
dG11_3_MuSub = transpose(muVec)*[dG11_3_by1; dG11_3_by2; dG11_3_by3];
dG11_3_conjMuSub = transpose(conjMuVec)*[dG11_3_by1; dG11_3_by2; dG11_3_by3];
dG11_3_vnormvSub = transpose(vnormvVec)*[dG11_3_by1; dG11_3_by2; dG11_3_by3];
dG22_1_MuSub = transpose(muVec)*[dG22_1_by1; dG22_1_by2; dG22_1_by3];
dG22_1_conjMuSub = transpose(conjMuVec)*[dG22_1_by1; dG22_1_by2; dG22_1_by3];
dG22_1_vnormvSub = transpose(vnormvVec)*[dG22_1_by1; dG22_1_by2; dG22_1_by3];
dG22_3_MuSub = transpose(muVec)*[dG22_3_by1; dG22_3_by2; dG22_3_by3];
dG22_3_conjMuSub = transpose(conjMuVec)*[dG22_3_by1; dG22_3_by2; dG22_3_by3];
dG22_3_vnormvSub = transpose(vnormvVec)*[dG22_3_by1; dG22_3_by2; dG22_3_by3];
dG33_1_MuSub = transpose(muVec)*[dG33_1_by1; dG33_1_by2; dG33_1_by3];
dG33_1_conjMuSub = transpose(conjMuVec)*[dG33_1_by1; dG33_1_by2; dG33_1_by3];
dG33_1_vnormvSub = transpose(vnormvVec)*[dG33_1_by1; dG33_1_by2; dG33_1_by3];
dG33_2_MuSub = transpose(muVec)*[dG33_2_by1; dG33_2_by2; dG33_2_by3];
dG33_2_conjMuSub = transpose(conjMuVec)*[dG33_2_by1; dG33_2_by2; dG33_2_by3];
dG33_2_vnormvSub = transpose(vnormvVec)*[dG33_2_by1; dG33_2_by2; dG33_2_by3];

variableSetMu = [dG12_3_Mu, dG23_1_Mu, dG31_2_Mu,...
    dG11_2_Mu, dG11_3_Mu, dG22_1_Mu, dG22_3_Mu, dG33_1_Mu, dG33_2_Mu];
variableSetConjMu = [dG12_3_conjMu, dG23_1_conjMu, dG31_2_conjMu,...
    dG11_2_conjMu, dG11_3_conjMu, dG22_1_conjMu, dG22_3_conjMu,...
    dG33_1_conjMu, dG33_2_conjMu];
variableSetVnormv = [dG12_3_vnormv, dG23_1_vnormv, dG31_2_vnormv,...
    dG11_2_vnormv, dG11_3_vnormv, dG22_1_vnormv, dG22_3_vnormv,...
    dG33_1_vnormv, dG33_2_vnormv];

subSetMu = [dG12_3_MuSub, dG23_1_MuSub, dG31_2_MuSub,...
    dG11_2_MuSub, dG11_3_MuSub, dG22_1_MuSub, dG22_3_MuSub,...
    dG33_1_MuSub, dG33_2_MuSub];
subSetConjMu = [dG12_3_conjMuSub, dG23_1_conjMuSub, dG31_2_conjMuSub,...
    dG11_2_conjMuSub, dG11_3_conjMuSub, dG22_1_conjMuSub, dG22_3_conjMuSub,...
    dG33_1_conjMuSub, dG33_2_conjMuSub];
subSetVnormv = [dG12_3_vnormvSub, dG23_1_vnormvSub, dG31_2_vnormvSub,...
    dG11_2_vnormvSub, dG11_3_vnormvSub, dG22_1_vnormvSub, dG22_3_vnormvSub,...
    dG33_1_vnormvSub, dG33_2_vnormvSub];

d2G12_3N = [d2G12_3_MuMu, d2G12_3_MuconjMu, d2G12_3_Muvnormv,...
    d2G12_3_conjMuconjMu, d2G12_3_conjMuvnormv, d2G12_3_vnormvvnormv];
d2G12_3_MuMuSub = transpose(muVec)*d2G12_3*muVec;
d2G12_3_MuconjMuSub = transpose(muVec)*d2G12_3*conjMuVec;
d2G12_3_MuvnormvSub = transpose(muVec)*d2G12_3*vnormvVec;
d2G12_3_conjMuconjMuSub = transpose(conjMuVec)*d2G12_3*conjMuVec;
d2G12_3_conjMuvnormvSub = transpose(conjMuVec)*d2G12_3*vnormvVec;
d2G12_3_vnormvvnormvSub = transpose(vnormvVec)*d2G12_3*vnormvVec;
d2G12_3Vec = [d2G12_3_MuMuSub, d2G12_3_MuconjMuSub, d2G12_3_MuvnormvSub,...
    d2G12_3_conjMuconjMuSub, d2G12_3_conjMuvnormvSub, d2G12_3_vnormvvnormvSub];
% %
d2G23_1N = [d2G23_1_MuMu, d2G23_1_MuconjMu, d2G23_1_Muvnormv,...
    d2G23_1_conjMuconjMu, d2G23_1_conjMuvnormv, d2G23_1_vnormvvnormv];
d2G23_1_MuMuSub = transpose(muVec)*d2G23_1*muVec;
d2G23_1_MuconjMuSub = transpose(muVec)*d2G23_1*conjMuVec;
d2G23_1_MuvnormvSub = transpose(muVec)*d2G23_1*vnormvVec;
d2G23_1_conjMuconjMuSub = transpose(conjMuVec)*d2G23_1*conjMuVec;
d2G23_1_conjMuvnormvSub = transpose(conjMuVec)*d2G23_1*vnormvVec;
d2G23_1_vnormvvnormvSub = transpose(vnormvVec)*d2G23_1*vnormvVec;
d2G23_1Vec = [d2G23_1_MuMuSub, d2G23_1_MuconjMuSub, d2G23_1_MuvnormvSub,...
    d2G23_1_conjMuconjMuSub, d2G23_1_conjMuvnormvSub, d2G23_1_vnormvvnormvSub];
% %
d2G31_2N = [d2G31_2_MuMu, d2G31_2_MuconjMu, d2G31_2_Muvnormv,...
    d2G31_2_conjMuconjMu, d2G31_2_conjMuvnormv, d2G31_2_vnormvvnormv];
d2G31_2_MuMuSub = transpose(muVec)*d2G31_2*muVec;
d2G31_2_MuconjMuSub = transpose(muVec)*d2G31_2*conjMuVec;
d2G31_2_MuvnormvSub = transpose(muVec)*d2G31_2*vnormvVec;
d2G31_2_conjMuconjMuSub = transpose(conjMuVec)*d2G31_2*conjMuVec;
d2G31_2_conjMuvnormvSub = transpose(conjMuVec)*d2G31_2*vnormvVec;
d2G31_2_vnormvvnormvSub = transpose(vnormvVec)*d2G31_2*vnormvVec;
d2G31_2Vec = [d2G31_2_MuMuSub, d2G31_2_MuconjMuSub, d2G31_2_MuvnormvSub,...
    d2G31_2_conjMuconjMuSub, d2G31_2_conjMuvnormvSub, d2G31_2_vnormvvnormvSub];
% %
d2G11_2N = [d2G11_2_MuMu, d2G11_2_MuconjMu, d2G11_2_Muvnormv,...
    d2G11_2_conjMuconjMu, d2G11_2_conjMuvnormv, d2G11_2_vnormvvnormv];
d2G11_2_MuMuSub = transpose(muVec)*d2G11_2*muVec;
d2G11_2_MuconjMuSub = transpose(muVec)*d2G11_2*conjMuVec;
d2G11_2_MuvnormvSub = transpose(muVec)*d2G11_2*vnormvVec;
d2G11_2_conjMuconjMuSub = transpose(conjMuVec)*d2G11_2*conjMuVec;
d2G11_2_conjMuvnormvSub = transpose(conjMuVec)*d2G11_2*vnormvVec;
d2G11_2_vnormvvnormvSub = transpose(vnormvVec)*d2G11_2*vnormvVec;
d2G11_2Vec = [d2G11_2_MuMuSub, d2G11_2_MuconjMuSub, d2G11_2_MuvnormvSub,...
    d2G11_2_conjMuconjMuSub, d2G11_2_conjMuvnormvSub, d2G11_2_vnormvvnormvSub];
% %
d2G11_3N = [d2G11_3_MuMu, d2G11_3_MuconjMu, d2G11_3_Muvnormv,...
    d2G11_3_conjMuconjMu, d2G11_3_conjMuvnormv, d2G11_3_vnormvvnormv];
d2G11_3_MuMuSub = transpose(muVec)*d2G11_3*muVec;
d2G11_3_MuconjMuSub = transpose(muVec)*d2G11_3*conjMuVec;
d2G11_3_MuvnormvSub = transpose(muVec)*d2G11_3*vnormvVec;
d2G11_3_conjMuconjMuSub = transpose(conjMuVec)*d2G11_3*conjMuVec;
d2G11_3_conjMuvnormvSub = transpose(conjMuVec)*d2G11_3*vnormvVec;
d2G11_3_vnormvvnormvSub = transpose(vnormvVec)*d2G11_3*vnormvVec;
d2G11_3Vec = [d2G11_3_MuMuSub, d2G11_3_MuconjMuSub, d2G11_3_MuvnormvSub,...
    d2G11_3_conjMuconjMuSub, d2G11_3_conjMuvnormvSub, d2G11_3_vnormvvnormvSub];
% %
d2G22_1N = [d2G22_1_MuMu, d2G22_1_MuconjMu, d2G22_1_Muvnormv,...
    d2G22_1_conjMuconjMu, d2G22_1_conjMuvnormv, d2G22_1_vnormvvnormv];
d2G22_1_MuMuSub = transpose(muVec)*d2G22_1*muVec;
d2G22_1_MuconjMuSub = transpose(muVec)*d2G22_1*conjMuVec;
d2G22_1_MuvnormvSub = transpose(muVec)*d2G22_1*vnormvVec;
d2G22_1_conjMuconjMuSub = transpose(conjMuVec)*d2G22_1*conjMuVec;
d2G22_1_conjMuvnormvSub = transpose(conjMuVec)*d2G22_1*vnormvVec;
d2G22_1_vnormvvnormvSub = transpose(vnormvVec)*d2G22_1*vnormvVec;
d2G22_1Vec = [d2G22_1_MuMuSub, d2G22_1_MuconjMuSub, d2G22_1_MuvnormvSub,...
    d2G22_1_conjMuconjMuSub, d2G22_1_conjMuvnormvSub, d2G22_1_vnormvvnormvSub];
% %
d2G22_3N = [d2G22_3_MuMu, d2G22_3_MuconjMu, d2G22_3_Muvnormv,...
    d2G22_3_conjMuconjMu, d2G22_3_conjMuvnormv, d2G22_3_vnormvvnormv];
d2G22_3_MuMuSub = transpose(muVec)*d2G22_3*muVec;
d2G22_3_MuconjMuSub = transpose(muVec)*d2G22_3*conjMuVec;
d2G22_3_MuvnormvSub = transpose(muVec)*d2G22_3*vnormvVec;
d2G22_3_conjMuconjMuSub = transpose(conjMuVec)*d2G22_3*conjMuVec;
d2G22_3_conjMuvnormvSub = transpose(conjMuVec)*d2G22_3*vnormvVec;
d2G22_3_vnormvvnormvSub = transpose(vnormvVec)*d2G22_3*vnormvVec;
d2G22_3Vec = [d2G22_3_MuMuSub, d2G22_3_MuconjMuSub, d2G22_3_MuvnormvSub,...
    d2G22_3_conjMuconjMuSub, d2G22_3_conjMuvnormvSub, d2G22_3_vnormvvnormvSub];
% %
d2G33_1N = [d2G33_1_MuMu, d2G33_1_MuconjMu, d2G33_1_Muvnormv,...
    d2G33_1_conjMuconjMu, d2G33_1_conjMuvnormv, d2G33_1_vnormvvnormv];
d2G33_1_MuMuSub = transpose(muVec)*d2G33_1*muVec;
d2G33_1_MuconjMuSub = transpose(muVec)*d2G33_1*conjMuVec;
d2G33_1_MuvnormvSub = transpose(muVec)*d2G33_1*vnormvVec;
d2G33_1_conjMuconjMuSub = transpose(conjMuVec)*d2G33_1*conjMuVec;
d2G33_1_conjMuvnormvSub = transpose(conjMuVec)*d2G33_1*vnormvVec;
d2G33_1_vnormvvnormvSub = transpose(vnormvVec)*d2G33_1*vnormvVec;
d2G33_1Vec = [d2G33_1_MuMuSub, d2G33_1_MuconjMuSub, d2G33_1_MuvnormvSub,...
    d2G33_1_conjMuconjMuSub, d2G33_1_conjMuvnormvSub, d2G33_1_vnormvvnormvSub];
% %
d2G33_2N = [d2G33_2_MuMu, d2G33_2_MuconjMu, d2G33_2_Muvnormv,...
    d2G33_2_conjMuconjMu, d2G33_2_conjMuvnormv, d2G33_2_vnormvvnormv];
d2G33_2_MuMuSub = transpose(muVec)*d2G33_2*muVec;
d2G33_2_MuconjMuSub = transpose(muVec)*d2G33_2*conjMuVec;
d2G33_2_MuvnormvSub = transpose(muVec)*d2G33_2*vnormvVec;
d2G33_2_conjMuconjMuSub = transpose(conjMuVec)*d2G33_2*conjMuVec;
d2G33_2_conjMuvnormvSub = transpose(conjMuVec)*d2G33_2*vnormvVec;
d2G33_2_vnormvvnormvSub = transpose(vnormvVec)*d2G33_2*vnormvVec;
d2G33_2Vec = [d2G33_2_MuMuSub, d2G33_2_MuconjMuSub, d2G33_2_MuvnormvSub,...
    d2G33_2_conjMuconjMuSub, d2G33_2_conjMuvnormvSub, d2G33_2_vnormvvnormvSub];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bianchi Identities
dG11_2_by3Sub = dG23_1_by1 + dG11_3_by2 + dG31_2_by1...
    - G11_3*G33_2 + G11_2*G22_3 - (G22_1+G33_1)*(G31_2+G23_1); %ok

d2G11_2_by31Sub = d2G23_1_by11 +  d2G11_3_by21 + d2G31_2_by11...
    - dG11_3_by1*G33_2 - G11_3*dG33_2_by1...
    + dG11_2_by1*G22_3 + G11_2*dG22_3_by1...
    - (dG22_1_by1+dG33_1_by1)*(G31_2+G23_1)...
    - (G22_1+G33_1)*(dG31_2_by1+dG23_1_by1);
d2G11_2_by13Sub = d2G11_2_by31Sub + G11_3*dG11_2_by1 - G33_1*dG11_2_by3...
    + dG11_2_by2*(G12_3 + G31_2); %ok

d2G11_2_by32Sub = d2G23_1_by12 +  d2G11_3_by22 + d2G31_2_by12...
    - dG11_3_by2*G33_2 - G11_3*dG33_2_by2...
    + dG11_2_by2*G22_3 + G11_2*dG22_3_by2...
    - (dG22_1_by2+dG33_1_by2)*(G31_2+G23_1)...
    - (G22_1+G33_1)*(dG31_2_by2+dG23_1_by2);
d2G11_2_by23Sub = d2G11_2_by32Sub + G22_3*dG11_2_by2 - G33_2*dG11_2_by3...
    - dG11_2_by1*(G23_1 + G31_2); %ok

d2G11_2_by33Sub = d2G23_1_by13 +  d2G11_3_by23 + d2G31_2_by13...
    - dG11_3_by3*G33_2 - G11_3*dG33_2_by3...
    + dG11_2_by3*G22_3 + G11_2*dG22_3_by3...
    - (dG22_1_by3+dG33_1_by3)*(G31_2+G23_1)...
    - (G22_1+G33_1)*(dG31_2_by3+dG23_1_by3); %ok
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dG22_3_by1Sub = dG12_3_by2 + dG31_2_by2 + dG22_1_by3...
    + G33_1*G22_3 - G22_1*G11_3 - (G11_2+G33_2)*(G31_2+G12_3); %ok

d2G22_3_by11Sub = d2G12_3_by21 + d2G31_2_by21 + d2G22_1_by31...
    + dG33_1_by1*G22_3 + G33_1*dG22_3_by1...
    - dG22_1_by1*G11_3 - G22_1*dG11_3_by1...
    - (dG11_2_by1+dG33_2_by1)*(G31_2+G12_3)...
    - (G11_2+G33_2)*(dG31_2_by1+dG12_3_by1); %ok

d2G22_3_by12Sub = d2G12_3_by22 + d2G31_2_by22 + d2G22_1_by32...
    + dG33_1_by2*G22_3 + G33_1*dG22_3_by2...
    - dG22_1_by2*G11_3 - G22_1*dG11_3_by2...
    - (dG11_2_by2+dG33_2_by2)*(G31_2+G12_3)...
    - (G11_2+G33_2)*(dG31_2_by2+dG12_3_by2); %ok

d2G22_3_by13Sub = d2G12_3_by23 + d2G31_2_by23 + d2G22_1_by33...
    + dG33_1_by3*G22_3 + G33_1*dG22_3_by3...
    - dG22_1_by3*G11_3 - G22_1*dG11_3_by3...
    - (dG11_2_by3+dG33_2_by3)*(G31_2+G12_3)...
    - (G11_2+G33_2)*(dG31_2_by3+dG12_3_by3); %ok
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dG23_1_by3Sub = dG33_1_by2 - dG33_2_by1 - dG12_3_by3...
    - G11_2*G33_1 + G22_1*G33_2 + (G11_3+G22_3)*(G12_3+G23_1); %ok

d2G23_1_by31Sub = d2G33_1_by21 - d2G33_2_by11 - d2G12_3_by31...
    - dG11_2_by1*G33_1 - G11_2*dG33_1_by1...
    + dG22_1_by1*G33_2 + G22_1*dG33_2_by1...
    + (dG11_3_by1+dG22_3_by1)*(G12_3+G23_1)...
    + (G11_3+G22_3)*(dG12_3_by1+dG23_1_by1);
d2G23_1_by13Sub = d2G23_1_by31Sub + G11_3*dG23_1_by1 - G33_1*dG23_1_by3...
    + dG23_1_by2*(G12_3 + G31_2); %ok

d2G23_1_by32Sub = d2G33_1_by22 - d2G33_2_by12 - d2G12_3_by32...
    - dG11_2_by2*G33_1 - G11_2*dG33_1_by2...
    + dG22_1_by2*G33_2 + G22_1*dG33_2_by2...
    + (dG11_3_by2+dG22_3_by2)*(G12_3+G23_1)...
    + (G11_3+G22_3)*(dG12_3_by2+dG23_1_by2);
d2G23_1_by23Sub = d2G23_1_by32Sub + G22_3*dG23_1_by2 - G33_2*dG23_1_by3...
    - dG23_1_by1*(G23_1 + G31_2); %ok

d2G23_1_by33Sub = d2G33_1_by23 - d2G33_2_by13 - d2G12_3_by33...
    - dG11_2_by3*G33_1 - G11_2*dG33_1_by3...
    + dG22_1_by3*G33_2 + G22_1*dG33_2_by3...
    + (dG11_3_by3+dG22_3_by3)*(G12_3+G23_1)...
    + (G11_3+G22_3)*(dG12_3_by3+dG23_1_by3); %ok

variableSetBianchi1 = [dG11_2_by3, dG22_3_by1, dG23_1_by3];
subSetBianchi1 = [dG11_2_by3Sub, dG22_3_by1Sub, dG23_1_by3Sub];

variableSetBianchi2 = [d2G11_2_by13, d2G11_2_by23, d2G11_2_by33,...
    d2G22_3_by11, d2G22_3_by12, d2G22_3_by13,...
    d2G23_1_by13, d2G23_1_by23, d2G23_1_by33];
subSetBianchi2 = [d2G11_2_by13Sub, d2G11_2_by23Sub, d2G11_2_by33Sub,...
    d2G22_3_by11Sub, d2G22_3_by12Sub, d2G22_3_by13Sub,...
    d2G23_1_by13Sub, d2G23_1_by23Sub, d2G23_1_by33Sub];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
