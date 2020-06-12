% variable_NatTwoR_W725_CR_rmW.m
Gset = [Gijk, dGijkE];
d2Gset = [d2G12_3_by11, d2G12_3_by12, d2G12_3_by13, d2G12_3_by22,...
    d2G12_3_by23, d2G12_3_by33, ...
    d2G23_1_by11, d2G23_1_by12, d2G23_1_by13, d2G23_1_by22,...
    d2G23_1_by23, d2G23_1_by33, ...
    d2G31_2_by11, d2G31_2_by12, d2G31_2_by13, d2G31_2_by22,...
    d2G31_2_by23, d2G31_2_by33, ...
    d2G11_2_by11, d2G11_2_by12, d2G11_2_by13, d2G11_2_by22,...
    d2G11_2_by23, d2G11_2_by33, ...
    d2G11_3_by11, d2G11_3_by12, d2G11_3_by13, d2G11_3_by22,...
    d2G11_3_by23, d2G11_3_by33, ...
    d2G22_1_by11, d2G22_1_by12, d2G22_1_by13, d2G22_1_by22,...
    d2G22_1_by23, d2G22_1_by33, ...
    d2G22_3_by11, d2G22_3_by12, d2G22_3_by13, d2G22_3_by22,...
    d2G22_3_by23, d2G22_3_by33, ...
    d2G33_1_by11, d2G33_1_by12, d2G33_1_by13, d2G33_1_by22,...
    d2G33_1_by23, d2G33_1_by33, ...
    d2G33_2_by11, d2G33_2_by12, d2G33_2_by13, d2G33_2_by22,...
    d2G33_2_by23, d2G33_2_by33];
% % % %
assumeAlso(Gset,'real');
assumeAlso(d2Gset,'real');
% df_Gijk_CR_rmW(FUN, Gset,GijkDict);
GijkDict.G12_3 = [dG12_3_by1; dG12_3_by2; dG12_3_by3];
GijkDict.G23_1 = [dG23_1_by1; dG23_1_by2; dG23_1_by3];
GijkDict.G31_2 = [dG31_2_by1; dG31_2_by2; dG31_2_by3];
GijkDict.G11_2 = [dG11_2_by1; dG11_2_by2; dG11_2_by3];
GijkDict.G11_3 = [dG11_3_by1; dG11_3_by2; dG11_3_by3];
GijkDict.G22_1 = [dG22_1_by1; dG22_1_by2; dG22_1_by3];
GijkDict.G22_3 = [dG22_3_by1; dG22_3_by2; dG22_3_by3];
GijkDict.G33_1 = [dG33_1_by1; dG33_1_by2; dG33_1_by3];
GijkDict.G33_2 = [dG33_2_by1; dG33_2_by2; dG33_2_by3];
GijkDict.dG12_3_by1 = [d2G12_3_by11; d2G12_3_by12; d2G12_3_by13];
GijkDict.dG12_3_by2 = [d2G12_3_by21; d2G12_3_by22; d2G12_3_by23];
GijkDict.dG12_3_by3 = [d2G12_3_by31; d2G12_3_by32; d2G12_3_by33];
GijkDict.dG23_1_by1 = [d2G23_1_by11; d2G23_1_by12; d2G23_1_by13];
GijkDict.dG23_1_by2 = [d2G23_1_by21; d2G23_1_by22; d2G23_1_by23];
GijkDict.dG23_1_by3 = [d2G23_1_by31; d2G23_1_by32; d2G23_1_by33];
GijkDict.dG31_2_by1 = [d2G31_2_by11; d2G31_2_by12; d2G31_2_by13];
GijkDict.dG31_2_by2 = [d2G31_2_by21; d2G31_2_by22; d2G31_2_by23];
GijkDict.dG31_2_by3 = [d2G31_2_by31; d2G31_2_by32; d2G31_2_by33];
GijkDict.dG11_2_by1 = [d2G11_2_by11; d2G11_2_by12; d2G11_2_by13];
GijkDict.dG11_2_by2 = [d2G11_2_by21; d2G11_2_by22; d2G11_2_by23];
GijkDict.dG11_2_by3 = [d2G11_2_by31; d2G11_2_by32; d2G11_2_by33];
GijkDict.dG11_3_by1 = [d2G11_3_by11; d2G11_3_by12; d2G11_3_by13];
GijkDict.dG11_3_by2 = [d2G11_3_by21; d2G11_3_by22; d2G11_3_by23];
GijkDict.dG11_3_by3 = [d2G11_3_by31; d2G11_3_by32; d2G11_3_by33];
GijkDict.dG22_1_by1 = [d2G22_1_by11; d2G22_1_by12; d2G22_1_by13];
GijkDict.dG22_1_by2 = [d2G22_1_by21; d2G22_1_by22; d2G22_1_by23];
GijkDict.dG22_1_by3 = [d2G22_1_by31; d2G22_1_by32; d2G22_1_by33];
GijkDict.dG22_3_by1 = [d2G22_3_by11; d2G22_3_by12; d2G22_3_by13];
GijkDict.dG22_3_by2 = [d2G22_3_by21; d2G22_3_by22; d2G22_3_by23];
GijkDict.dG22_3_by3 = [d2G22_3_by31; d2G22_3_by32; d2G22_3_by33];
GijkDict.dG33_1_by1 = [d2G33_1_by11; d2G33_1_by12; d2G33_1_by13];
GijkDict.dG33_1_by2 = [d2G33_1_by21; d2G33_1_by22; d2G33_1_by23];
GijkDict.dG33_1_by3 = [d2G33_1_by31; d2G33_1_by32; d2G33_1_by33];
GijkDict.dG33_2_by1 = [d2G33_2_by11; d2G33_2_by12; d2G33_2_by13];
GijkDict.dG33_2_by2 = [d2G33_2_by21; d2G33_2_by22; d2G33_2_by23];
GijkDict.dG33_2_by3 = [d2G33_2_by31; d2G33_2_by32; d2G33_2_by33];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct G in a 3x3x3 matrix; Gm(i,j,k) = G_{ij}^k
Gm = sym('Gm',[3 3 3]);
Gm(1,1,1) = 0;
Gm(1,1,2) = G11_2;
Gm(1,1,3) = G11_3;
Gm(1,2,1) = -G11_2;
Gm(1,2,2) = 0;
Gm(1,2,3) = G12_3;
Gm(1,3,1) = -G11_3;
Gm(1,3,2) = -G12_3;
Gm(1,3,3) = 0;
Gm(2,1,1) = 0;
Gm(2,1,2) = -G22_1;
Gm(2,1,3) = -G23_1;
Gm(2,2,1) = G22_1;
Gm(2,2,2) = 0;
Gm(2,2,3) = G22_3;
Gm(2,3,1) = G23_1;
Gm(2,3,2) = -G22_3;
Gm(2,3,3) = 0;
Gm(3,1,1) = 0;
Gm(3,1,2) = G31_2;
Gm(3,1,3) = -G33_1;
Gm(3,2,1) = -G31_2;
Gm(3,2,2) = 0;
Gm(3,2,3) = -G33_2;
Gm(3,3,1) = G33_1;
Gm(3,3,2) = G33_2;
Gm(3,3,3) = 0;
% rm: Riemannian curvature tensor on M
rm1212 = -dG22_1_by1 -dG11_2_by2 +G23_1*G12_3 +G11_3*G22_3...
    -G23_1*G31_2 -G12_3*G31_2 +G22_1*G22_1 + G11_2*G11_2; %ok
rm1213 = -dG23_1_by1 -dG11_3_by2 -G22_1*G12_3 -G11_2*G22_3...
    +G22_1*G23_1 +G23_1*G33_1 + G11_2*G11_3 +G12_3*G33_1; %ok
rm1223 = dG31_2_by2 +dG22_1_by3 +G33_1*G22_3 +G23_1*G33_2...
    -G31_2*G11_2 -G33_2*G31_2 -G23_1*G11_2 -G22_3*G22_1; %ok
rm1313 = -dG33_1_by1 -dG11_3_by3 +G31_2*G12_3 +G11_2*G33_2...
    -G31_2*G23_1 -G12_3*G23_1 + G33_1*G33_1 +G11_3*G11_3; %ok
rm1323 = -dG33_2_by1 -dG12_3_by3 -G31_2*G11_3 -G11_2*G33_1...
    +G31_2*G22_3 +G33_1*G33_2 +G11_3*G12_3 +G12_3*G22_3;
rm2323 = -dG33_2_by2 -dG22_3_by3 +G31_2*G23_1 + G22_1*G33_1...
    -G31_2*G12_3 -G23_1*G12_3 +G33_2*G33_2 +G22_3*G22_3;
% ricci curvature
ricm = sym('ricM',[3,3]);
ricm(1,1) = -rm1212-rm1313;
ricm(1,2) = -rm1323;
ricm(1,3) = rm1223;
ricm(2,1) = ricm(1,2);
ricm(2,2) = -rm1212-rm2323;
ricm(2,3) = -rm1213;
ricm(3,1) = ricm(1,3);
ricm(3,2) = ricm(2,3);
ricm(3,3) = -rm1313 -rm2323;
% scalar curvature
scalar = -2*(rm1212 +rm1313 +rm2323);
% schouten tensor
P = sym('schouten',[3,3]);
for j=1:3
    for k=1:3
        P(j,k) = ricm(j,k)-myDelta(j,k)*(scalar/4);
    end
end
% P(1,1) = 1/2*(rm2323 -rm1212 -rm1313);
% P(1,2) = -rm1323;
% P(1,3) = rm1223;
% P(2,1) = P(1,2);
% P(2,2) = 1/2*(rm1313 -rm1212 -rm2323);
% P(2,3) = -rm1213;
% P(3,1) = P(1,3);
% P(3,2) = P(2,3);
% P(3,3) = 1/2*(rm1212 -rm1313 -rm2323);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covP(m,n,k) = (\nabla_{e_m}P)(e_n,e_k);
covP = sym('cotton',[3,3,3]);
dscalar = df_Gijk_CR_rmW(scalar,Gset,GijkDict);
for m=1:3
    for n=1:3
        for k=1:3
            ricTemp = ricm(n,k);
            dricTemp = df_Gijk_CR_rmW(ricTemp,Gset,GijkDict);
            part1 = dricTemp(m)-1/4*dscalar(m)*myDelta(n,k);
            part2 = 0;
            for ll=1:3
                part2 = part2 -Gm(m,n,ll)*P(ll,k)...
                    -Gm(m,k,ll)*P(ll,n);
            end
            covP(m,n,k) = part1 + part2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NC substitution
% nCijk = (\nabla_{e_i} C)(e_j,e_k)
syms covP111 covP112 covP113 covP122 covP123 covP133 real
syms covP211 covP212 covP213 covP222 covP223 covP233 real
syms covP311 covP312 covP313 covP322 covP323 covP333 real

% 15 replacement of cotton terms
 %d2G33_2_by11 OK
[term112, gVec112] = coeffs(covP(1,1,2),d2Gset);
d2G33_2_by11Sub = covP112 -d2G12_3_by13 -term112(3);

%d2G33_2_by22 OK
[term211, gVec211] = coeffs(covP(2,1,1),d2Gset);
d2G33_2_by22Sub = d2G11_2_by22 + d2G11_3_by23 + d2G22_1_by12...
    -d2G22_3_by23 +d2G33_1_by12 -2*covP211 + 2*term211(7); 

%d2G31_2_by12 OK
[term113, gVec113] = coeffs(covP(1,1,3),d2Gset);
d2G31_2_by12Sub = covP113 - d2G22_1_by13 - term113(3); 

%d2G33_2_by23 OK
[term311, gVec311] = coeffs(covP(3,1,1),d2Gset);
d2G33_2_by23Sub = d2G11_2_by23Sub +d2G11_3_by33 +d2G22_1_by13...
    -d2G22_3_by33 +d2G33_1_by13 -2*covP311 +2*term311(7); 

%d2G33_2_by12 updated
[term122, gVec122] = coeffs(covP(1,2,2),d2Gset);
d2G33_2_by12Sub = 2*covP122 -d2G11_2_by12 +d2G11_3_by13 - d2G22_1_by11...
    -d2G22_3_by13Sub +d2G33_1_by11 -2*term122(7); 

%d2G12_3_by23 OK
[term212, gVec212] = coeffs(covP(2,1,2),d2Gset);
d2G12_3_by23Sub = covP212 - d2G33_2_by12 - term212(3); 

%d2G33_1_by11 OK
[term133, gVec133] =coeffs(covP(1,3,3),d2Gset);
d2G33_1_by11Sub = covP133 +d2G11_2_by12 -d2G11_3_by13 +d2G22_1_by11...
    -covP122 -term133(7) +term122(7);

%d2G31_2_by23 OK
[term313, gVec313] =coeffs(covP(3,1,3),d2Gset);
d2G31_2_by23Sub = covP313 -d2G22_1_by33 - term313(3); 

%d2G23_1_by12 OK
[term223, gVec223] =coeffs(covP(2,2,3),d2Gset);
d2G23_1_by12Sub = covP223 - d2G11_3_by22 - term223(3);

%d2G22_1_by13 OK
[term322, gVec322] =coeffs(covP(3,2,2),d2Gset);
d2G22_1_by13Sub = covP322+covP311-term311(7)-term322(7)-d2G11_2_by23Sub;

%d2G11_3_by23 OK
[term323, gVec323] =coeffs(covP(3,2,3),d2Gset);
d2G11_3_by23Sub = covP323 -d2G23_1_by13Sub -term323(3); 

%d2G33_1_by12 OK
[term233, gVec233] =coeffs(covP(2,3,3),d2Gset);
d2G33_1_by12Sub = covP233 +covP211 -term233(7) -term211(7) -d2G11_3_by23; 

% d2G23_1_by11 OK
[term123, gVec123] =coeffs(covP(1,2,3),d2Gset);
d2G23_1_by11Sub = covP123 - d2G11_3_by12 - term123(3); 

%d2G31_2_by22 OK
[term213, gVec213] =coeffs(covP(2,1,3),d2Gset);
d2G31_2_by22Sub = covP213 - d2G22_1_by23 - term213(3); 

%d2G33_2_by13 OK
[term312, gVec312] =coeffs(covP(3,1,2),d2Gset);
d2G33_2_by13Sub = covP312 -d2G12_3_by33 - term312(3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covPVariableSet1 = [d2G33_1_by11, d2G31_2_by23, d2G23_1_by12,...
    d2G22_1_by13, d2G11_3_by23, d2G33_2_by11, d2G23_1_by11,...
    d2G31_2_by22, d2G33_2_by13];
covPVariableSet2 = [d2G31_2_by12, d2G33_2_by23,...
    d2G33_2_by12, d2G33_1_by12];
covPVariableSet3 = [d2G12_3_by23, d2G33_2_by22];

covPSubSet1 = [d2G33_1_by11Sub, d2G31_2_by23Sub, d2G23_1_by12Sub,...
    d2G22_1_by13Sub, d2G11_3_by23Sub, d2G33_2_by11Sub, d2G23_1_by11Sub,...
    d2G31_2_by22Sub, d2G33_2_by13Sub];
covPSubSet2 = [d2G31_2_by12Sub, d2G33_2_by23Sub,...
    d2G33_2_by12Sub, d2G33_1_by12Sub];
covPSubSet3 = [d2G12_3_by23Sub, d2G33_2_by22Sub];

covPSet = [covP112, covP113, covP122, covP123, covP133,...
    covP211, covP212, covP213, covP223, covP233,...
    covP311, covP312, covP313, covP322, covP323];

totalSet = [Gset, d2Gset, covPSet];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




