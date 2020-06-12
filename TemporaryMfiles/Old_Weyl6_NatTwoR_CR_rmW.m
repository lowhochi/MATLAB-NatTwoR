% Weyl6_NatTwoR_CR_rmW.m
% M is flat in the Weyl6 files.
load('DataWeylTwo_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myNum = length(symSetWeylTwo);
subSetFlat = sym('subSet',[1,myNum]);
wSet6 =[];
for j=1:myNum
    myVar = symSetWeylTwo(j);
    myChar = char(myVar);
    if strfind(myChar,'dw')==1
       subSetFlat(j)=myVar;
   elseif strfind(myChar,'d2w')==1
       subSetFlat(j)=myVar;
   elseif strfind(myChar,'d3w')==1
       subSetFlat(j)=myVar;
    elseif strfind(myChar,'d4w')==1
       subSetFlat(j)=myVar;
   elseif myVar==u
       subSetFlat(j)=myVar;
   elseif myVar==w
       subSetFlat(j)=myVar;
   elseif myVar==rho
       subSetFlat(j)=myVar; 
   else
       subSetFlat(j)=0;
       continue
   end
   wSet6 = union(wSet6,[myVar]);
end

WeylFlat = sym('WeylFlat',[120,5]);
for j=1:120
    WeylFlat(j,1) = WeylTwo(j,1);
    WeylFlat(j,2) = WeylTwo(j,2);
    WeylFlat(j,3) = WeylTwo(j,3);
    WeylFlat(j,4)= WeylTwo(j,4);
    temp = WeylTwo(j,5);
    temp = subs(temp, Y, 1+u*conj(u));
    temp = subs(temp, symSetWeylTwo, subSetFlat);
    WeylFlat(j,5) = temp;
end
clearvars j temp myNum myVar myChar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(1) w(x,u) = lamba0 +lambda1*u +K*u^2 -conj(lambda1)*u^3
%   + conj(lambda0)*u^4;
syms lambda0 lambda1 K
syms dlambda0_x dlambda0_y dlambda0_z
syms dlambda1_x dlambda1_y dlambda1_z
syms dK_x dK_y dK_z

dlambda0Vec = [dlambda0_x; dlambda0_y; dlambda0_z];
dlambda0_Mu = [mu1,mu2,mu3]*dlambda0Vec;
dlambda0_conjMu = [conj(mu1),conj(mu2),conj(mu3)]*dlambda0Vec;
dlambda0_vnormv = [v1normv,v2normv,v3normv]*dlambda0Vec;

dlambda1Vec = [dlambda1_x; dlambda1_y; dlambda1_z];
dlambda1_Mu = [mu1,mu2,mu3]*dlambda1Vec;
dlambda1_conjMu = [conj(mu1),conj(mu2),conj(mu3)]*dlambda1Vec;
dlambda1_vnormv = [v1normv,v2normv,v3normv]*dlambda1Vec;

dKVec = [dK_x; dK_y; dK_z];
dK_Mu = [mu1,mu2,mu3]*dKVec;
dK_conjMu = [conj(mu1),conj(mu2),conj(mu3)]*dKVec;
dK_vnormv = [v1normv,v2normv,v3normv]*dKVec;

muVec = [mu1;mu2;mu3];
conjMuVec = [conj(mu1);conj(mu2);conj(mu3)];
vnormv = [v1normv; v2normv; v3normv];

% % lambda0
syms d2lambda0_xx d2lambda0_xy d2lambda0_xz
syms d2lambda0_yy d2lambda0_yz d2lambda0_zz
d2lambda0Mat = [d2lambda0_xx, d2lambda0_xy, d2lambda0_xz;
    d2lambda0_xy, d2lambda0_yy, d2lambda0_yz;
    d2lambda0_xz, d2lambda0_yz, d2lambda0_zz]; 

d2lambda0_MuMu = transpose(muVec)*d2lambda0Mat*muVec;
d2lambda0_MuconjMu = transpose(muVec)*d2lambda0Mat*conjMuVec;
d2lambda0_Muvnormv = transpose(muVec)*d2lambda0Mat*vnormv;
d2lambda0_conjMuMu = d2lambda0_MuconjMu;
d2lambda0_conjMuconjMu = transpose(conjMuVec)*d2lambda0Mat*conjMuVec;
d2lambda0_conjMuvnormv = transpose(conjMuVec)*d2lambda0Mat*vnormv;
d2lambda0_vnormvMu = d2lambda0_Muvnormv;
d2lambda0_vnormvconjMu = d2lambda0_conjMuvnormv;
d2lambda0_vnormvvnormv = transpose(vnormv)*d2lambda0Mat*vnormv;

d2lambda0_Muu = 2*conj(u)/(1+u*conj(u))*dlambda0_Mu...
    + 2*dlambda0_vnormv;
d2lambda0_vnormvu = -1/(1+u*conj(u))^2*dlambda0_conjMu;
d2lambda0_conjMuconju = 2*u/(1+u*conj(u))*dlambda0_conjMu...
    + 2*dlambda0_vnormv;
d2lambda0_vnormvconju = -1/(1+u*conj(u))^2*dlambda0_Mu;

% % lambda1
syms d2lambda1_xx d2lambda1_xy d2lambda1_xz
syms d2lambda1_yy d2lambda1_yz d2lambda1_zz
d2lambda1Mat = [d2lambda1_xx, d2lambda1_xy, d2lambda1_xz;
    d2lambda1_xy, d2lambda1_yy, d2lambda1_yz;
    d2lambda1_xz, d2lambda1_yz, d2lambda1_zz]; 

d2lambda1_MuMu = transpose(muVec)*d2lambda1Mat*muVec;
d2lambda1_MuconjMu = transpose(muVec)*d2lambda1Mat*conjMuVec;
d2lambda1_Muvnormv = transpose(muVec)*d2lambda1Mat*vnormv;
d2lambda1_conjMuMu = d2lambda1_MuconjMu;
d2lambda1_conjMuconjMu = transpose(conjMuVec)*d2lambda1Mat*conjMuVec;
d2lambda1_conjMuvnormv = transpose(conjMuVec)*d2lambda1Mat*vnormv;
d2lambda1_vnormvMu = d2lambda1_Muvnormv;
d2lambda1_vnormvconjMu = d2lambda1_conjMuvnormv;
d2lambda1_vnormvvnormv = transpose(vnormv)*d2lambda1Mat*vnormv;

d2lambda1_Muu = 2*conj(u)/(1+u*conj(u))*dlambda1_Mu...
    + 2*dlambda1_vnormv;
d2lambda1_vnormvu = -1/(1+u*conj(u))^2*dlambda1_conjMu;
d2lambda1_conjMuconju = 2*u/(1+u*conj(u))*dlambda1_conjMu...
    + 2*dlambda1_vnormv;
d2lambda1_vnormvconju = -1/(1+u*conj(u))^2*dlambda1_Mu;

% % K (K is real-valued.)
syms d2K_xx d2K_xy d2K_xz
syms d2K_yy d2K_yz d2K_zz
d2KMat = [d2K_xx, d2K_xy, d2K_xz;
    d2K_xy, d2K_yy, d2K_yz;
    d2K_xz, d2K_yz, d2K_zz]; 

d2K_MuMu = transpose(muVec)*d2KMat*muVec;
d2K_MuconjMu = transpose(muVec)*d2KMat*conjMuVec;
d2K_Muvnormv = transpose(muVec)*d2KMat*vnormv;
d2K_conjMuMu = d2K_MuconjMu;
d2K_conjMuconjMu = transpose(conjMuVec)*d2KMat*conjMuVec;
d2K_conjMuvnormv = transpose(conjMuVec)*d2KMat*vnormv;
d2K_vnormvMu = d2K_Muvnormv;
d2K_vnormvconjMu = d2K_conjMuvnormv;
d2K_vnormvvnormv = transpose(vnormv)*d2KMat*vnormv;

d2K_Muu = 2*conj(u)/(1+u*conj(u))*dK_Mu...
    + 2*dK_vnormv;
d2K_vnormvu = -1/(1+u*conj(u))^2*dK_conjMu;
d2K_conjMuconju = 2*u/(1+u*conj(u))*dK_conjMu...
    + 2*dK_vnormv;
d2K_vnormvconju = -1/(1+u*conj(u))^2*dK_Mu;

% derivativeDict in the format [Mu,conjMu,vnormv,u,conju];
CVarWeyl6 = [u, lambda0, lambda1, K, ...
    dlambda0_x, dlambda0_y, dlambda0_z, ...
    dlambda1_x, dlambda1_y, dlambda1_z, ...
    dK_x, dK_y, dK_z];

derivativeDict.u = [0; 0; 0; 1; 0];
derivativeDict.lambda0 = [dlambda0_Mu; dlambda0_conjMu; dlambda0_vnormv;
    0; 0];
derivativeDict.lambda1 = [dlambda1_Mu; dlambda1_conjMu; dlambda1_vnormv;
    0; 0];
derivativeDict.K = [dK_Mu; dK_conjMu; dK_vnormv; 0; 0];
% % lambda0
derivativeDict.dlambda0_x = [d2lambda0Mat(1,:)*muVec;
    d2lambda0Mat(1,:)*conjMuVec;
    d2lambda0Mat(1,:)*vnormv; 0; 0];
derivativeDict.dlambda0_y = [d2lambda0Mat(2,:)*muVec;
    d2lambda0Mat(2,:)*conjMuVec;
    d2lambda0Mat(2,:)*vnormv; 0; 0];
derivativeDict.dlambda0_z = [d2lambda0Mat(3,:)*muVec;
    d2lambda0Mat(3,:)*conjMuVec;
    d2lambda0Mat(3,:)*vnormv; 0; 0];
% % lambda1
derivativeDict.dlambda1_x = [d2lambda1Mat(1,:)*muVec;
    d2lambda1Mat(1,:)*conjMuVec;
    d2lambda1Mat(1,:)*vnormv; 0; 0];
derivativeDict.dlambda1_y = [d2lambda1Mat(2,:)*muVec;
    d2lambda1Mat(2,:)*conjMuVec;
    d2lambda1Mat(2,:)*vnormv; 0; 0];
derivativeDict.dlambda1_z = [d2lambda1Mat(3,:)*muVec;
    d2lambda1Mat(3,:)*conjMuVec;
    d2lambda1Mat(3,:)*vnormv; 0; 0];
% % K
derivativeDict.dK_x = [d2KMat(1,:)*muVec;
    d2KMat(1,:)*conjMuVec;
    d2KMat(1,:)*vnormv; 0; 0];
derivativeDict.dK_y = [d2KMat(2,:)*muVec;
    d2KMat(2,:)*conjMuVec;
    d2KMat(2,:)*vnormv; 0; 0];
derivativeDict.dK_z = [d2KMat(3,:)*muVec;
    d2KMat(3,:)*conjMuVec;
    d2KMat(3,:)*vnormv; 0; 0];

Kset = [K, dK_x, dK_y, dK_z, d2K_xx, d2K_xy, d2K_xz, d2K_yy,...
    d2K_yz, d2K_zz];
conjKset = [conj(K), conj(dK_x), conj(dK_y), conj(dK_z),...
    conj(d2K_xx), conj(d2K_xy), conj(d2K_xz), conj(d2K_yy),...
    conj(d2K_yz), conj(d2K_zz)];

wSub = lambda0 +lambda1*u+ K*u^2 -conj(lambda1)*u^3 +conj(lambda0)*u^4;
dwVec = df_main4_CR_funWv2(wSub ,CVarWeyl6, derivativeDict, gamma);
dw_MuSub = dwVec(1);
dw_conjMuSub = dwVec(2);
dw_vnormvSub = dwVec(3);
dw_uSub = dwVec(4);

variableSet1 = [w, dw_Mu, dw_conjMu, dw_vnormv, dw_u];
subSet1 = [wSub, dw_MuSub, dw_conjMuSub, dw_vnormvSub, dw_uSub];
% %
d2w_uVec = df_main4_CR_funWv2(dw_uSub ,CVarWeyl6, derivativeDict, gamma);
d2w_uMuSub = d2w_uVec(1);
d2w_uconjMuSub = d2w_uVec(2);
d2w_uvnormvSub = d2w_uVec(3);
d2w_uuSub = d2w_uVec(4);

d2w_MuVec = df_main4_CR_funWv2(dw_MuSub ,CVarWeyl6, derivativeDict, gamma);
d2w_MuMuSub = d2w_MuVec(1);
d2w_MuconjMuSub = d2w_MuVec(2);
d2w_MuvnormvSub = d2w_MuVec(3);
d2w_conjMuVec = df_main4_CR_funWv2(dw_conjMuSub ,CVarWeyl6, derivativeDict, gamma);
d2w_conjMuconjMuSub = d2w_conjMuVec(2);
d2w_conjMuvnormvSub = d2w_conjMuVec(3);
d2w_vnormvVec = df_main4_CR_funWv2(dw_vnormvSub ,CVarWeyl6, derivativeDict, gamma);
d2w_vnormvvnormvSub = d2w_vnormvVec(3);

variableSet2 = [d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_uu,...
    d2w_MuMu, d2w_MuconjMu, d2w_Muvnormv, d2w_conjMuconjMu,...
    d2w_conjMuvnormv, d2w_vnormvvnormv];
subSet2 = [d2w_uMuSub, d2w_uconjMuSub, d2w_uvnormvSub, d2w_uuSub,...
    d2w_MuMuSub, d2w_MuconjMuSub, d2w_MuvnormvSub, d2w_conjMuconjMuSub,...
    d2w_conjMuvnormvSub, d2w_vnormvvnormvSub];
% %
d3w_uuVec = df_main4_CR_funWv2(d2w_uuSub ,CVarWeyl6, derivativeDict, gamma);
d3w_uuMuSub = d3w_uuVec(1);
d3w_uuconjMuSub = d3w_uuVec(2);
d3w_uuvnormvSub = d3w_uuVec(3);
d3w_uuuSub = d3w_uuVec(4);

d3w_uMuVec = df_main4_CR_funWv2(d2w_uMuSub ,CVarWeyl6, derivativeDict, gamma);
d3w_uMuMuSub = d3w_uMuVec(1);
d3w_uMuconjMuSub = d3w_uMuVec(2);
d3w_uMuvnormvSub = d3w_uMuVec(3);

d3w_uconjMuVec = df_main4_CR_funWv2(d2w_uconjMuSub ,CVarWeyl6, derivativeDict, gamma);
d3w_uconjMuconjMuSub = d3w_uconjMuVec(2);
d3w_uconjMuvnormvSub = d3w_uconjMuVec(3);

variableSet3 = [d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv, d3w_uuu,...
    d3w_uMuMu, d3w_uMuconjMu, d3w_uMuvnormv, d3w_uconjMuconjMu,...
    d3w_uconjMuvnormv];

subSet3 = [d3w_uuMuSub, d3w_uuconjMuSub, d3w_uuvnormvSub, d3w_uuuSub,...
    d3w_uMuMuSub, d3w_uMuconjMuSub, d3w_uMuvnormvSub, d3w_uconjMuconjMuSub,...
    d3w_uconjMuvnormvSub];
% %
d4w_uuuVec = df_main4_CR_funWv2(d3w_uuuSub ,CVarWeyl6, derivativeDict, gamma);
d4w_uuuMuSub = d4w_uuuVec(1);
d4w_uuuconjMuSub = d4w_uuuVec(2);
d4w_uuuuSub = d4w_uuuVec(4);

d4w_uuMuVec = df_main4_CR_funWv2(d3w_uuMuSub ,CVarWeyl6, derivativeDict, gamma);
d4w_uuconjMuVec = df_main4_CR_funWv2(d3w_uuconjMuSub ,CVarWeyl6, derivativeDict, gamma);
d4w_uuMuMuSub = d4w_uuMuVec(1);
d4w_uuMuconjMuSub = d4w_uuMuVec(2);
d4w_uuconjMuconjMuSub = d4w_uuconjMuVec(2);

variableSet4 = [d4w_uuuMu, d4w_uuuconjMu, d4w_uuuu,...
    d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuconjMuconjMu];
subSet4 = [d4w_uuuMuSub, d4w_uuuconjMuSub, d4w_uuuuSub,...
    d4w_uuMuMuSub, d4w_uuMuconjMuSub, d4w_uuconjMuconjMuSub];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MVarW6 = [];
for j=1:120
    temp = WeylFlat(j,5);
    temp = subs(temp, variableSet4, subSet4);
    temp = subs(temp, variableSet3, subSet3);
    temp = subs(temp, variableSet2, subSet2);
    temp = subs(temp, variableSet1, subSet1);
    temp = subs(temp, conjKset, Kset);
    tempSet = symvar(temp);
    MVarW6 = union(MVarW6, tempSet);
    WeylFlat(j,5) = temp;
end

for j=1:120
    temp = WeylFlat(j,5);
    WeylFlat(j,5) = complex_simple3(temp, MVarW6);
end

clearvars j temp tempSet 
clearvars variableSet1 variableSet2 variableSet3 variableSet4
clearvars subSet1 subSet2 subSet3 subSet4
clearvars MVarWeyls1 MVarWeyls2 

clearvars d4w_uuconjMuVec d3aV_conjuconjuVec d2w_MuVec d2w_vnormvVec 
clearvars wSetTwo d4w_uuMuVec d3w_uconjMuVec d3w_uMuVec d2w_conjMuVec
save('DataWeyl6_NatTwoR_CR_rmW.mat');

%% Test on WeylFlat
W1212 = WeylFlat(1,5); %C1111, Weyl(u1,u2,u1,u2)
W1214 = WeylFlat(3,5); %C1112, Weyl(u1,u2,u1,u4)
W1223 = WeylFlat(6,5); %-C1121, -Weyl(u1,u2,u3,u2)
W1234 = WeylFlat(10,5); %C1122, Weyl(u1,u2,u3,u4)
W1414 = WeylFlat(30,5); %C1212, Weyl(u1,u4,u1,u4)
W1423 = WeylFlat(33,5); %-C1221, -Weyl(u1,u4,u3,u2)
W1434 = WeylFlat(37,5); %C1222, Weyl(u1,u4,u3,u4)
W2323 = WeylFlat(66,5); %C2121, Weyl(u3,u2,u3,u2)
W2334 = WeylFlat(70,5); %-C2122, -Weyl(u3,u2,u3,u4)
W3434 = WeylFlat(100,5); %C2222, Weyl(u3,u4,u3,u4)

% rho-drho terms
W1435 = WeylFlat(38,5);
W1445 = WeylFlat(40,5);
W1534 = WeylFlat(49,5);
W2335 = WeylFlat(71,5);
W2345 = WeylFlat(73,5);
W2534 = WeylFlat(87,5);
W3535 = WeylFlat(106,5);
W3545 = WeylFlat(108,5);
W4545 = WeylFlat(115,5);
