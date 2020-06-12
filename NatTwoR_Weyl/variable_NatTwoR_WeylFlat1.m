%w(x,u) = lamba0+lambda1*u+K*u^2-conj(lambda1)*u^3+conj(lambda0)*u^4;
syms lambda0 lambda1 K
syms dlambda0_x dlambda0_y dlambda0_z dlambda1_x dlambda1_y dlambda1_z
syms dK_x dK_y dK_z
muVec = [mu1;mu2;mu3];
conjMuVec = [conj(mu1);conj(mu2);conj(mu3)];
vnormvVec = [v1normv; v2normv; v3normv];
syms d2lambda0_xx d2lambda0_xy d2lambda0_xz
syms d2lambda0_yy d2lambda0_yz d2lambda0_zz
syms d2lambda1_xx d2lambda1_xy d2lambda1_xz
syms d2lambda1_yy d2lambda1_yz d2lambda1_zz
syms d2K_xx d2K_xy d2K_xz d2K_yy d2K_yz d2K_zz
d2lambda0Mat = [d2lambda0_xx, d2lambda0_xy, d2lambda0_xz;
    d2lambda0_xy, d2lambda0_yy, d2lambda0_yz;
    d2lambda0_xz, d2lambda0_yz, d2lambda0_zz]; 
d2lambda1Mat = [d2lambda1_xx, d2lambda1_xy, d2lambda1_xz;
    d2lambda1_xy, d2lambda1_yy, d2lambda1_yz;
    d2lambda1_xz, d2lambda1_yz, d2lambda1_zz]; 
d2KMat = [d2K_xx, d2K_xy, d2K_xz;
    d2K_xy, d2K_yy, d2K_yz;
    d2K_xz, d2K_yz, d2K_zz];
MVarFlat1 = [u, lambda0, lambda1, K,...
    dlambda0_x, dlambda0_y, dlambda0_z, dlambda1_x, dlambda1_y, dlambda1_z,...
    dK_x, dK_y, dK_z, d2lambda0_xx, d2lambda0_xy, d2lambda0_xz,...
    d2lambda0_yy, d2lambda0_yz, d2lambda0_zz,...
    d2lambda1_xx, d2lambda1_xy, d2lambda1_xz,...
    d2lambda1_yy, d2lambda1_yz, d2lambda1_zz,...
    d2K_xx, d2K_xy, d2K_xz, d2K_yy, d2K_yz, d2K_zz];
Kset = [K, dK_x, dK_y, dK_z, d2K_xx, d2K_xy, d2K_xz, d2K_yy, d2K_yz, d2K_zz];
conjKset = [conj(K), conj(dK_x), conj(dK_y), conj(dK_z),...
    conj(d2K_xx), conj(d2K_xy), conj(d2K_xz), conj(d2K_yy),...
    conj(d2K_yz), conj(d2K_zz)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlambda0Vec = [dlambda0_x; dlambda0_y; dlambda0_z];
dlambda0_Mu = transpose(muVec)*dlambda0Vec;
dlambda0_conjMu = transpose(conjMuVec)*dlambda0Vec;
dlambda0_vnormv = transpose(vnormvVec)*dlambda0Vec;
dlambda1Vec = [dlambda1_x; dlambda1_y; dlambda1_z];
dlambda1_Mu = transpose(muVec)*dlambda1Vec;
dlambda1_conjMu = transpose(conjMuVec)*dlambda1Vec;
dlambda1_vnormv = transpose(vnormvVec)*dlambda1Vec;
dKVec = [dK_x; dK_y; dK_z];
dK_Mu = transpose(muVec)*dKVec;
dK_conjMu = conj(dK_Mu);
dK_vnormv = transpose(vnormvVec)*dKVec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2lambda1_MuMu = transpose(muVec)*d2lambda1Mat*muVec;
d2lambda1_MuconjMu = transpose(muVec)*d2lambda1Mat*conjMuVec;
d2lambda1_Muvnormv = transpose(muVec)*d2lambda1Mat*vnormvVec;
d2lambda1_conjMuMu = d2lambda1_MuconjMu;
d2lambda1_conjMuconjMu = transpose(conjMuVec)*d2lambda1Mat*conjMuVec;
d2lambda1_conjMuvnormv = transpose(conjMuVec)*d2lambda1Mat*vnormvVec;
d2lambda1_vnormvMu = d2lambda1_Muvnormv;
d2lambda1_vnormvconjMu = d2lambda1_conjMuvnormv;
d2lambda1_vnormvvnormv = transpose(vnormvVec)*d2lambda1Mat*vnormvVec;
d2lambda1_Muu = 2*conj(u)/(1+u*conj(u))*dlambda1_Mu...
    + 2*dlambda1_vnormv;
d2lambda1_vnormvu = -1/(1+u*conj(u))^2*dlambda1_conjMu;
d2lambda1_conjMuconju = 2*u/(1+u*conj(u))*dlambda1_conjMu...
    + 2*dlambda1_vnormv;
d2lambda1_vnormvconju = -1/(1+u*conj(u))^2*dlambda1_Mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2lambda0_MuMu = transpose(muVec)*d2lambda0Mat*muVec;
d2lambda0_MuconjMu = transpose(muVec)*d2lambda0Mat*conjMuVec;
d2lambda0_Muvnormv = transpose(muVec)*d2lambda0Mat*vnormvVec;
d2lambda0_conjMuMu = d2lambda0_MuconjMu;
d2lambda0_conjMuconjMu = transpose(conjMuVec)*d2lambda0Mat*conjMuVec;
d2lambda0_conjMuvnormv = transpose(conjMuVec)*d2lambda0Mat*vnormvVec;
d2lambda0_vnormvMu = d2lambda0_Muvnormv;
d2lambda0_vnormvconjMu = d2lambda0_conjMuvnormv;
d2lambda0_vnormvvnormv = transpose(vnormvVec)*d2lambda0Mat*vnormvVec;
d2lambda0_Muu = 2*conj(u)/(1+u*conj(u))*dlambda0_Mu...
    + 2*dlambda0_vnormv;
d2lambda0_vnormvu = -1/(1+u*conj(u))^2*dlambda0_conjMu;
d2lambda0_conjMuconju = 2*u/(1+u*conj(u))*dlambda0_conjMu...
    + 2*dlambda0_vnormv;
d2lambda0_vnormvconju = -1/(1+u*conj(u))^2*dlambda0_Mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2K_MuMu = transpose(muVec)*d2KMat*muVec;
d2K_MuconjMu = transpose(muVec)*d2KMat*conjMuVec;
d2K_Muvnormv = transpose(muVec)*d2KMat*vnormvVec;
d2K_conjMuMu = d2K_MuconjMu;
d2K_conjMuconjMu = transpose(conjMuVec)*d2KMat*conjMuVec;
d2K_conjMuvnormv = transpose(conjMuVec)*d2KMat*vnormvVec;
d2K_vnormvMu = d2K_Muvnormv;
d2K_vnormvconjMu = d2K_conjMuvnormv;
d2K_vnormvvnormv = transpose(vnormvVec)*d2KMat*vnormvVec;
d2K_Muu = 2*conj(u)/(1+u*conj(u))*dK_Mu +2*dK_vnormv;
d2K_vnormvu = -1/(1+u*conj(u))^2*dK_conjMu;
d2K_conjMuconju = conj(d2K_Muu);
d2K_vnormvconju = conj(d2K_vnormvu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVarFlat1 = [u, lambda0, lambda1, K, dlambda0_x, dlambda0_y, dlambda0_z, ...
    dlambda1_x, dlambda1_y, dlambda1_z, dK_x, dK_y, dK_z];
derivativeDict.u = [0; 0; 0; 1; 0];
derivativeDict.lambda0 = [dlambda0_Mu; dlambda0_conjMu; dlambda0_vnormv;
    0; 0];
derivativeDict.lambda1 = [dlambda1_Mu; dlambda1_conjMu; dlambda1_vnormv;
    0; 0];
derivativeDict.K = [dK_Mu; dK_conjMu; dK_vnormv; 0; 0];
derivativeDict.dlambda0_x = [d2lambda0Mat(1,:)*muVec;
    d2lambda0Mat(1,:)*conjMuVec;
    d2lambda0Mat(1,:)*vnormvVec; 0; 0];
derivativeDict.dlambda0_y = [d2lambda0Mat(2,:)*muVec;
    d2lambda0Mat(2,:)*conjMuVec;
    d2lambda0Mat(2,:)*vnormvVec; 0; 0];
derivativeDict.dlambda0_z = [d2lambda0Mat(3,:)*muVec;
    d2lambda0Mat(3,:)*conjMuVec;
    d2lambda0Mat(3,:)*vnormvVec; 0; 0];
derivativeDict.dlambda1_x = [d2lambda1Mat(1,:)*muVec;
    d2lambda1Mat(1,:)*conjMuVec;
    d2lambda1Mat(1,:)*vnormvVec; 0; 0];
derivativeDict.dlambda1_y = [d2lambda1Mat(2,:)*muVec;
    d2lambda1Mat(2,:)*conjMuVec;
    d2lambda1Mat(2,:)*vnormvVec; 0; 0];
derivativeDict.dlambda1_z = [d2lambda1Mat(3,:)*muVec;
    d2lambda1Mat(3,:)*conjMuVec;
    d2lambda1Mat(3,:)*vnormvVec; 0; 0];
derivativeDict.dK_x = [d2KMat(1,:)*muVec;
    d2KMat(1,:)*conjMuVec;
    d2KMat(1,:)*vnormvVec; 0; 0];
derivativeDict.dK_y = [d2KMat(2,:)*muVec;
    d2KMat(2,:)*conjMuVec;
    d2KMat(2,:)*vnormvVec; 0; 0];
derivativeDict.dK_z = [d2KMat(3,:)*muVec;
    d2KMat(3,:)*conjMuVec;
    d2KMat(3,:)*vnormvVec; 0; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
