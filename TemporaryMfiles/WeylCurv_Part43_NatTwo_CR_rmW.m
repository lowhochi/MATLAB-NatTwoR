% WeylCurv_Part43_NatTwo_CR_rmW.m
load('Data_WeylChern_Part4_May30.mat');
assumeAlso([x y z u1 u2 gamma K], 'real');
countLinear = [1, 2, 1, 4;
    1, 2, 2, 3;
    1, 2, 3, 5;
    1, 2, 4, 5;
    1, 4, 1, 5;
    1, 4, 2, 5;
    1, 5, 2, 3;
    1, 5, 3, 5;
    1, 5, 4, 5;
    2, 3, 2, 5;
    2, 5, 3, 5;
    2, 5, 4, 5];
countNonLinear = [1, 2, 1, 2;
    1, 2, 1, 5;
    1, 2, 2, 5;
    1, 5, 1, 5;
    1, 5, 2, 5;
    2, 5, 2, 5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms dK_y d2K_yy d3K_yyy real
syms Psi dPsi_p d2Psi_pp d3Psi_ppp
syms Zeta dZeta_p d2Zeta_pp d3Zeta_ppp
syms conju
MVarSPsi = [u, Psi, Zeta, Psi, dPsi_p, d2Psi_pp, d3Psi_ppp, ...
    dZeta_p, d2Zeta_pp, d3Zeta_ppp];
subSetPsi = sym('subSetPsi', [33,1]); % length(symSetWeylS)=33
muVector = [mu1, mu2, mu3];
conjMuVector = [conj(mu1), conj(mu2), conj(mu3)];
vnormv = [(u+conj(u))/(1+u*conj(u)), (1-u*conj(u))/(1+u*conj(u)), ...
    i*(u-conj(u))/(1+u*conj(u))];

dlambda0_x = -y/2*d2Psi_pp + dZeta_p;
dlambda0_y = -1/2*dPsi_p;
dlambda0_z = -i/2*y*d2Psi_pp + i*dZeta_p;
dlambda1_x = dPsi_p;
dlambda1_y = 0;
dlambda1_z = i*dPsi_p;
dK_x = 0; dK_z = 0; % dK_y

d2lambda1_xx = d2Psi_pp;
d2lambda1_xy = 0;
d2lambda1_xz = i*d2Psi_pp;
d2lambda1_yy = 0;
d2lambda1_yz = 0; 
d2lambda1_zz = -d2Psi_pp;
d2lambda1Temp = [d2lambda1_xx, d2lambda1_xy, d2lambda1_xz;
    d2lambda1_xy, d2lambda1_yy, d2lambda1_yz;
    d2lambda1_xz, d2lambda1_yz, d2lambda1_zz];
d2lambda0_xx = -y/2*d3Psi_ppp + d2Zeta_pp;
d2lambda0_xy = -1/2*d2Psi_pp; 
d2lambda0_xz = -i/2*y*d3Psi_ppp + i*d2Zeta_pp;
d2lambda0_yy = 0;
d2lambda0_yz = -i/2*d2Psi_pp;
d2lambda0_zz = 1/2*y*d3Psi_ppp - d2Zeta_pp;
d2lambda0Temp = [d2lambda0_xx, d2lambda0_xy, d2lambda0_xz;
    d2lambda0_xy, d2lambda0_yy, d2lambda0_yz;
    d2lambda0_xz, d2lambda0_yz, d2lambda0_zz];
d2K_xx = 0; d2K_xy = 0; d2K_xz = 0; % d2K_yy
d2K_yz = 0; d2K_zz = 0;
d2KTemp = [0, 0, 0;
    0, d2K_yy, 0; 
    0, 0, 0];

for j=1:33
    variable = string(symSetWeylS(j));
    switch variable
        case string(u)
            subSetPsi(j) = u;
        case string(lambda0)
            subSetPsi(j) = -y/2*dPsi_p;
        case string(lambda1)
            subSetPsi(j) = Psi;
        case string(K)
            subSetPsi(j) = K;
        case string(dlambda0_Mu)
            subSetPsi(j) = muVector*[dlambda0_x;dlambda0_y;dlambda0_z];
        case string(dlambda0_conjMu)
            subSetPsi(j) = conjMuVector*[dlambda0_x;dlambda0_y;dlambda0_z];
        case string(dlambda0_vnormv)
            subSetPsi(j) = vnormv*[dlambda0_x;dlambda0_y;dlambda0_z];
        case string(dlambda1_Mu)
            subSetPsi(j) = muVector*[dlambda1_x;dlambda1_y;dlambda1_z];
        case string(dlambda1_conjMu)
            subSetPsi(j) = conjMuVector*[dlambda1_x;dlambda1_y;dlambda1_z];
        case string(dlambda1_vnormv)
            subSetPsi(j) = vnormv*[dlambda1_x;dlambda1_y;dlambda1_z];
        case string(dK_Mu)
            subSetPsi(j) = muVector*[dK_x;dK_y;dK_z];
        case string(dK_conjMu)
            subSetPsi(j) = conjMuVector*[dK_x;dK_y;dK_z];
        case string(dK_vnormv)
            subSetPsi(j) = vnormv*[dK_x;dK_y;dK_z];  
        case string(d2lambda0_MuMu)
            subSetPsi(j) = muVector*d2lambda0Temp*transpose(muVector);
        case string(d2lambda0_MuconjMu)
            subSetPsi(j) = muVector*d2lambda0Temp*transpose(conjMuVector);
        case string(d2lambda0_Muvnormv)
            subSetPsi(j) = muVector*d2lambda0Temp*transpose(vnormv);
        case string(d2lambda0_conjMuconjMu)
            subSetPsi(j) = conjMuVector*d2lambda0Temp*transpose(conjMuVector);
        case string(d2lambda0_conjMuvnormv)
            subSetPsi(j) = conjMuVector*d2lambda0Temp*transpose(vnormv);
        case string(d2lambda0_vnormvvnormv)
            subSetPsi(j) = vnormv*d2lambda0Temp*transpose(vnormv);
        case string(d2lambda1_MuMu)
             subSetPsi(j) = muVector*d2lambda1Temp*transpose(muVector);
        case string(d2lambda1_MuconjMu)
            subSetPsi(j) = muVector*d2lambda1Temp*transpose(conjMuVector);
        case string(d2lambda1_Muvnormv)
            subSetPsi(j) = muVector*d2lambda1Temp*transpose(vnormv);
        case string(d2lambda1_conjMuconjMu)
            subSetPsi(j) = conjMuVector*d2lambda1Temp*transpose(conjMuVector);
        case string(d2lambda1_conjMuvnormv)
            subSetPsi(j) = conjMuVector*d2lambda1Temp*transpose(vnormv);
        case string(d2lambda1_vnormvvnormv)    
            subSetPsi(j) = vnormv*d2lambda1Temp*transpose(vnormv);
        case string(d2K_MuMu)
            subSetPsi(j) = muVector*d2KTemp*transpose(muVector);
        case string(d2K_MuconjMu)
            subSetPsi(j) = muVector*d2KTemp*transpose(conjMuVector);
        case string(d2K_Muvnormv)
            subSetPsi(j) = muVector*d2KTemp*transpose(vnormv);
        case string(d2K_conjMuconjMu)
            subSetPsi(j) = conjMuVector*d2KTemp*transpose(conjMuVector);
        case string(d2K_conjMuvnormv)
            subSetPsi(j) = conjMuVector*d2KTemp*transpose(vnormv);
        case string(d2K_vnormvvnormv)        
            subSetPsi(j) = vnormv*d2KTemp*transpose(vnormv);
        otherwise
            subSetPsi(j) = 0;
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WeylPsiL = sym('WeylPsiL',[12, 5]);
WeylPsiN = sym('WeylPsiN',[6, 5]);

for j=1:12 
    m = countLinear(j,1);
    n = countLinear(j,2);
    k = countLinear(j,3);
    ll = countLinear(j,4);
    temp = WeylS(m,n,k,ll);
    temp = subs(temp, symSetWeylS, subSetPsi);
    temp = complex_simple3(temp, MVarSPsi);
    WeylPsiL(j,:) = [m,n,k,ll, temp];
end

for j=1:6 
    m = countNonLinear(j,1);
    n = countNonLinear(j,2);
    k = countNonLinear(j,3);
    ll = countNonLinear(j,4);
    temp = WeylS(m,n,k,ll);
    temp = subs(temp, symSetWeylS, subSetPsi);
    temp = complex_simple3(temp, MVarSPsi);
    WeylPsiN(j,:) = [m,n,k,ll, temp];
end
clearvars m n k ll temp j variable
clearvars d2lambda0Temp d2lambda1Temp d2KTemp
save('Data_WeylChern_Part43_May31.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% section 2// Print WeylPsiN and WeylPsiL
load('Data_WeylChern_Part43_May31.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
assumeAlso([K, dK_y, d2K_yy, d3K_yyy], 'real');
latex_WeylCurv43 = fopen('latex_WeylCurv43.txt','w');
WP1212 = WeylPsiN(1,5)/(1+u*conj(u))^2;
WP1215 = WeylPsiN(2,5)/(1+u*conj(u));
WP1225 = WeylPsiN(3,5)/(1+u*conj(u));
WP1515 = WeylPsiN(4,5);
WP1525 = WeylPsiN(5,5);
WP2525 = WeylPsiN(6,5);

WP1212 = subs(WP1212, conj(u), conju);
WP1215 = subs(WP1215, conj(u), conju);
WP1225 = subs(WP1225, conj(u), conju);
WP1515 = subs(WP1515, conj(u), conju);
WP1525 = subs(WP1525, conj(u), conju);
WP2525 = subs(WP2525, conj(u), conju);
[c1212, u1212] = coeffs(WP1212, [u,conju,y]);
[c1215, u1215] = coeffs(WP1215, [u,conju,y]);
[c1225, u1225] = coeffs(WP1225, [u,conju,y]);
[c1515, u1515] = coeffs(WP1515, [u,conju,y]);
[c1525, u1525] = coeffs(WP1525, [u,conju,y]);
[c2525, u2525] = coeffs(WP2525, [u,conju,y]);

fprintf(latex_WeylCurv43,'WP1212 = $(1+u*conj(u))^2$*Weyl(1,2,1,2)\n');
for j=1:length(c1212)
  fprintf(latex_WeylCurv43,'$%s$ : $%s $\n', latex(u1212(j)), latex(c1212(j)));  
  fprintf(latex_WeylCurv43,'%s\n', ' ');
end

fprintf(latex_WeylCurv43,'WP1215 = $(1+u*conj(u))$*Weyl(1,2,1,5)\n');
for j=1:length(c1215)
  fprintf(latex_WeylCurv43,'$%s$ : $%s $ \n', latex(u1215(j)), latex(c1215(j)));  
  fprintf(latex_WeylCurv43,'%s\n', ' ');
end

fprintf(latex_WeylCurv43,'WP1225 = $(1+u*conj(u))$*Weyl(1,2,2,5)\n');
for j=1:length(c1225)
  fprintf(latex_WeylCurv43,'$%s$ : $%s $ \n', latex(u1225(j)), latex(c1225(j)));  
  fprintf(latex_WeylCurv43,'%s\n', ' ');
end

fprintf(latex_WeylCurv43,'WP1515 = Weyl(1,5,1,5)\n');
for j=1:length(c1515)
  fprintf(latex_WeylCurv43,'$%s$ : $%s $ \n', latex(u1515(j)), latex(c1515(j)));  
  fprintf(latex_WeylCurv43,'%s\n', ' ');
end

fprintf(latex_WeylCurv43,'WP1525 = Weyl(1,5,2,5)\n');
for j=1:length(c1525)
  fprintf(latex_WeylCurv43,'$%s$ : $%s $ \n', latex(u1525(j)), latex(c1525(j)));  
  fprintf(latex_WeylCurv43,'%s\n', ' ');
end

fprintf(latex_WeylCurv43,'WP2525 = Weyl(2,5,2,5)\n');
for j=1:length(c2525)
  fprintf(latex_WeylCurv43,'$%s$ : $%s $ \n', latex(u2525(j)), latex(c2525(j)));  
  fprintf(latex_WeylCurv43,'%s\n', ' ');
end

fclose(latex_WeylCurv43);