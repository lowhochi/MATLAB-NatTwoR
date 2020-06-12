load('DataWeyl6_NatTwoR_CR_rmW.mat');
for j=1:120
    temp = WeylFlat(j,5) ;
    temp = subs(temp,[d2theta_Muvnormv, d2theta_vnormvvnormv],[0,0]);
    WeylFlat(j,5) = temp;
end
W1212 = WeylFlat(1,5);
W1215 = WeylFlat(4,5);
W1515 = WeylFlat(43,5);
W1525 = WeylFlat(47,5);
nonLinearVec =[W1212, W1215, W1515, W1525];
nonLinearIndex = [1, 4, 43, 47];

variableSet = [lambda1, lambda0, K,...
    dlambda1_x, dlambda1_y, dlambda1_z, dlambda0_x, dlambda0_y, dlambda0_z,...
    dK_x, dK_y, dK_z, d2lambda1_xx, d2lambda1_xy, d2lambda1_xz, d2lambda1_yy,...
    d2lambda1_yz, d2lambda1_zz,  d2lambda0_xx, d2lambda0_xy, d2lambda0_xz,... 
    d2lambda0_yy, d2lambda0_yz, d2lambda0_zz,...
    d2K_xx, d2K_xy, d2K_xz, d2K_yy, d2K_yz, d2K_zz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms p y % this t means real variable x2 or y
syms alpha zeta1 zeta2 R % alpha and R are real-valued.

% derivatives of alpha
syms dalpha_p dalpha_y
dalpha_conjp = conj(dalpha_p);
syms d2alpha_pp d2alpha_pconjp d2alpha_py d2alpha_yy
syms d3alpha_ppp d3alpha_ppconjp d3alpha_yyy 
syms d3alpha_ppy d3alpha_pconjpy d3alpha_pyy
syms d4alpha_pppy d4alpha_ppconjpy d4alpha_ppyy d4alpha_pconjpyy d4alpha_pyyy 
syms d4alpha_pppp d4alpha_pppconjp d4alpha_ppconjpconjp
syms d4alpha_yyyy
% derivatives of zeta1
syms dzeta1_p d2zeta1_pp d3zeta1_ppp
% derivatives of zeta2
syms dzeta2_p d2zeta2_pp
% derivatives of R
syms dR_y d2R_yy

realVariable = [alpha, R, dalpha_y, d2alpha_pconjp, d2alpha_yy,...
    d3alpha_pconjpy, d3alpha_yyy, d4alpha_ppconjpconjp,...
    d4alpha_pconjpyy, d4alpha_yyyy, dR_y, d2R_yy];

MVarW63 = [u, p, alpha, zeta1, zeta2, R,...
    dalpha_p, dalpha_y, d2alpha_pp, d2alpha_pconjp, d2alpha_py, d2alpha_yy,...
    d3alpha_ppp, d3alpha_ppconjp, d3alpha_yyy,...
    d3alpha_ppy, d3alpha_pconjpy, d3alpha_pyy,...
    d4alpha_pppy, d4alpha_ppconjpy, d4alpha_ppyy, d4alpha_pconjpyy,...
    d4alpha_pyyy, d4alpha_pppp, d4alpha_pppconjp, d4alpha_ppconjpconjp,...
    d4alpha_yyyy, dzeta1_p, d2zeta1_pp, d3zeta1_ppp,...
    dzeta2_p, d2zeta2_pp, dR_y d2R_yy] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda1Sub = d2alpha_py;
lambda0Sub = -1/2*d2alpha_pp -i/2*(dzeta1_p*conj(p)+zeta2);
KSub = d2alpha_pconjp -1/2*d2alpha_yy - R +i*(zeta1-conj(zeta1));

dlambda1_p = d3alpha_ppy;
dlambda1_conjp = d3alpha_pconjpy;
dlambda1_xSub = dlambda1_p + dlambda1_conjp;
dlambda1_zSub = i*(dlambda1_p - dlambda1_conjp);
dlambda1_ySub = d3alpha_pyy;

dlambda0_p = -1/2*d3alpha_ppp -i/2*(d2zeta1_pp*conj(p)+dzeta2_p);
dlambda0_conjp = -1/2*d3alpha_ppconjp -i/2*dzeta1_p;
dlambda0_xSub = dlambda0_p + dlambda0_conjp;
dlambda0_zSub = i*(dlambda0_p-dlambda0_conjp);
dlambda0_ySub = -1/2*d3alpha_ppy;

dK_p = d3alpha_ppconjp -1/2*d3alpha_pyy + i*dzeta1_p;
dK_conjp = conj(dK_p);
dK_xSub = dK_p + dK_conjp;
dK_zSub = i*(dK_p - dK_conjp);
dK_ySub = d3alpha_pconjpy -1/2*d3alpha_yyy -dR_y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2lambda1_pp = d4alpha_pppy; 
d2lambda1_pconjp = d4alpha_ppconjpy;
d2lambda1_conjpconjp = conj(d4alpha_ppconjpy);
d2lambda1_py = d4alpha_ppyy;
d2lambda1_conjpy = d4alpha_pconjpyy;
d2lambda1_yySub = d4alpha_pyyy;

d2lambda1_xxSub = d2lambda1_pp +2*d2lambda1_pconjp +d2lambda1_conjpconjp; 
d2lambda1_xzSub = i*(d2lambda1_pp -d2lambda1_conjpconjp);
d2lambda1_zzSub = 2*d2lambda1_pconjp -d2lambda1_pp -d2lambda1_conjpconjp;
d2lambda1_xySub = d2lambda1_py + d2lambda1_conjpy;
d2lambda1_yzSub = i*(d2lambda1_py -d2lambda1_conjpy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2lambda0_pp = -1/2*d4alpha_pppp -i/2*(d3zeta1_ppp*conj(p) +d2zeta2_pp);
d2lambda0_pconjp = -1/2*d4alpha_pppconjp -i/2*d2zeta1_pp;
d2lambda0_conjpconjp = -1/2*d4alpha_ppconjpconjp;
d2lambda0_py = -1/2*d4alpha_pppy;
d2lambda0_conjpy = -1/2*d4alpha_ppconjpy;
d2lambda0_yySub = -1/2*d4alpha_ppyy;

d2lambda0_xxSub = d2lambda0_pp +2*d2lambda0_pconjp +d2lambda0_conjpconjp; 
d2lambda0_xzSub = i*(d2lambda0_pp -d2lambda0_conjpconjp);
d2lambda0_zzSub = 2*d2lambda0_pconjp -d2lambda0_pp -d2lambda0_conjpconjp;
d2lambda0_xySub = d2lambda0_py + d2lambda0_conjpy;
d2lambda0_yzSub = i*(d2lambda0_py -d2lambda0_conjpy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2K_pp = d4alpha_pppconjp -1/2*d4alpha_ppyy +i*d2zeta1_pp;
d2K_pconjp = d4alpha_ppconjpconjp -1/2*d4alpha_pconjpyy;
d2K_conjpconjp = conj(d2K_pp);
d2K_py = d4alpha_ppconjpy -1/2*d4alpha_pyyy;
d2K_conjpy = conj(d2K_py);
d2K_yySub = d4alpha_pconjpyy -1/2*d4alpha_yyyy -d2R_yy;

d2K_xxSub = d2K_pp +2*d2K_pconjp +d2K_conjpconjp; 
d2K_xzSub = i*(d2K_pp -d2K_conjpconjp);
d2K_zzSub = 2*d2K_pconjp -d2K_pp -d2K_conjpconjp;
d2K_xySub = d2K_py + d2K_conjpy;
d2K_yzSub = i*(d2K_py -d2K_conjpy);

subSet = [lambda1Sub, lambda0Sub, KSub,...
    dlambda1_xSub, dlambda1_ySub, dlambda1_zSub, dlambda0_xSub,...
    dlambda0_ySub, dlambda0_zSub, dK_xSub, dK_ySub, dK_zSub,...
    d2lambda1_xxSub, d2lambda1_xySub, d2lambda1_xzSub, d2lambda1_yySub,...
    d2lambda1_yzSub, d2lambda1_zzSub, d2lambda0_xxSub, d2lambda0_xySub,...
    d2lambda0_xzSub, d2lambda0_yySub, d2lambda0_yzSub, d2lambda0_zzSub,...
    d2K_xxSub, d2K_xySub, d2K_xzSub, d2K_yySub, d2K_yzSub, d2K_zzSub];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W1212 = WeylFlat(1,5);
% W1215 = WeylFlat(4,5);
% W1515 = WeylFlat(43,5);
% W1525 = WeylFlat(47,5);
% nonLinearVec =[W1212, W1215, W1515, W1525];
% nonLinearIndex = [1, 4, 43, 47];
for j=1:4
    temp = nonLinearVec(j);
    temp = subs(temp, variableSet, subSet);
    for q=1:12
        variable = realVariable(q);
        temp = subs(temp, conj(variable), variable);
    end
    clearvars q variable
    temp = complex_simple3(temp, MVarW63);
    nonLinearVec(j) = temp;
    clear temp
end
clear j m n k ll
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('DataWeyl6_Part3_NatTwoR_CR_rmW.mat');

%%
load('DataWeyl6_Part3_NatTwoR_CR_rmW.mat');
latexFile = fopen('latex_Nov15.txt','w');

for j=1:4
    temp = nonLinearVec(j);
    temp = subs(temp, u*conj(u), Y-1);
    temp = complex_simple3(temp, [MVarW63, Y]);
    nonLinearVec(j) = temp;
    clear temp
end

Wy1212 = nonLinearVec(1)/Y^2;
Wy1215 = nonLinearVec(2)/Y;
Wy1515 = nonLinearVec(3);
Wy1525 = nonLinearVec(4);

[term1212, u1212] = coeffs(Wy1212, [u,conj(u),Y]);
[term1215, u1215] = coeffs(Wy1215, [u,conj(u),Y]);
[term1515, u1515] = coeffs(Wy1515, [u,conj(u),Y]);
[term1525, u1525] = coeffs(Wy1525, [u,conj(u),Y]);

fprintf(latexFile, '$1/Y^2 \\cdot W_{1212}$ is \n');
fprintf(latexFile, '%s\n', ' '); 
for j=1:length(u1212)
    uLatex = latex(u1212(j));
    termLatex = latex(term1212(j));
    fprintf(latexFile, 'term $ %s $: \\\\ \n', uLatex);
    fprintf(latexFile, '%s\n', ' '); 
    fprintf(latexFile, '$ %s $ \\\\[0.1in] \n', termLatex);
    fprintf(latexFile, '%s\n', ' ');
end

fprintf(latexFile, '\\newpage \n');

fprintf(latexFile, '$1/Y \\cdot W_{1215}$ is \n');
fprintf(latexFile, '%s\n', ' '); 
for j=1:length(u1215)
    uLatex = latex(u1215(j));
    termLatex = latex(term1215(j));
    fprintf(latexFile, 'term $ %s $: \\\\ \n', uLatex);
    fprintf(latexFile, '%s\n', ' '); 
    fprintf(latexFile, '$ %s $ \\\\[0.1in] \n', termLatex);
    fprintf(latexFile, '%s\n', ' ');
end

fprintf(latexFile, '\\newpage \n');

fprintf(latexFile, '$W_{1515}$ is \n');
fprintf(latexFile, '%s\n', ' '); 
for j=1:length(u1515)
    uLatex = latex(u1515(j));
    termLatex = latex(term1515(j));
    fprintf(latexFile, 'term $ %s $: \\\\ \n', uLatex);
    fprintf(latexFile, '%s\n', ' '); 
    fprintf(latexFile, '$ %s $ \\\\[0.1in] \n', termLatex);
    fprintf(latexFile, '%s\n', ' ');
end

fprintf(latexFile, '\\newpage \n');

fprintf(latexFile, '$W_{1525}$ is \n');
fprintf(latexFile, '%s\n', ' '); 
for j=1:length(u1525)
    uLatex = latex(u1525(j));
    termLatex = latex(term1525(j));
    fprintf(latexFile, 'term $ %s $: \\\\ \n', uLatex);
    fprintf(latexFile, '%s\n', ' '); 
    fprintf(latexFile, '$ %s $ \\\\[0.1in] \n', termLatex);
    fprintf(latexFile, '%s\n', ' ');
end

fclose(latexFile);













