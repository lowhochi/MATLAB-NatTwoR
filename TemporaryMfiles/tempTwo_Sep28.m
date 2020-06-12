load('Data_ChernThree_Section2_Sep24_2019.mat')
% Check that Chern(1,1,1,1) is real
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H0= (2*theta + i*aV - i*conj(aV));
Y = (1+u*conj(u));

Rcheck1 = H0*Y^6*(daV_u*daV_conju);

Rcheck2 = Y^4*(-8*i*aV*theta^2 + 8*aV^2*theta + 2*i*aV^3 - 2*i*aV^2*conj(aV));

Rcheck3 = Y^4*H0*(-4*i*dtheta_vnormv)+Y^4*daV_u*(-2)*dtheta_Mu...
    +Y^4*H0*(-d2aV_conjuconjMu + 2*daV_vnormv);

Rcheck4 = Y^2*(2*d2theta_MuconjMu + i*d2aV_MuconjMu);

RcheckSum = (Rcheck1 + Rcheck2 + Rcheck3 + Rcheck4)/Y^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace variables in Gijk and derivatives.
variableRCheck = [aV, d2aV_MuconjMu, d2aV_conjuconjMu,...
    d2theta_MuconjMu, daV_conju, daV_u, daV_vnormv, dtheta_Mu,...
    dtheta_vnormv, theta]; % plus u

temp_Sep24_CR_rmW

thetaSub = G12_3 + G23_1 + G31_2;

dtheta_vnormvSub = dG12_3_vnormv + dG23_1_vnormv + dG31_2_vnormv;
dtheta_MuSub = dG12_3_Mu + dG23_1_Mu + dG31_2_Mu;
d2theta_MuconjMuSub = d2G12_3_MuconjMu + d2G23_1_MuconjMu + d2G31_2_MuconjMu;

aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);

daVVec = df_NatTwo_updated_CR_rmW(aVSub, CVarTest, dCVarTest);
daV_MuSub = daVVec(1);
daV_vnormvSub = daVVec(3);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);

d2aV_conjuVec = df_NatTwo_updated_CR_rmW(daV_conjuSub, CVarTest, dCVarTest);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);
d2aV_MuVec = df_NatTwo_updated_CR_rmW(daV_MuSub, CVarTest, dCVarTest);
d2aV_MuconjMuSub = d2aV_MuVec(2);

subSetRCheck = [aVSub, d2aV_MuconjMuSub, d2aV_conjuconjMuSub,...
    d2theta_MuconjMuSub, daV_conjuSub, daV_uSub, daV_vnormvSub, dtheta_MuSub,...
    dtheta_vnormvSub, thetaSub];

RcheckSum = subs(RcheckSum, variableRCheck, subSetRCheck);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variable2RCheck=[d2G11_2_MuconjMu, d2G11_3_MuconjMu, d2G12_3_MuconjMu,...
    d2G22_1_MuconjMu, d2G22_3_MuconjMu, d2G23_1_MuconjMu, d2G31_2_MuconjMu,...
    d2G33_1_MuconjMu, d2G33_2_MuconjMu, dG11_2_conjMu, dG11_2_vnormv,...
    dG11_3_conjMu, dG11_3_vnormv, dG12_3_Mu, dG12_3_conjMu, dG12_3_vnormv,...
    dG22_1_conjMu, dG22_1_vnormv, dG22_3_conjMu, dG22_3_vnormv, dG23_1_Mu,...
    dG23_1_conjMu, dG23_1_vnormv, dG31_2_Mu, dG31_2_conjMu, dG31_2_vnormv,...
    dG33_1_conjMu, dG33_1_vnormv, dG33_2_conjMu, dG33_2_vnormv];

subSet2RCheck = [d2G11_2_MuconjMuSub, d2G11_3_MuconjMuSub, d2G12_3_MuconjMuSub,...
    d2G22_1_MuconjMuSub, d2G22_3_MuconjMuSub, d2G23_1_MuconjMuSub, d2G31_2_MuconjMuSub,...
    d2G33_1_MuconjMuSub, d2G33_2_MuconjMuSub, dG11_2_conjMuSub, dG11_2_vnormvSub,...
    dG11_3_conjMuSub, dG11_3_vnormvSub, dG12_3_MuSub, dG12_3_conjMuSub, dG12_3_vnormvSub,...
    dG22_1_conjMuSub, dG22_1_vnormvSub, dG22_3_conjMuSub, dG22_3_vnormvSub, dG23_1_MuSub,...
    dG23_1_conjMuSub, dG23_1_vnormvSub, dG31_2_MuSub, dG31_2_conjMuSub, dG31_2_vnormvSub,...
    dG33_1_conjMuSub, dG33_1_vnormvSub, dG33_2_conjMuSub, dG33_2_vnormvSub];

Rcheck2Sum = subs(RcheckSum, variable2RCheck, subSet2RCheck);
Rcheck2Sum = complex_simple3(Rcheck2Sum,[u]);
save('DataTemp_Sep28.mat');
% latex_file = fopen('latex_Rcheck2_Sep29.txt','w');
% tempLatex = latex(Rcheck2Sum);
% fprintf(latex_file, '%s\n', tempLatex);
% fclose(latex_file);

%%
load('DataTemp_Sep28.mat');

dG11_2_by3Sub = dG23_1_by1 + dG11_3_by2 + dG31_2_by1...
    - G11_3*G33_2 + G11_2*G22_3 - (G22_1+G33_1)*(G31_2+G23_1);
% %
d2G11_2_by31Sub = d2G23_1_by11 + d2G11_3_by21 + d2G31_2_by11...
    - (dG11_3_by1*G33_2 + G11_3*dG33_2_by1)...
    + (dG11_2_by1*G22_3 + G11_2*dG22_3_by1)...
    - (dG22_1_by1 + dG33_1_by1)*(G31_2+G23_1)...
    - (G22_1+G33_1)*(dG31_2_by1+dG23_1_by1);
d2G11_2_by13Sub = d2G11_2_by31Sub + G11_3*dG11_2_by1...
    - G33_1*dG11_2_by3 + dG11_2_by2*(G12_3 + G31_2);

d2G11_2_by32Sub = d2G23_1_by12 + d2G11_3_by22 + d2G31_2_by12...
    - (dG11_3_by2*G33_2 + G11_3*dG33_2_by2)...
    + (dG11_2_by2*G22_3 + G11_2*dG22_3_by2)...
    - (dG22_1_by2 + dG33_1_by2)*(G31_2+G23_1)...
    - (G22_1+G33_1)*(dG31_2_by2+dG23_1_by2);
d2G11_2_by23Sub = d2G11_2_by32Sub + G22_3*dG11_2_by2 - G33_2*dG11_2_by3...
    - dG11_2_by1*(G23_1 + G31_2);

d2G11_2_by33Sub = d2G23_1_by13 + d2G11_3_by23 + d2G31_2_by13...
    - (dG11_3_by3*G33_2 + G11_3*dG33_2_by3)...
    + (dG11_2_by3*G22_3 + G11_2*dG22_3_by3)...
    - (dG22_1_by3 + dG33_1_by3)*(G31_2+G23_1)...
    - (G22_1+G33_1)*(dG31_2_by3+dG23_1_by3);

% %
dG22_3_by1Sub = dG12_3_by2 + dG31_2_by2 + dG22_1_by3...
    + G33_1*G22_3 - G22_1*G11_3 - (G11_2+G33_2)*(G31_2+G12_3);

d2G22_3_by11Sub = d2G12_3_by21 + d2G31_2_by21 + d2G22_1_by31...
    + dG33_1_by1*G22_3 + G33_1*dG22_3_by1...
    - dG22_1_by1*G11_3 - G22_1*dG11_3_by1...
    -(dG11_2_by1+dG33_2_by1)*(G31_2+G12_3)...
    -(G11_2+G33_2)*(dG31_2_by1 + dG12_3_by1);

d2G22_3_by12Sub = d2G12_3_by22 + d2G31_2_by22 + d2G22_1_by32...
    + dG33_1_by2*G22_3 + G33_1*dG22_3_by2...
    - dG22_1_by2*G11_3 - G22_1*dG11_3_by2...
    -(dG11_2_by2+dG33_2_by2)*(G31_2+G12_3)...
    -(G11_2+G33_2)*(dG31_2_by2 + dG12_3_by2);
 
d2G22_3_by13Sub = d2G12_3_by23 + d2G31_2_by23 + d2G22_1_by33...
    + dG33_1_by3*G22_3 + G33_1*dG22_3_by3...
    - dG22_1_by3*G11_3 - G22_1*dG11_3_by3...
    -(dG11_2_by3+dG33_2_by3)*(G31_2+G12_3)...
    -(G11_2+G33_2)*(dG31_2_by3 + dG12_3_by3);

% %
dG23_1_by3Sub = dG33_1_by2 - dG33_2_by1 - dG12_3_by3...
    - G11_2*G33_1 + G22_1*G33_2 + (G11_3+G22_3)*(G12_3+G23_1);

d2G23_1_by31Sub = d2G33_1_by21 - d2G33_2_by11 - d2G12_3_by31...
    - dG11_2_by1*G33_1 - G11_2*dG33_1_by1...
    + dG22_1_by1*G33_2 + G22_1*dG33_2_by1...
    + (dG11_3_by1 + dG22_3_by1)*(G12_3+G23_1)...
    + (G11_3 + G22_3)*(dG12_3_by1 + dG23_1_by1);
d2G23_1_by13Sub = d2G23_1_by31Sub + G11_3*dG23_1_by1 - G33_1*dG23_1_by3...
    + dG23_1_by2*(G12_3 + G31_2);


d2G23_1_by32Sub = d2G33_1_by22 - d2G33_2_by12 - d2G12_3_by32...
    - dG11_2_by2*G33_1 - G11_2*dG33_1_by2...
    + dG22_1_by2*G33_2 + G22_1*dG33_2_by2...
    + (dG11_3_by2 + dG22_3_by2)*(G12_3+G23_1)...
    + (G11_3 + G22_3)*(dG12_3_by2 + dG23_1_by2);

d2G23_1_by23Sub = d2G23_1_by32Sub + G22_3*dG23_1_by2 - G33_2*dG23_1_by3...
    - dG23_1_by1*(G23_1 + G31_2);

d2G23_1_by33Sub = d2G33_1_by23 - d2G33_2_by13 - d2G12_3_by33...
    - dG11_2_by3*G33_1 - G11_2*dG33_1_by3...
    + dG22_1_by3*G33_2 + G22_1*dG33_2_by3...
    + (dG11_3_by3 + dG22_3_by3)*(G12_3+G23_1)...
    + (G11_3 + G22_3)*(dG12_3_by3 + dG23_1_by3);

variableD2G = [d2G11_2_by13, d2G11_2_by23, d2G11_2_by33,...
    d2G22_3_by11, d2G22_3_by12, d2G22_3_by13,...
    d2G23_1_by13, d2G23_1_by23, d2G23_1_by33];

variableDG = [dG11_2_by3, dG22_3_by1, dG23_1_by3];

subSetD2G = [d2G11_2_by13Sub, d2G11_2_by23Sub, d2G11_2_by33Sub,...
    d2G22_3_by11Sub, d2G22_3_by12Sub, d2G22_3_by13Sub,...
    d2G23_1_by13Sub, d2G23_1_by23Sub, d2G23_1_by33Sub];

subSetDG = [dG11_2_by3Sub, dG22_3_by1Sub, dG23_1_by3Sub];

Rcheck2Sum = subs(Rcheck2Sum, variableD2G, subSetD2G);
Rcheck2Sum = subs(Rcheck2Sum, variableDG, subSetDG);

Rcheck2Sum = complex_simple3(Rcheck2Sum,[u]);

syms u1 u2 real
Rcheck2Sum = subs(Rcheck2Sum, u, u1+i*u2);
Rcheck2Sum = simplify(Rcheck2Sum);


d2GTempSet = [d2G11_2_by11, d2G11_2_by12, d2G11_2_by13,...
    d2G11_2_by22, d2G11_2_by23, d2G11_2_by33, d2G11_3_by11,...
    d2G11_3_by12, d2G11_3_by13, d2G11_3_by22, d2G11_3_by23,...
    d2G11_3_by33, d2G12_3_by11, d2G12_3_by12, d2G12_3_by13,...
    d2G12_3_by22, d2G12_3_by23, d2G12_3_by33, d2G22_1_by11,...
    d2G22_1_by12, d2G22_1_by13, d2G22_1_by22, d2G22_1_by23,...
    d2G22_1_by33, d2G22_3_by11, d2G22_3_by12, d2G22_3_by13,...
    d2G22_3_by22, d2G22_3_by23, d2G22_3_by33, d2G23_1_by11,...
    d2G23_1_by12, d2G23_1_by13, d2G23_1_by22, d2G23_1_by23,...
    d2G23_1_by33, d2G31_2_by11, d2G31_2_by12, d2G31_2_by13,...
    d2G31_2_by22, d2G31_2_by23, d2G31_2_by33, d2G33_1_by11,...
    d2G33_1_by12, d2G33_1_by13, d2G33_1_by22, d2G33_1_by23,...
    d2G33_1_by33, d2G33_2_by11, d2G33_2_by12, d2G33_2_by13,...
    d2G33_2_by22, d2G33_2_by23, d2G33_2_by33];

[termVec, gVec] = coeffs(Rcheck2Sum , d2GTempSet);

% latexDoc= fopen('latex_C1111real_Sep29_2019.txt','w');
% N = length(gVec);
% for j=1:N
%     temp0 = termVec(j);
%     temp0 = simplify(temp0);
%     temp0 = simplify(temp0);
%     fprintf(latexDoc,'The d2G-term is: $ %s $ \\\\ \n', latex(gVec(j)) );
%     fprintf(latexDoc,'%s \n', ' ');
%     fprintf(latexDoc, 'The term is: \\\\ \n');
%     fprintf(latexDoc,'$ %s $ \\\\[0.1in] \n', latex(temp0));
%     fprintf(latexDoc,'%s \n', ' ');
% end
% fclose(latexDoc);
N = length(gVec);
latexDocTwo= fopen('latex_C1111Two_Sep29_2019.txt','w');
imVec = sym('imVec',[1,N]);
for j=1:(N-1)
    imVec(j) = termVec(j) - conj(termVec(j));
    imVec(j) = expand(imVec(j));
    imVec(j) = simplify(imVec(j));
    fprintf(latexDocTwo,'The d2G-term is: $ %s $ \\\\ \n', latex(gVec(j)) );
    fprintf(latexDocTwo,'%s \n', ' ');
    fprintf(latexDocTwo, 'The term is: \\\\ \n');
    fprintf(latexDocTwo,'$ %s $ \\\\[0.1in] \n', latex(imVec(j)));
    fprintf(latexDocTwo,'%s \n', ' ');
end
fclose(latexDocTwo);

termN = termVec(N);
termN = expand(termN);