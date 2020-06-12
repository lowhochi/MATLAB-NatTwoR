% temp_NatTwo_checkRho_Apr23.m
load('DataMain15_Nat_CR_rmW_Apr23.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phiW = d2w_uu - 6*conj(u)*dw_u/(1+u*conj(u)) + 12*conj(u)^2*w/(1+u*conj(u))^2;
rhoTerm_h11 = -d2h11_uconju  + 2*u/(1+u*conj(u))*dh11_u...
    + 2*conj(u)/(1+u*conj(u))*dh11_conju - 4*u*conj(u)/(1+u*conj(u))^2*h11;
rhoTerm_T4 = -6*i*(dT4_u - 2*conj(u)/(1+u*conj(u))*T4)...
    +4*i*(conj(dT4_u)- 2*u/(1+u*conj(u))*conj(T4));

rhoTest1 = i*(phiW-conj(phiW)) + rhoTerm_h11 + rhoTerm_T4 - 2*i*conj(aV);
rhoDiff = rho - rhoTest1;
rhoDiff = complex_simple3(rhoDiff, rhoSet); % rhoDiff = 0 OK.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms G12_3 G23_1 G31_2 G11_2 G11_3 G22_1 G22_3 G33_1 G33_2 real
h11Sub = -mu1*conj(mu1)*G12_3 - mu2*conj(mu2)*G23_1 - mu3*conj(mu3)*G31_2...
    + 1/(2*norm_of_v)*(G11_2*v1*v3 - G11_3*v1*v2 - G22_1*v3*v2...
    + G22_3*v2*v1 + G33_1*v2*v3 - G33_2*v1*v3);
%
T4Sub = -i/(2*norm_of_v)*v1*(mu3*G11_2 - mu2*G11_3 + mu1*G12_3)...
    -i/(2*norm_of_v)*v2*(-mu3*G22_1 + mu2*G23_1 + mu1*G22_3)...
    -i/(2*norm_of_v)*v3*(mu3*G31_2 + mu2*G33_1 -mu1*G33_2);
%
aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);
%
dT4_uSub = complexdiff3(T4Sub,u,0);
dh11_uSub = complexdiff3(h11Sub,u,0);
dh11_conjuSub = complexdiff3(h11Sub,u,1);
d2h11_uconjuSub = complexdiff3(dh11_uSub,u,1);
% 
VarSubSet = [h11, T4, aV, dh11_u, dh11_conju, d2h11_uconju, dT4_u];
rhoSubSet = [h11Sub, T4Sub, aVSub, dh11_uSub, dh11_conjuSub,...
    d2h11_uconjuSub, dT4_uSub];
rhoTest2 = subs(rhoTest1, VarSubSet, rhoSubSet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that rho = I(w) - I(f)
f = -i/2*mu1*(mu3*G11_2 + mu1*G12_3 - mu2*G11_3)...
    -i/2*mu2*(-mu3*G22_1 + mu1*G22_3 + mu2*G23_1)...
    -i/2*mu3*(mu3*G31_2 - mu1*G33_2 + mu2*G33_1);
dfdu = complexdiff3(f,u,0);
d2fdu2 = complexdiff3(dfdu,u,0);
phiF = d2fdu2 - 6*conj(u)/(1+u*conj(u))*dfdu + 12*conj(u)^2/(1+u*conj(u))^2*f;
%
rhoDiff2 = rhoTest2 - i*(phiW-conj(phiW)) + i*(phiF-conj(phiF));
rhoDiff2 = complex_simple3(rhoDiff2, [u,w,dw_u,d2w_uu]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aMuSub0 = -i*(conj(u)^2-1)*(v1*(G23_1+G31_2) + v2*G11_3 - v3*G11_2)...
    -2*i*conj(u)*(-v1*G22_3 + v2*(G12_3+G31_2) + v3*G22_1)...
    -(conj(u)^2+1)*(v1*G33_2 - v2*G33_1 + v3*(G12_3+G23_1));
aMuSub = 1/(2*(1+u*conj(u))^2)*aMuSub0;
daMu_conjuSub = complexdiff3(aMuSub,u,1);
%
VarSubSet2 = [h11, T4, aV, aMu, dh11_u, dh11_conju, d2h11_uconju, dT4_u, daMu_conju];
rhoSubSet2 = [h11Sub, T4Sub, aVSub, aMuSub, dh11_uSub, dh11_conjuSub,...
    d2h11_uconjuSub, dT4_uSub, daMu_conjuSub];
% The Rm-Scalar curvature S. (S1)
STerm_h11 = 2*d2h11_uconju - 4*u/(1+u*conj(u))*dh11_u...
    - 4*conj(u)/(1+u*conj(u))*dh11_conju + 8*h11/(1+u*conj(u))...
    -14*h11/(1+u*conj(u))^2;

STerm_T4 = -3*i*dT4_u + 5*i*conj(dT4_u) + 6*i*conj(u)/(1+u*conj(u))*T4...
    - 10*i*u/(1+u*conj(u))*conj(T4);

STerm_aVaMu = 2*i*aV - 3*i/2*daMu_conju + 3*i/2*conj(daMu_conju)...
    + 3*i*u/(1+u*conj(u))*aMu - 3*i*conj(u)/(1+u*conj(u))*conj(aMu);

S1 = 1/6*rhoTest1 + 3*i/2*(phiW-conj(phiW))...
    + STerm_h11 + STerm_T4 + STerm_aVaMu;
Stest = subs(S1, VarSubSet2, rhoSubSet2);
SminusRho = Stest - 5/3*rhoTest2; 
SminusRho = complex_simple3(SminusRho, [u,w,dw_u,d2w_uu]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%