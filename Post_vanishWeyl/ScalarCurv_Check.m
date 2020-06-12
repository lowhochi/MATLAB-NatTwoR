load('DataWeyl2_NatTwoR_CR_rmW.mat');
variableS = symvar(S); 
Stemp = subs(S, d2h11_uconju, d2h11_uconjuSub);
Stemp = subs(Stemp, dh11_u, dh11_uSub);
Stemp = subs(Stemp, dh11_conju, dh11_conjuSub);
Stemp = subs(Stemp, daM_conju, daM_conjuSub);
Stemp = subs(Stemp, dT4_u, dT4_uSub);
Stemp = subs(Stemp, rho, rhoSub);
Stemp = subs(Stemp, h11, h11Sub);
Stemp = subs(Stemp, T4, T4Sub);
Stemp = subs(Stemp, aM, aMSub);
% symvar(Stemp) = [aV, d2aV_uconju, d2w_uu, d3aV_uconjuconju,...
%   d3aV_uuconju, d4aV_uuconjuconju, daV_conju, daV_u,...
%   dw_u, theta, u, w]
Stemp = subs(Stemp, d4aV_uuconjuconju, d4aV_uuconjuconjuR);
Stemp = subs(Stemp, d3aV_uconjuconju, d3aV_uconjuconjuR);
Stemp = subs(Stemp, d3aV_uuconju, d3aV_uuconjuR);
Stemp = subs(Stemp, d2aV_uconju, d2aV_uconjuR);
Stemp = complex_simple3(Stemp, MVarWeylTwo);
test = Stemp - 5/3*rhoSub;
test = complex_simple3(test, MVarWeylTwo);