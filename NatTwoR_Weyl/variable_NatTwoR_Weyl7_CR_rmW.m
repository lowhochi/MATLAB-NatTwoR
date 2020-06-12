% variable_NatTwoR_Weyl7_CR_rmW.m
d2aV_vnormvMuR = d2aV_Muvnormv + 2*T4Sub*daV_vnormv + aV*daV_Mu + bVSub*daV_conjMu;
d2aV_vnormvconjMuR = d2aV_conjMuvnormv...
    + 2*daV_vnormv*conj(T4Sub) + daV_conjMu*conj(aV) + daV_Mu*conj(bVSub);
d2aV_conjMuMuR = d2aV_MuconjMu + aMSub*daV_Mu...
    + daV_vnormv*h11SubTwo*2i - daV_conjMu*conj(aMSub);
d3aV_uconjMuMuR = d3aV_uMuconjMu + aMSub*d2aV_uMu...
    + d2aV_uvnormv*h11SubTwo*2i - d2aV_uconjMu*conj(aMSub);
d3aV_uconjMuconjuR = 2*d2aV_uvnormv + d3aV_uconjuconjMuR...
    +(2*d2aV_uconjMu*u)/(u*conj(u)+1);
d3aV_uvnormvMuR = d3aV_uMuvnormv +2*T4Sub*d2aV_uvnormv...
    +aV*d2aV_uMu +bVSub*d2aV_uconjMu;
d3aV_uvnormvconjMuR =d3aV_uconjMuvnormv +2*d2aV_uvnormv*conj(T4Sub)...
    +d2aV_uconjMu*conj(aV) +d2aV_uMu*conj(bVSub);
d3aV_uvnormvconjuR = d3aV_uconjuvnormvR - d2aV_uMu/(u*conj(u) + 1)^2;
d3aV_conjuMuuR = d3aV_uconjuMuR + 2*d2aV_conjuvnormv...
    + (2*d2aV_conjuMu*conj(u))/(u*conj(u) + 1);
d3aV_conjuconjMuMuR = d3aV_conjuMuconjMu +aMSub*d2aV_conjuMu ...
    + d2aV_conjuvnormv*h11SubTwo*2i - d2aV_conjuconjMu*conj(aMSub);
d3aV_conjuconjMuconjuR = 2*d2aV_conjuvnormv + d3aV_conjuconjuconjMuR...
    + (2*d2aV_conjuconjMu*u)/(u*conj(u) + 1);
d2theta_vnormvMuR = d2theta_Muvnormv...
    + 2*T4Sub*dtheta_vnormv + aV*dtheta_Mu + bVSub*dtheta_conjMu;
d2theta_vnormvconjMuR = conj(d2theta_Muvnormv)...
    + 2*conj(T4Sub)*conj(dtheta_vnormv)+conj(aV)*conj(dtheta_Mu)...
    + conj(bVSub)*conj(dtheta_conjMu);

% define d4aV_uuxx
syms d4aV_uuuu d4aV_uuuMu d4aV_uuuconjMu d4aV_uuuvnormv
syms d4aV_uuMuMu d4aV_uuMuconjMu d4aV_uuconjMuconjMu
syms d4aV_uuMuvnormv d4aV_uuconjMuvnormv d4aV_uuvnormvvnormv
d4aV_uuuuR = 12*conj(u)^3/(1+u*conj(u))^3*daV_u...
    -6*conj(u)/(1+u*conj(u))*d3aV_uuuR;
% d3aV_uuuR=-(6*d2aV_uu*conj(u))/(u*conj(u)+1)-(6*daV_u*conj(u)^2)/(u*conj(u) + 1)^2;
d4aV_uuuMuR= -(6*d3aV_uuMu*conj(u))/(u*conj(u)+1)...
    -(6*d2aV_uMu*conj(u)^2)/(u*conj(u)+1)^2;
d4aV_uuuconjMuR = -(6*d3aV_uuconjMu*conj(u))/(u*conj(u)+1)...
    -(6*d2aV_uconjMu*conj(u)^2)/(u*conj(u)+1)^2;
d4aV_uuuvnormvR = -(6*d3aV_uuvnormv*conj(u))/(u*conj(u)+1)...
    -(6*d2aV_uvnormv*conj(u)^2)/(u*conj(u)+1)^2;

d3aV_uuRow = [d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv];
lie_Mu_conjMuSub = [aMSub; -conj(aMSub); 2*i*h11SubTwo];
lie_Mu_vnormvSub = [aV; bVSub; 2*T4Sub]; 
lie_conjMu_vnormvSub = [conj(bVSub); conj(aV); 2*conj(T4Sub)];
d4aV_uuconjMuMu = d4aV_uuMuconjMu+ d3aV_uuRow*lie_Mu_conjMuSub;
d4aV_uuvnormvMu = d4aV_uuMuvnormv+ d3aV_uuRow*lie_Mu_vnormvSub;
d4aV_uuvnormvconjMu = d4aV_uuconjMuvnormv +d3aV_uuRow*lie_conjMu_vnormvSub;
d4aV_uuMuu = d4aV_uuuMuR -d3aV_uuRow*lie_Mu_ddu;
d4aV_uuconjMuconju = d4aV_uuconjuconjMuR - d3aV_uuRow*lie_conjMu_ddconju;
d4aV_uuvnormvu = d4aV_uuuvnormvR -d3aV_uuRow*lie_vnormv_ddu;
d4aV_uuvnormvconju = d4aV_uuconjuvnormvR -d3aV_uuRow*lie_vnormv_ddconju;


CVarW7 = [u, aV, daV_u, d2aV_uu, daV_conju, theta,...
    daV_Mu, daV_conjMu, dtheta_Mu, dtheta_conjMu,...
    d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv, d2aV_conjuMu, d2aV_conjuconjMu,...
    d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv];

derivativeDict.u =[0; 0; 0; 1; 0];

derivativeDict.aV =[daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju];

derivativeDict.daV_u = [d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv;
    d2aV_uu; d2aV_uconjuR];

derivativeDict.d2aV_uu = [d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv;
    d3aV_uuuR; d3aV_uuconjuR];

derivativeDict.daV_conju = [d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;
    d2aV_uconjuR; d2aV_conjuconjuR];

derivativeDict.theta = [dtheta_Mu; dtheta_conjMu; dtheta_vnormv; 0; 0];

derivativeDict.d2aV_uMu = [d3aV_uMuMu; d3aV_uMuconjMu; d3aV_uMuvnormv;
    d3aV_uMuu; d3aV_uconjuMuR];

derivativeDict.d2aV_uconjMu = [d3aV_uconjMuMuR; d3aV_uconjMuconjMu;
    d3aV_uconjMuvnormv; d3aV_uuconjMu; d3aV_uconjMuconjuR];

derivativeDict.d2aV_uvnormv = [d3aV_uvnormvMuR; d3aV_uvnormvconjMuR;
    d3aV_uvnormvvnormv; d3aV_uvnormvu; d3aV_uvnormvconjuR];

derivativeDict.d3aV_uuconjMu = [d4aV_uuconjMuMu; d4aV_uuconjMuconjMu;
    d4aV_uuconjMuvnormv; d4aV_uuuconjMuR; d4aV_uuconjMuconju];

derivativeDict.d3aV_uuMu = [d4aV_uuMuMu; d4aV_uuMuconjMu; d4aV_uuMuvnormv;
    d4aV_uuMuu; d4aV_uuconjuMuR];

derivativeDict.d3aV_uuvnormv = [d4aV_uuvnormvMu; d4aV_uuvnormvconjMu;
    d4aV_uuvnormvvnormv; d4aV_uuvnormvu; d4aV_uuvnormvconju];

derivativeDict.d2aV_conjuMu = [d3aV_conjuMuMu; d3aV_conjuMuconjMu;
    d3aV_conjuMuvnormv; d3aV_conjuMuuR; d3aV_conjuconjuMuR];

derivativeDict.d2aV_conjuconjMu = [d3aV_conjuconjMuMuR;
    d3aV_conjuconjMuconjMu; d3aV_conjuconjMuvnormv;
    d3aV_uconjuconjMuR; d3aV_conjuconjMuconjuR];

derivativeDict.daV_Mu = [d2aV_MuMu; d2aV_MuconjMu; d2aV_Muvnormv;
    d2aV_Muu; d2aV_conjuMu];

derivativeDict.daV_conjMu = [d2aV_conjMuMuR; d2aV_conjMuconjMu;
    d2aV_conjMuvnormv; d2aV_uconjMu; d2aV_conjMuconju];

derivativeDict.dtheta_Mu = [d2theta_MuMu; d2theta_MuconjMu; d2theta_Muvnormv;
    d2theta_Muu; 0];

derivativeDict.dtheta_conjMu = [d2theta_conjMuMu; d2theta_conjMuconjMu;
    d2theta_conjMuvnormv; 0; d2theta_conjMuconju];

