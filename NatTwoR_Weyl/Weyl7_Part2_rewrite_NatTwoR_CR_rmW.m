% Weyl7_Part2_rewrite_NatTwoR_CR_rmW.m
load('DataWeyl7new_NatTwoR_CR_rmW.mat');
Wf1214 = Weylf(3,5);
Wf1235 = Weylf(11,5);
Wf1535 = Weylf(50,5);
Wf1545 = Weylf(52,5);
group2Vec = [Wf1214, Wf1235, Wf1535, Wf1545];

MVar72 = symvar(group2Vec);
Ric12Sub = d2aV_conjuconjMu...
    -(1+u*conj(u))^2*(daV_u*daV_conju)...
    +(1+u*conj(u))^2*conj(daV_u)*conj(daV_conju)...
    -2*daV_vnormv +2*conj(daV_vnormv) +4*i*dtheta_vnormv...
    +4*i*theta*(aV+conj(aV)) +2*conj(aV)^2 - 2*aV^2;

for j=1:4
   temp01 = group2Vec(j);
   temp01 = subs(temp01, conj(d2aV_conjuconjMu),Ric12Sub);
   temp01 = complex_simple3(temp01, MVar72);
   group2Vec(j) = temp01;
end

variable_NatTwoR_Weyl7Part2_CR_rmW % variable file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVarW72 = [u, G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3,...
    G33_1, G33_2, dG11_2_Mu, dG11_3_Mu, dG12_3_Mu, dG22_1_Mu, dG22_3_Mu,...
    dG23_1_Mu, dG31_2_Mu, dG33_1_Mu, dG33_2_Mu,...
    dG11_2_conjMu, dG11_3_conjMu, dG12_3_conjMu, dG22_1_conjMu,...
    dG22_3_conjMu, dG23_1_conjMu, dG31_2_conjMu, dG33_1_conjMu,...
    dG33_2_conjMu, dG11_2_vnormv, dG11_3_vnormv, dG12_3_vnormv,...
    dG22_1_vnormv, dG22_3_vnormv, dG23_1_vnormv, dG31_2_vnormv,...
    dG33_1_vnormv, dG33_2_vnormv];

gDerivDict.u = [0; 0; 0; 1; 0];
gDerivDict.G12_3 = [dG12_3_Mu; dG12_3_conjMu; dG12_3_vnormv; 0; 0];
gDerivDict.G23_1 = [dG23_1_Mu; dG23_1_conjMu; dG23_1_vnormv; 0; 0];
gDerivDict.G31_2 = [dG31_2_Mu; dG31_2_conjMu; dG31_2_vnormv; 0; 0];
gDerivDict.G11_2 = [dG11_2_Mu; dG11_2_conjMu; dG11_2_vnormv; 0; 0];
gDerivDict.G11_3 = [dG11_3_Mu; dG11_3_conjMu; dG11_3_vnormv; 0; 0];
gDerivDict.G22_1 = [dG22_1_Mu; dG22_1_conjMu; dG22_1_vnormv; 0; 0];
gDerivDict.G22_3 = [dG22_3_Mu; dG22_3_conjMu; dG22_3_vnormv; 0; 0];
gDerivDict.G33_1 = [dG33_1_Mu; dG33_1_conjMu; dG33_1_vnormv; 0; 0];
gDerivDict.G33_2 = [dG33_2_Mu; dG33_2_conjMu; dG33_2_vnormv; 0; 0];

gDerivDict.dG12_3_Mu = [d2G12_3_MuMu; d2G12_3_MuconjMu; d2G12_3_Muvnormv;
    d2G12_3_Muu; 0];
gDerivDict.dG12_3_conjMu = [d2G12_3_conjMuMu; d2G12_3_conjMuconjMu;
    d2G12_3_conjMuvnormv; 0; d2G12_3_conjMuconju];
gDerivDict.dG12_3_vnormv = [d2G12_3_vnormvMu; d2G12_3_vnormvconjMu;
    d2G12_3_vnormvvnormv; d2G12_3_vnormvu; d2G12_3_vnormvconju];

gDerivDict.dG23_1_Mu = [d2G23_1_MuMu; d2G23_1_MuconjMu; d2G23_1_Muvnormv;
    d2G23_1_Muu; 0];
gDerivDict.dG23_1_conjMu = [d2G23_1_conjMuMu; d2G23_1_conjMuconjMu;
    d2G23_1_conjMuvnormv; 0; d2G23_1_conjMuconju];
gDerivDict.dG23_1_vnormv = [d2G23_1_vnormvMu; d2G23_1_vnormvconjMu;
    d2G23_1_vnormvvnormv; d2G23_1_vnormvu; d2G23_1_vnormvconju];

gDerivDict.dG31_2_Mu = [d2G31_2_MuMu; d2G31_2_MuconjMu; d2G31_2_Muvnormv;
    d2G31_2_Muu; 0];
gDerivDict.dG31_2_conjMu = [d2G31_2_conjMuMu; d2G31_2_conjMuconjMu;
    d2G31_2_conjMuvnormv; 0; d2G31_2_conjMuconju];
gDerivDict.dG31_2_vnormv = [d2G31_2_vnormvMu; d2G31_2_vnormvconjMu;
    d2G31_2_vnormvvnormv; d2G31_2_vnormvu; d2G31_2_vnormvconju];

gDerivDict.dG11_2_Mu = [d2G11_2_MuMu; d2G11_2_MuconjMu; d2G11_2_Muvnormv;
    d2G11_2_Muu; 0];
gDerivDict.dG11_2_conjMu = [d2G11_2_conjMuMu; d2G11_2_conjMuconjMu;
    d2G11_2_conjMuvnormv; 0; d2G11_2_conjMuconju];
gDerivDict.dG11_2_vnormv = [d2G11_2_vnormvMu; d2G11_2_vnormvconjMu;
    d2G11_2_vnormvvnormv; d2G11_2_vnormvu; d2G11_2_vnormvconju];

gDerivDict.dG11_3_Mu = [d2G11_3_MuMu; d2G11_3_MuconjMu; d2G11_3_Muvnormv;
    d2G11_3_Muu; 0];
gDerivDict.dG11_3_conjMu = [d2G11_3_conjMuMu; d2G11_3_conjMuconjMu;
    d2G11_3_conjMuvnormv; 0; d2G11_3_conjMuconju];
gDerivDict.dG11_3_vnormv = [d2G11_3_vnormvMu; d2G11_3_vnormvconjMu;
    d2G11_3_vnormvvnormv; d2G11_3_vnormvu; d2G11_3_vnormvconju];

gDerivDict.dG22_1_Mu = [d2G22_1_MuMu; d2G22_1_MuconjMu; d2G22_1_Muvnormv;
    d2G22_1_Muu; 0];
gDerivDict.dG22_1_conjMu = [d2G22_1_conjMuMu; d2G22_1_conjMuconjMu;
    d2G22_1_conjMuvnormv; 0; d2G22_1_conjMuconju];
gDerivDict.dG22_1_vnormv = [d2G22_1_vnormvMu; d2G22_1_vnormvconjMu;
    d2G22_1_vnormvvnormv; d2G22_1_vnormvu; d2G22_1_vnormvconju];

gDerivDict.dG22_3_Mu = [d2G22_3_MuMu; d2G22_3_MuconjMu; d2G22_3_Muvnormv;
    d2G22_3_Muu; 0];
gDerivDict.dG22_3_conjMu = [d2G22_3_conjMuMu; d2G22_3_conjMuconjMu;
    d2G22_3_conjMuvnormv; 0; d2G22_3_conjMuconju];
gDerivDict.dG22_3_vnormv = [d2G22_3_vnormvMu; d2G22_3_vnormvconjMu;
    d2G22_3_vnormvvnormv; d2G22_3_vnormvu; d2G22_3_vnormvconju];

gDerivDict.dG33_1_Mu = [d2G33_1_MuMu; d2G33_1_MuconjMu; d2G33_1_Muvnormv;
    d2G33_1_Muu; 0];
gDerivDict.dG33_1_conjMu = [d2G33_1_conjMuMu; d2G33_1_conjMuconjMu;
    d2G33_1_conjMuvnormv; 0; d2G33_1_conjMuconju];
gDerivDict.dG33_1_vnormv = [d2G33_1_vnormvMu; d2G33_1_vnormvconjMu;
    d2G33_1_vnormvvnormv; d2G33_1_vnormvu; d2G33_1_vnormvconju];

gDerivDict.dG33_2_Mu = [d2G33_2_MuMu; d2G33_2_MuconjMu; d2G33_2_Muvnormv;
    d2G33_2_Muu; 0];
gDerivDict.dG33_2_conjMu = [d2G33_2_conjMuMu; d2G33_2_conjMuconjMu;
    d2G33_2_conjMuvnormv; 0; d2G33_2_conjMuconju];
gDerivDict.dG33_2_vnormv = [d2G33_2_vnormvMu; d2G33_2_vnormvconjMu;
    d2G33_2_vnormvvnormv; d2G33_2_vnormvu; d2G33_2_vnormvconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSetG1 = [theta, aV, daV_u, daV_conju, d2aV_uu];

variableSetG2=[daV_Mu, daV_conjMu, daV_vnormv, ...
    d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv,...
    d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv,...
    d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv, dtheta_Mu, dtheta_vnormv];

aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);

thetaSub = G12_3 + G23_1 + G31_2;

daVVec = df_main4_CR_funWv2(aVSub, CVarW72, gDerivDict, gamma);
daV_MuSub = daVVec(1);
daV_conjMuSub = daVVec(2);
daV_vnormvSub = daVVec(3);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);

d2aV_uVec = df_main4_CR_funWv2(daV_uSub, CVarW72, gDerivDict, gamma);
d2aV_uMuSub = d2aV_uVec(1);
d2aV_uconjMuSub = d2aV_uVec(2);
d2aV_uvnormvSub = d2aV_uVec(3);
d2aV_uuSub = d2aV_uVec(4);

d2aV_conjuVec = df_main4_CR_funWv2(daV_conjuSub, CVarW72, gDerivDict, gamma);
d2aV_conjuMuSub = d2aV_conjuVec(1);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);
d2aV_conjuvnormvSub = d2aV_conjuVec(3);

d3aV_uuVec = df_main4_CR_funWv2(d2aV_uuSub, CVarW72, gDerivDict, gamma);
d3aV_uuMuSub = d3aV_uuVec(1);
d3aV_uuconjMuSub = d3aV_uuVec(2);
d3aV_uuvnormvSub = d3aV_uuVec(3);

% % break point % %
variableSetG3 = [d2aV_MuMu, d2aV_MuconjMu, d2aV_Muvnormv,...
    d2aV_conjMuconjMu, d2aV_conjMuvnormv, d2aV_vnormvvnormv,...
    d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv, d3aV_uconjMuconjMu,...
    d3aV_uconjMuvnormv, d3aV_uvnormvvnormv,...
    d3aV_conjuMuMu, d3aV_conjuMuconjMu, d3aV_conjuMuvnormv,...
    d3aV_conjuconjMuconjMu, d3aV_conjuconjMuvnormv];

d2aV_MuVec = df_main4_CR_funWv2(daV_MuSub, CVarW72, gDerivDict, gamma);
d2aV_conjMuVec = df_main4_CR_funWv2(daV_conjMuSub, CVarW72, gDerivDict, gamma);
d2aV_vnormvVec = df_main4_CR_funWv2(daV_vnormvSub, CVarW72, gDerivDict, gamma);
d2aV_MuMuSub = d2aV_MuVec(1);
d2aV_MuconjMuSub = d2aV_MuVec(2);
d2aV_MuvnormvSub = d2aV_MuVec(3);
d2aV_conjMuconjMuSub = d2aV_conjMuVec(2);
d2aV_conjMuvnormvSub = d2aV_conjMuVec(3);
d2aV_vnormvvnormvSub = d2aV_vnormvVec(3);

d3aV_uMuVec = df_main4_CR_funWv2(d2aV_uMuSub, CVarW72, gDerivDict, gamma);
d3aV_uconjMuVec = df_main4_CR_funWv2(d2aV_uconjMuSub, CVarW72, gDerivDict, gamma);
d3aV_uvnormvVec = df_main4_CR_funWv2(d2aV_uvnormvSub, CVarW72, gDerivDict, gamma);
d3aV_uMuMuSub = d3aV_uMuVec(1);
d3aV_uMuconjMuSub = d3aV_uMuVec(2);
d3aV_uMuvnormvSub = d3aV_uMuVec(3);
d3aV_uconjMuconjMuSub = d3aV_uconjMuVec(2); %add
d3aV_uconjMuvnormvSub = d3aV_uconjMuVec(3);
d3aV_uvnormvvnormvSub = d3aV_uvnormvVec(3);

d3aV_conjuMuVec = df_main4_CR_funWv2(d2aV_conjuMuSub, CVarW72, gDerivDict, gamma);
d3aV_conjuconjMuVec = df_main4_CR_funWv2(d2aV_conjuconjMuSub, CVarW72, gDerivDict, gamma);
d3aV_conjuMuMuSub = d3aV_conjuMuVec(1);
d3aV_conjuMuconjMuSub = d3aV_conjuMuVec(2);
d3aV_conjuMuvnormvSub = d3aV_conjuMuVec(3);
d3aV_conjuconjMuconjMuSub = d3aV_conjuconjMuVec(2);
d3aV_conjuconjMuvnormvSub = d3aV_conjuconjMuVec(3);

% % break point % %
variableSetG4 = [d4aV_uuMuMu, d4aV_uuMuconjMu, d4aV_uuMuvnormv, ...
    d4aV_uuvnormvvnormv, d4aV_uuconjMuconjMu, d4aV_uuconjMuvnormv,...
    d2theta_MuMu, d2theta_MuconjMu, d2theta_Muvnormv, d2theta_vnormvvnormv];

d4aV_uuMuVec = df_main4_CR_funWv2(d3aV_uuMuSub, CVarW72, gDerivDict, gamma);
d4aV_uuvnormvVec = df_main4_CR_funWv2(d3aV_uuvnormvSub, CVarW72, gDerivDict, gamma);
d4aV_uuconjMuVec = df_main4_CR_funWv2(d3aV_uuconjMuSub, CVarW72, gDerivDict, gamma);

d4aV_uuMuMuSub = d4aV_uuMuVec(1);
d4aV_uuMuconjMuSub = d4aV_uuMuVec(2);
d4aV_uuMuvnormvSub = d4aV_uuMuVec(3);
d4aV_uuconjMuconjMuSub = d4aV_uuconjMuVec(2);
d4aV_uuconjMuvnormvSub = d4aV_uuconjMuVec(3);
d4aV_uuvnormvvnormvSub = d4aV_uuvnormvVec(3);

dthetaVec = df_main4_CR_funWv2(thetaSub, CVarW72, gDerivDict, gamma);
dtheta_MuSub = dthetaVec(1);
dtheta_vnormvSub = dthetaVec(3);

d2theta_MuVec = df_main4_CR_funWv2(dtheta_MuSub, CVarW72, gDerivDict, gamma);
d2theta_vnormvVec = df_main4_CR_funWv2(dtheta_vnormvSub, CVarW72, gDerivDict, gamma);
d2theta_MuMuSub = d2theta_MuVec(1);
d2theta_MuconjMuSub = d2theta_MuVec(2);
d2theta_MuvnormvSub = d2theta_MuVec(3);
d2theta_vnormvvnormvSub = d2theta_vnormvVec(3);

for j=1:4
    temp = group2Vec(j);
    for k=1:length(variableSetG4)
        myChar = char(variableSetG4(k));
        subChar = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar, ',', subChar, ');']);
    end
    for k=1:length(variableSetG3)
        myChar = char(variableSetG3(k));
        subChar = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar, ',', subChar, ');']);
    end
    for k=1:length(variableSetG2)
        myChar = char(variableSetG2(k));
        subChar = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar, ',', subChar, ');']);
    end
    for k=1:length(variableSetG1)
        myChar = char(variableSetG1(k));
        subChar = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar, ',', subChar, ');']);
    end
    group2Vec(j) = temp; 
end

% Second substitution
for j=1:4
    temp = group2Vec(j);
    temp = subs(temp, d2G12_3N, d2G12_3Vec);
    temp = subs(temp, d2G23_1N, d2G23_1Vec);
    temp = subs(temp, d2G31_2N, d2G31_2Vec);
    temp = subs(temp, d2G11_2N, d2G11_2Vec);    
    temp = subs(temp, d2G11_3N, d2G11_3Vec);
    temp = subs(temp, d2G22_1N, d2G22_1Vec);
    temp = subs(temp, d2G22_3N, d2G22_3Vec);
    temp = subs(temp, d2G33_1N, d2G33_1Vec);
    temp = subs(temp, d2G33_2N, d2G33_2Vec);
    temp = subs(temp, variableSetMu, subSetMu);
    temp = subs(temp, variableSetConjMu, subSetConjMu);
    temp = subs(temp, variableSetVnormv, subSetVnormv);
    temp = subs(temp, variableSetBianchi2, subSetBianchi2);
    temp = subs(temp, variableSetBianchi1, subSetBianchi1);
    group2Vec(j) = complex_simple3(temp, u);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



