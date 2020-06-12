% Weyl3_NatTwoR_CR_rmW.m
load('DataWeyl2_NatTwoR_CR_rmW.mat');
load('Data_Chern_in_aV_Oct10.mat');
% work with variable WeylTwo in [m,n,k,ll,W(m,n,k,ll)].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(1) check that W(u6,.,.,.) = 0
indexWithSix = []; % size(indexWithSix) =[65 5];
for j=1:120
    m = countIndex120(j,1);
    n = countIndex120(j,2);
    k = countIndex120(j,3);
    ll = countIndex120(j,4);
    rowTemp = [m,n,k,ll,j]; %j=row index in WeylTwo
    if (m==6)||(n==6)
        indexWithSix =[indexWithSix; rowTemp];
        continue
    end
    if (k==6)||(ll==6)
        indexWithSix =[indexWithSix; rowTemp];
    end
end

WeylSix = sym('WeylSix',[65,5]); 
for myNumber=1:65
    j = indexWithSix(myNumber,5);
    m = indexWithSix(myNumber,1);
    n = indexWithSix(myNumber,2);
    k = indexWithSix(myNumber,3);
    ll = indexWithSix(myNumber,4);
    temp = WeylTwo(j,5);
    WeylSix(myNumber,:) =[m,n,k,ll,temp];
end
W1256 = WeylSix(5,5);
W1526 = WeylSix(17,5);
W1556 = WeylSix(20,5);
W1625 = WeylSix(24,5);
W2556 = WeylSix(43,5);

testW6 = [W1256, W1526, W1556];
variableSetw6 = symvar([W1256, W1526, W1556]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variable_NatTwoR_Chern2_CR_rmW

variableSetG1 = [theta, aV, daV_u, daV_conju, d2aV_uu];
variableSetG2 = [daV_vnormv, daV_conjMu, dtheta_vnormv];
variableSetG3 = [d2aV_uMu, d2aV_uvnormv, d2aV_conjuconjMu, d3aV_uuMu];

variableSet0 = [aM, bV, h11, T4];
subSet0 = [aMSub, bVSub, T4Sub, h11SubTwo];
% substitute variable0 BEFORE [theta, aV, daV_u, daV_conju, d2aV_uu].

thetaSub = G12_3 + G23_1 + G31_2;
dtheta_vnormvSub = dG12_3_vnormv + dG23_1_vnormv + dG31_2_vnormv;

aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);

daVVec = df_NatTwo_MuSet_CR_rmW(aVSub, CVarG1, dCVarG1);
daV_conjMuSub = daVVec(2);
daV_vnormvSub = daVVec(3);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);
d2aV_uVec = df_NatTwo_MuSet_CR_rmW(daV_uSub, CVarG1, dCVarG1);
d2aV_uMuSub = d2aV_uVec(1);
d2aV_uvnormvSub = d2aV_uVec(3);
d2aV_uuSub = d2aV_uVec(4);
d2aV_conjuVec = df_NatTwo_MuSet_CR_rmW(daV_conjuSub, CVarG1, dCVarG1);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);
d3aV_uuVec = df_NatTwo_MuSet_CR_rmW(d2aV_uuSub, CVarG1, dCVarG1);
d3aV_uuMuSub = d3aV_uuVec(1);

subSetG1 = [thetaSub, aVSub, daV_uSub, daV_conjuSub, d2aV_uuSub];
subSetG2 = [daV_vnormvSub, daV_conjMuSub, dtheta_vnormvSub];
subSetG3 = [d2aV_uMuSub, d2aV_uvnormvSub, d2aV_conjuconjMuSub, d3aV_uuMuSub];

for j=1:3
    temp = testW6(j);
    temp = subs(temp, variableSetG3, subSetG3);
    temp = subs(temp, variableSetG2, subSetG2);
    temp = subs(temp, variableSet0, subSet0);
    temp = subs(temp, variableSetG1, subSetG1);
    % Second substitution
    temp = subs(temp, variableSetMu, subSetMu);
    temp = subs(temp, variableSetConjMu, subSetConjMu);
    temp = subs(temp, variableSetVnormv, subSetVnormv);
    % Bianchi Identities
    temp = subs(temp, variableSetBianchi2, subSetBianchi2);
    temp = subs(temp, variableSetBianchi1, subSetBianchi1);
    temp = complex_simple3(temp,[u]);
    testW6(j)=temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars m n k ll j myNumber rowTemp temp