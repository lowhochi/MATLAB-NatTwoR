load('DataTest_Chern2_NatTwoR_CR_rmW.mat');
C1111 = Chern_in_aV(1,1,1,1);

wSet = [w, dw_u, dw_conjMu, dw_vnormv, d2w_uu, d2w_uconjMu,...
    d2w_conjMuconjMu];
 
C1111 = subs(C1111,wSet,zeros(1,7));

testPart1 = Y0^2*(2*d2theta_MuconjMu + 2*i*d2aV_MuconjMu...
    -i*conj(d2aV_MuconjMu));
testPart2 = -Y0^4*(2*theta+i*aV-i*conj(aV))...
    *(4*i*dtheta_vnormv + d2aV_conjuconjMu - 2*daV_vnormv);
testPart3 = -Y0^4*(4*dtheta_Mu*daV_u + 2*dtheta_conjMu*conj(daV_u));
testPart4 = -i*Y0^4*((2*daV_Mu - conj(daV_conjMu))*daV_u...
    - (2*conj(daV_Mu) - daV_conjMu)*conj(daV_u));
testPart5 = Y0^6*(2*theta+i*aV-i*conj(aV))...
    *(daV_u*daV_conju + daV_u*conj(daV_u));
testPart6 = Y0^4*(-16/3*theta^3 - 8*i/3*theta^2*(5*aV-2*conj(aV))...
    +4/3*theta*(7*aV-conj(aV))*(aV-conj(aV)) + 2*i*aV*(aV-conj(aV))^2);

test1111 = testPart1 + testPart2 + testPart3 + testPart4...
    + testPart5 + testPart6;


checkdifference = test1111- C1111;
checkdifference = subs(checkdifference, Y0, u*conj(u)+1);
checkdifference = complex_simple3(checkdifference, MVarChern1S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% testTwo = Y0^2*(2*d2theta_MuconjMu + i*d2aV_MuconjMu)...
%     - Y0^4*(2*theta+i*aV-i*conj(aV))...
%         *(4*i*dtheta_vnormv + d2aV_conjuconjMu - 2*daV_vnormv)...
%     -2*Y0^4*dtheta_Mu*daV_u...
%     + Y0^6*(2*theta+i*aV-i*conj(aV))*(daV_u*daV_conju)...
%     + Y0^4*(-8*i*theta^2*aV+8*theta*aV*(aV-conj(aV))+2*i*aV*(aV-conj(aV))^2);

testThree = 2*Y0^2*d2theta_MuconjMu + i*Y0^2*d2aV_MuconjMu...
    -2*i*Y0^4*(2*theta+i*aV-i*conj(aV))*dtheta_vnormv - 2*Y0^4*dtheta_Mu*daV_u;

test4 =   2*Y0^2*d2theta_MuconjMu...
    -2*i*Y0^4*(2*theta+i*aV-i*conj(aV))*dtheta_vnormv...
    - 2*Y0^4*dtheta_Mu*daV_u;

test5 = -Y0^4*(2*theta+i*aV-i*conj(aV))...
    *(2*i*dtheta_vnormv + d2aV_conjuconjMu - 2*daV_vnormv)...
    +Y0^6*(2*theta+i*aV-i*conj(aV))*(daV_u*daV_conju)...
    +Y0^4*(-16/3*theta^3 - 8*i/3*theta^2*(5*aV-2*conj(aV))...
        +4/3*theta*(7*aV-conj(aV))*(aV-conj(aV)) + 2*i*aV*(aV-conj(aV))^2);

variable_NatTwoR_Chern2_CR_rmW

variableSet0 = [aM, bV, h11, T4];
subSet0 = [aMSubTwo, bVSubTwo, T4SubTwo, h11SubTwo];
% substitute variable0 BEFORE [theta, aV, daV_u, daV_conju, d2aV_uu].
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

daVVec = df_NatTwo_MuSet_CR_rmW(aVSub, CVarG1, dCVarG1);
daV_MuSub = daVVec(1);
% symvar(daV_MuSub)=[dG11_2_Mu, dG11_3_Mu, dG12_3_Mu, dG22_1_Mu,...
%   dG22_3_Mu, dG23_1_Mu, dG31_2_Mu, dG33_1_Mu, dG33_2_Mu, u];
daV_conjMuSub = daVVec(2);
daV_vnormvSub = daVVec(3);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);

d2aV_uVec = df_NatTwo_MuSet_CR_rmW(daV_uSub, CVarG1, dCVarG1);
d2aV_uMuSub = d2aV_uVec(1);
d2aV_uconjMuSub = d2aV_uVec(2);
d2aV_uuSub = d2aV_uVec(4);

d2aV_conjuVec = df_NatTwo_MuSet_CR_rmW(daV_conjuSub, CVarG1, dCVarG1);
d2aV_conjuMuSub = d2aV_conjuVec(1);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);

d2aV_MuVec = df_NatTwo_MuSet_CR_rmW(daV_MuSub, CVarG1, dCVarG1);
d2aV_MuconjMuSub = d2aV_MuVec(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Substitution on Chern_in_aV
variableSetG1 = [theta, aV, daV_u, daV_conju, d2aV_uu];
subSetG1 = [thetaSub, aVSub, daV_uSub, daV_conjuSub, d2aV_uuSub]; %latest

variableSetG2 = [daV_Mu, daV_conjMu, daV_vnormv, dtheta_Mu, dtheta_vnormv,...
    d2theta_MuconjMu];
subSetG2 = [daV_MuSub, daV_conjMuSub, daV_vnormvSub, dtheta_MuSub,...
    dtheta_vnormvSub, d2theta_MuconjMuSub]; %middle

variableSetG3 = [d2aV_uMu, d2aV_uconjMu, d2aV_conjuMu, d2aV_conjuconjMu,...
    d2aV_MuconjMu];
subSetG3 =  [d2aV_uMuSub, d2aV_uconjMuSub, d2aV_conjuMuSub, d2aV_conjuconjMuSub,...
    d2aV_MuconjMuSub]; %first

test5  = subs(test5 , variableSetG3, subSetG3);
test5  = subs(test5 , variableSetG2, subSetG2);
test5  = subs(test5 , variableSet0, subSet0);
test5  = subs(test5 , variableSetG1, subSetG1);

test5  = subs(test5 , d2G12_3N, d2G12_3Vec);
test5  = subs(test5 , d2G23_1N, d2G23_1Vec);
test5  = subs(test5 , d2G31_2N, d2G31_2Vec);
test5  = subs(test5 , d2G11_2N, d2G11_2Vec);
test5  = subs(test5 , d2G11_3N, d2G11_3Vec);
test5  = subs(test5 , d2G22_1N, d2G22_1Vec);
test5  = subs(test5 , d2G22_3N, d2G22_3Vec);
test5  = subs(test5 , d2G33_1N, d2G33_1Vec);
test5  = subs(test5 , d2G33_2N, d2G33_2Vec);
test5  = subs(test5 , variableSetMu, subSetMu);
test5  = subs(test5 , variableSetConjMu, subSetConjMu);
test5  = subs(test5 , variableSetVnormv, subSetVnormv);
test5  = subs(test5 , Y0, u*conj(u)+1);
% Bianchi Identities
test5  = subs(test5 , variableSetBianchi2, subSetBianchi2);
test5  = subs(test5 , variableSetBianchi1, subSetBianchi1);
test5  = complex_simple3(test5 ,[u]);
test5  = subs(test5 , variableSetBianchi2, subSetBianchi2);
test5  = subs(test5 , variableSetBianchi1, subSetBianchi1);
test5  = complex_simple3(test5 ,[u]);

variabletest4= symvar(test5);

Gset = [G11_2, G11_3, G12_3, G22_1, G22_3, G23_1, G31_2, G33_1, G33_2,...
    d2G11_2_by11, d2G11_2_by12, d2G11_2_by22, d2G11_3_by11, d2G11_3_by12,...
    d2G11_3_by13, d2G11_3_by22, d2G11_3_by23, d2G11_3_by33, d2G12_3_by11,...
    d2G12_3_by12, d2G12_3_by13, d2G12_3_by22, d2G12_3_by23, d2G12_3_by33,...
    d2G22_1_by11, d2G22_1_by12, d2G22_1_by13, d2G22_1_by22, d2G22_1_by23,...
    d2G22_1_by33, d2G22_3_by22, d2G22_3_by23, d2G22_3_by33, d2G23_1_by11,...
    d2G23_1_by12, d2G23_1_by13, d2G23_1_by22, d2G31_2_by11, d2G31_2_by12,...
    d2G31_2_by13, d2G31_2_by22, d2G31_2_by23, d2G31_2_by33, d2G33_1_by11,...
    d2G33_1_by12, d2G33_1_by13, d2G33_1_by22, d2G33_1_by23, d2G33_1_by33,...
    d2G33_2_by11, d2G33_2_by12, d2G33_2_by13, d2G33_2_by22, d2G33_2_by23,...
    d2G33_2_by33, dG11_2_by1, dG11_2_by2, dG11_3_by1, dG11_3_by2,...
    dG11_3_by3, dG12_3_by1, dG12_3_by2, dG12_3_by3, dG22_1_by1,...
    dG22_1_by2, dG22_1_by3, dG22_3_by2, dG22_3_by3, dG23_1_by1,...
    dG23_1_by2, dG31_2_by1, dG31_2_by2, dG31_2_by3, dG33_1_by1,...
    dG33_1_by2, dG33_1_by3, dG33_2_by1, dG33_2_by2, dG33_2_by3];
 
[termVec, gVec] = coeffs(test5, Gset); %length(gVec) =315
save('DataTempV2_Oct2_2019.mat');
% 
% %%
% load('DataTemp_Oct2_2019.mat');
% assumeAlso(Gset,'real');
myCount = 0;
differenceVec = sym('diff',[1, length(gVec)]);
for j=1:length(gVec)
    temp = termVec(j)-conj(termVec(j));
    temp = complex_simple3(temp,[u]);
    differenceVec(j) = temp;
    if temp==0
        disp('good');
        myCount = myCount+1;
    else
        disp('sorry');
    end
end

%%
load('DataTempV2_Oct2_2019.mat');
assumeAlso(Gset,'real');
test6 = d2aV_MuconjMuSub + conj(d2aV_MuconjMuSub);

test6  = subs(test6 , d2G12_3N, d2G12_3Vec);
test6  = subs(test6 , d2G23_1N, d2G23_1Vec);
test6  = subs(test6 , d2G31_2N, d2G31_2Vec);
test6  = subs(test6 , d2G11_2N, d2G11_2Vec);
test6  = subs(test6 , d2G11_3N, d2G11_3Vec);
test6  = subs(test6 , d2G22_1N, d2G22_1Vec);
test6  = subs(test6 , d2G22_3N, d2G22_3Vec);
test6  = subs(test6 , d2G33_1N, d2G33_1Vec);
test6  = subs(test6 , d2G33_2N, d2G33_2Vec);
test6  = subs(test6 , variableSetMu, subSetMu);
test6  = subs(test6 , variableSetConjMu, subSetConjMu);
test6  = subs(test6 , variableSetVnormv, subSetVnormv);
test6  = subs(test6 , Y0, u*conj(u)+1);
% Bianchi Identities
test6  = subs(test6 , variableSetBianchi2, subSetBianchi2);
test6  = subs(test6 , variableSetBianchi1, subSetBianchi1);
test6  = complex_simple3(test6 ,[u]);
test6  = subs(test6 , variableSetBianchi2, subSetBianchi2);
test6  = subs(test6 , variableSetBianchi1, subSetBianchi1);
test6  = complex_simple3(test6 ,[u]);

myCount2 = 0;
[termVec2, gVec2] = coeffs(test6, Gset);
differenceVec2 = sym('diff',[1, length(gVec2)]);
for j=1:length(gVec2)
    temp = termVec2(j);
    temp = complex_simple3(temp,[u]);
    differenceVec2(j) = temp;
    if temp==0
        disp('good');
        myCount2 = myCount2+1;
    else
        disp('sorry');
    end
end
