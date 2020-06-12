% Weyl3Part2_NatTwoR_CR_rmW.m
load('DataWeyl2_NatTwoR_CR_rmW.mat');
indexChern = [1,1,1,1; 1,1,1,2; 1,1,2,1; 1,1,2,2;
    1,2,1,1; 1,2,1,2; 1,2,2,1; 1,2,2,2;
    2,1,1,1; 2,1,1,2; 2,1,2,1; 2,1,2,2;
    2,2,1,1; 2,2,1,2; 2,2,2,1; 2,2,2,2];

%load('Data_Chern_in_aV_Oct10.mat');
%Wijkl = Weyl(ui,uj,uk,ul), represents the coefficient 'Cxxxx'.
W1212 = WeylTwo(1,5); %C1111, Weyl(u1,u2,u1,u2)
W1214 = WeylTwo(3,5); %C1112, Weyl(u1,u2,u1,u4)
W1223 = WeylTwo(6,5); %-C1121, -Weyl(u1,u2,u3,u2)
W1234 = WeylTwo(10,5); %C1122, Weyl(u1,u2,u3,u4)
W1412 = W1214; %C1211, Weyl(u1,u4,u1,u2)=Weyl(u1,u2,u1,u4)
W1414 = WeylTwo(30,5); %C1212, Weyl(u1,u4,u1,u4)
W1423 = WeylTwo(33,5); %-C1221, -Weyl(u1,u4,u3,u2)
W1434 = WeylTwo(37,5); %C1222, Weyl(u1,u4,u3,u4)
W2312 = W1223; %-C2111, -Weyl(u3,u2,u1,u2)
W2314 = W1423; %-C2112, -Weyl(u3,u2,u1,u4)
W2323 = WeylTwo(66,5); %C2121, Weyl(u3,u2,u3,u2)
W2334 = WeylTwo(70,5); %-C2122, -Weyl(u3,u2,u3,u4)
W3412 = W1234; %C2211, Weyl(u3,u4,u1,u2)
W3414 = W1434; %C2212, Weyl(u3,u4,u1,u4)
W3423 =  W2334; %-C2221, -Weyl(u3,u4,u3,u2)
W3434 = WeylTwo(100,5); %C2222, Weyl(u3,u4,u3,u4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test(1) C2221, C2212, C2122, C1222, C2222 = 0;
% WtestVec1 = [-W3423, W3414, -W2334, W1434, W3434];
% disp(WtestVec1);

% Test(2) C1122, C1212, C1221, C2112, C2121, C2211 = -rho/6;
% WtestVec2 = [W1234, W1414, -W1423, -W2314, W2323, W3412];
% diffVec2 = sym('diffVec2',[1,6]);
% for j=1:6
%     temp = WtestVec2(j)+rhoSub/6;
%     diffVec2(j) = complex_simple3(temp,MVarWeylTwo);
% end

% Test(3) C1112 = C1211, C1121 = C2111 and C1112 = conj(C1121);
% WtestVec3 = [W1214-W1412, -W1223+W2312, W1214+conj(W1223)];
% diffVec3 = sym('diffVec3',[1,3]);
% for j=1:3
%     temp = WtestVec3(j);
%     temp = subs(temp, conj(theta),theta);
%     temp = subs(temp, conj(dtheta_vnormv), dtheta_vnormv);
%     diffVec3(j) = complex_simple3(temp,MVarWeylTwo);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test(4) C1112 = W1214 and C1111 = W1212
load('Data_Chern_in_aV_Oct10.mat');
C1111 = Chern_in_aV(1,1,1,1);
C1112 = Chern_in_aV(1,1,1,2);

symSetTest4 = symvar([C1111, C1112, W1214, W1212]);
WtestVec4 = [W1212-C1111, W1214-C1112];

% 'Ric12'-substitution
KtermSub = d2aV_conjuconjMu...
    -(1+u*conj(u))^2*(daV_u*daV_conju)...
    +(1+u*conj(u))^2*conj(daV_u)*conj(daV_conju)...
    -2*daV_vnormv +2*conj(daV_vnormv) +4*i*dtheta_vnormv...
    +4*i*theta*(aV+conj(aV)) +2*conj(aV)^2 - 2*aV^2;
WtestVec4(1) = subs(WtestVec4(1),conj(d2aV_conjuconjMu),KtermSub);
WtestVec4(1) = complex_simple3(WtestVec4(1),symSetWeylTwo);

variable_NatTwoR_Chern2_CR_rmW % Use the Chern2 file

variableSet0 = [aM, bV, h11, T4];
subSet0 = [aMSub, bVSub, T4Sub, h11SubTwo];
% substitute variable0 BEFORE [theta, aV, daV_u, daV_conju, d2aV_uu].

wSetTest4 = [u, w, dw_conjMu, dw_u, dw_vnormv, d2w_conjMuconjMu,...
    d2w_uconjMu, d2w_uu];

variableSetG1 = [theta, aV, daV_u, daV_conju, d2aV_uu];
variableSetG2 = [daV_Mu, daV_conjMu, daV_vnormv, dtheta_Mu, dtheta_vnormv,...
    d2aV_uMu, d2aV_uconjMu, d2aV_conjuMu, d2aV_conjuconjMu];
variableSetG3 = [d2aV_MuconjMu, d2theta_MuconjMu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% follow from the m-file Chern2_NatTwoR_CR_rmW.m
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

subSetG1 = [thetaSub, aVSub, daV_uSub, daV_conjuSub, d2aV_uuSub];

subSetG2 = [daV_MuSub, daV_conjMuSub, daV_vnormvSub, dtheta_MuSub, dtheta_vnormvSub,...
    d2aV_uMuSub, d2aV_uconjMuSub, d2aV_conjuMuSub, d2aV_conjuconjMuSub];

subSetG3 = [d2aV_MuconjMuSub, d2theta_MuconjMuSub];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:2
    temp = WtestVec4(j);
    temp = subs(temp, variableSetG3, subSetG3);
    temp = subs(temp, variableSetG2, subSetG2);
    temp = subs(temp, variableSet0, subSet0);
    temp = subs(temp, variableSetG1, subSetG1);
    % Second substitution
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
    % Bianchi Identity
    temp = subs(temp, variableSetBianchi2, subSetBianchi2);
    temp = subs(temp, variableSetBianchi1, subSetBianchi1);
    temp = complex_simple3(temp, wSetTest4);
    WtestVec4(j)=temp;
end

[termVec, gVec]=coeffs(WtestVec4(1),Gijk);
for j=1:length(gVec)
    temp = termVec(j);
    temp = complex_simple3(temp,u);
    termVec(j)=temp;
end

clearvars wSetTest4 variableSet0 variableSetG1 variableSetG2 variableSetG3
clearvars subSet0 subSetG1 subSetG2 subSetG3
clearvars symSetTest4 termVec temp j gVec daVVec
clearvars dG33_1Row dG33_2Row dG12_3Row dG23_1Row dG31_2Row 
clearvars dG11_2Row dG11_3Row dG22_1Row dG22_3Row

clearvars dG12_3ERow dG23_1ERow dG31_2ERow dG11_2ERow dG11_3ERow
clearvars dG22_1ERow dG22_3ERow dG33_1ERow dG33_2ERow

clearvars d3aV_conjuconjuVec d2h11_MuVec
clearvars KtermSub 

% save('DataTemp_Oct11.mat');
% save('DataWeyl3_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Focus on W1212 - C1111.
% load('DataWeylTwo_NatTwoR_CR_rmW.mat');
% load('Data_Chern_in_aV_Oct10.mat');
% C1111 = Chern_in_aV(1,1,1,1);
% W1212 = WeylTwo(1,5);
% W1212 = subs(W1212, Y, 1+u*conj(u));
% 
% H0= (2*theta + i*aV - i*conj(aV));
% Yk = (1+u*conj(u));
% term1 = H0*Yk^6*(daV_u*daV_conju + daV_u*conj(daV_u));
% 
% term2 = Yk^4*((-16/3)*theta^3 + (-40/3*i*aV + 16/3*i*conj(aV))*theta^2 ...
%     +(28/3*aV^2 - 32/3*aV*conj(aV) + 4/3*conj(aV)^2)*theta ...
%     + 2*i*aV^3 - 4*i*aV^2*conj(aV) + 2*i*aV*conj(aV)^2);
% 
% term3 = Yk^4*H0^2*(i/6*d2w_uu - i/6*conj(d2w_uu))...
%     + Yk^4*H0*(2*daV_u*dw_u + w*d2aV_uu + 2*conj(daV_u)*conj(dw_u)...
%     +conj(w)*conj(d2aV_uu));
% 
% term4 = Yk^4*(2*i*w*(-2*daV_u^2 + daV_u*conj(daV_conju))...
%     +2*i*conj(w)*(-daV_conju*conj(daV_u)+ 2*conj(daV_u)^2) );
% 
% term5 = Yk^4*H0*(-4*i*dtheta_vnormv)...
%     +Yk^4*(daV_u*(-4*dtheta_Mu -2*i*daV_Mu + i*conj(daV_conjMu))...
%         +conj(daV_u)*(-2*conj(dtheta_Mu) - i*daV_conjMu +2*i*conj(daV_Mu)) )...
%     +Yk^4*H0*(-d2aV_conjuconjMu + 2*daV_vnormv);
% 
% term6 = -i*Yk^3*H0^2*(conj(u)*dw_u - u*conj(dw_u))...
%     + Yk^3*H0*(-6*conj(u)*w*daV_u - 6*u*conj(w)*conj(daV_u));
% 
% term7 = Yk^2*H0*(-d2w_uconjMu - conj(d2w_uconjMu))...
%     +Yk^2*(dw_conjMu*(4*i*daV_u - i*conj(daV_conju))...
%         + conj(dw_conjMu)*(i*daV_conju - 4*i*conj(daV_u)))...
%     +Yk^2*dw_u*(-2*conj(dtheta_Mu) + i*conj(daV_Mu) - i*daV_conjMu)...
%     +Yk^2*conj(dw_u)*(-2*dtheta_Mu + i*conj(daV_conjMu) - i*daV_Mu)...
%     +Yk^2*dw_u*conj(w)*(-i*daV_conju + 2*i*conj(daV_u))...
%     +Yk^2*conj(dw_u)*w*(i*conj(daV_conju) - 2*i*daV_u);
% 
% term8 = Yk^2*H0*(dw_u*conj(dw_u))...
%     +Yk^2*w*(2*i*conj(u)^2*H0^2+2*i*d2aV_uconjMu-i*conj(d2aV_conjuMu))...
%     +Yk^2*conj(w)*(-2*i*u^2*H0^2+i*d2aV_conjuMu-2*i*conj(d2aV_uconjMu));
% 
% term9 = Yk^2*(2*d2theta_MuconjMu-i*conj(d2aV_MuconjMu)+2*i*d2aV_MuconjMu)...
%     + 4*Yk*H0*(conj(u)*dw_conjMu + u*conj(dw_conjMu) - u*conj(w)*dw_u...
%         - conj(u)*w*conj(dw_u) + 4*w*conj(w))...
%     + Yk*w*conj(w)*(8*i*(u*daV_u - conj(u)*conj(daV_u))...
%         +4*i*(conj(u)*daV_conju - u*conj(daV_conju)));
% 
% term10 = Yk*(4*u*conj(w)*(2*dtheta_Mu + i*daV_Mu - i*conj(daV_conjMu))...
%     +4*conj(u)*w*(2*conj(dtheta_Mu) + i*daV_conjMu - i*conj(daV_Mu)))...
%     +w*conj(w)*(-40*theta-24*i*aV+24*i*conj(aV));
% 
% term11 = i*(-d2w_conjMuconjMu + conj(d2w_conjMuconjMu))...
%     +2*i*(-conj(w)*dw_vnormv + w*conj(dw_vnormv))...
%     +i*(-dw_u*conj(dw_conjMu) + conj(dw_u)*dw_conjMu)...
%     -(4*i*u/Yk)*conj(w)*dw_conjMu + (4*i*conj(u)/Yk)*w*conj(dw_conjMu);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The correction term: term12
% term12 = d2aV_MuconjMu*(u*conj(u)*1i + 1i)^2*1i;
% 
% test1111 = term1 + term2 + term3 + term4 + term5 + term6...
%     + term7 + term8 + term9 + term10 + term11 + term12;
% 
% Kterm = d2aV_conjuconjMu -2*daV_vnormv +2*i*dtheta_vnormv...
%     +4*i*theta*conj(aV)+2*conj(aV)^2 -(1+u*conj(u))^2*daV_u*daV_conju;
% 
% diff1111 = (1+u*conj(u))^4/2*(2*i*theta+conj(aV)-aV)*(-i*Kterm+i*conj(Kterm));
% 
% test05 = W1212 - test1111 -diff1111;
% test05 = subs(test05,[conj(theta),conj(dtheta_vnormv)],[theta,dtheta_vnormv]);
% test05 = complex_simple3(test05, symSetWeylTwo);

