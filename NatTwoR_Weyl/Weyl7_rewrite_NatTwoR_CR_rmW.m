load('DataWeyl2_NatTwoR_CR_rmW.mat');
Weylf = sym('Weylf',[120,5]);
for j=1:120
    for k=1:5
        Weylf(j,k) = WeylTwo(j,k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define d4aV_uuXX
syms d4aV_uuMuMu d4aV_uuMuconjMu d4aV_uuconjMuconjMu
syms d4aV_uuMuvnormv d4aV_uuconjMuvnormv d4aV_uuvnormvvnormv
d4aV_uuuu = 12*conj(u)^3/(1+u*conj(u))^3*daV_u...
    -6*conj(u)/(1+u*conj(u))*d3aV_uuuR;
% d3aV_uuuR=-(6*d2aV_uu*conj(u))/(u*conj(u)+1)-(6*daV_u*conj(u)^2)/(u*conj(u) + 1)^2;
d4aV_uuuMu= -(6*d3aV_uuMu*conj(u))/(u*conj(u)+1)...
    -(6*d2aV_uMu*conj(u)^2)/(u*conj(u)+1)^2;
d4aV_uuuconjMu = -(6*d3aV_uuconjMu*conj(u))/(u*conj(u)+1)...
    -(6*d2aV_uconjMu*conj(u)^2)/(u*conj(u)+1)^2;
d4aV_uuuvnormv = -(6*d3aV_uuvnormv*conj(u))/(u*conj(u)+1)...
    -(6*d2aV_uvnormv*conj(u)^2)/(u*conj(u)+1)^2;

d3aV_uuRow = [d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv];
lie_Mu_conjMuSub = [aMSub; -conj(aMSub); 2*i*h11SubTwo];
lie_Mu_vnormvSub = [aV; bVSub; 2*T4Sub]; 
lie_conjMu_vnormvSub = [conj(bVSub); conj(aV); 2*conj(T4Sub)];
d4aV_uuconjMuMu = d4aV_uuMuconjMu+ d3aV_uuRow*lie_Mu_conjMuSub;
d4aV_uuvnormvMu = d4aV_uuMuvnormv+ d3aV_uuRow*lie_Mu_vnormvSub;
d4aV_uuvnormvconjMu = d4aV_uuconjMuvnormv +d3aV_uuRow*lie_conjMu_vnormvSub;
d4aV_uuMuu = d4aV_uuuMu -d3aV_uuRow*lie_Mu_ddu;
d4aV_uuconjMuconju = d4aV_uuconjuconjMuR - d3aV_uuRow*lie_conjMu_ddconju;
d4aV_uuvnormvu = d4aV_uuuvnormv -d3aV_uuRow*lie_vnormv_ddu;
d4aV_uuvnormvconju = d4aV_uuconjuvnormvR -d3aV_uuRow*lie_vnormv_ddconju;

variableSet1 = [aM, h11, T4, bV];
subSet1 = [aMSub, h11SubTwo, T4Sub, bVSub];
d3aV_uconjMuMuR = subs(d3aV_uconjMuMu, variableSet1, subSet1);
d3aV_uvnormvMuR = subs(d3aV_uvnormvMu, variableSet1, subSet1);
d3aV_uvnormvconjMuR = subs(d3aV_uvnormvconjMu, variableSet1, subSet1);
d3aV_uvnormvconjuR = subs(d3aV_uvnormvconju, d3aV_uconjuvnormv,...
    d3aV_uconjuvnormvR);
d3aV_uconjMuconjuR = subs(d3aV_uconjMuconju, d3aV_uconjuconjMu, ...
    d3aV_uconjuconjMuR);

d3aV_conjuMuuR = subs(d3aV_conjuMuu, d3aV_uconjuMu, d3aV_uconjuMuR);
d3aV_conjuconjMuMuR = subs(d3aV_conjuconjMuMu, variableSet1, subSet1);
d3aV_conjuconjMuconjuR = subs(d3aV_conjuconjMuconju, d3aV_conjuconjuconjMu,...
    d3aV_conjuconjuconjMuR);
d2aV_conjMuMuR = subs(d2aV_conjMuMu, variableSet1, subSet1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVarW7 = [u, aV, daV_u, daV_conju, d2aV_uu, theta,...
    d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv,...
    d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv,...
    d2aV_conjuMu, d2aV_conjuconjMu,...
    daV_Mu, daV_conjMu, dtheta_Mu];

derivativeDict.u =[0; 0; 0; 1; 0];

derivativeDict.aV =[daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju];

derivativeDict.daV_u = [d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv;
    d2aV_uu; d2aV_uconjuR];

derivativeDict.daV_conju = [d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;
    d2aV_uconjuR; d2aV_conjuconjuR];

derivativeDict.d2aV_uu = [d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv;
    d3aV_uuuR; d3aV_uuconjuR];

derivativeDict.theta = [dtheta_Mu; dtheta_conjMu; dtheta_vnormv; 0; 0];

derivativeDict.d2aV_uMu = [d3aV_uMuMu; d3aV_uMuconjMu; d3aV_uMuvnormv;
    d3aV_uMuu; d3aV_uconjuMuR];

derivativeDict.d2aV_uconjMu =[d3aV_uconjMuMuR; d3aV_uconjMuconjMu;
    d3aV_uconjMuvnormv; d3aV_uuconjMu; d3aV_uconjMuconjuR]; %h11 aM

derivativeDict.d2aV_uvnormv = [d3aV_uvnormvMuR; d3aV_uvnormvconjMuR;
    d3aV_uvnormvvnormv; d3aV_uvnormvu; d3aV_uvnormvconjuR]; %T4 bV

derivativeDict.d3aV_uuMu = [d4aV_uuMuMu; d4aV_uuMuconjMu; d4aV_uuMuvnormv;
    d4aV_uuMuu; d4aV_uuconjuMuR];

derivativeDict.d3aV_uuconjMu = [d4aV_uuconjMuMu; d4aV_uuconjMuconjMu; 
    d4aV_uuconjMuvnormv; d4aV_uuuconjMu; d4aV_uuconjMuconju];

derivativeDict.d3aV_uuvnormv = [d4aV_uuvnormvMu; d4aV_uuvnormvconjMu; 
    d4aV_uuvnormvvnormv; d4aV_uuvnormvu; d4aV_uuvnormvconju];

derivativeDict.d2aV_conjuMu = [d3aV_conjuMuMu; d3aV_conjuMuconjMu;
    d3aV_conjuMuvnormv; d3aV_conjuMuuR; d3aV_conjuconjuMuR]; 

derivativeDict.d2aV_conjuconjMu = [d3aV_conjuconjMuMuR; d3aV_conjuconjMuconjMu;
    d3aV_conjuconjMuvnormv; d3aV_uconjuconjMuR; d3aV_conjuconjMuconjuR];

derivativeDict.daV_Mu = [d2aV_MuMu; d2aV_MuconjMu; d2aV_Muvnormv;
    d2aV_Muu; d2aV_conjuMu];

derivativeDict.daV_conjMu = [d2aV_conjMuMuR; d2aV_conjMuconjMu; 
    d2aV_conjMuvnormv; d2aV_uconjMu; d2aV_conjMuconju];

derivativeDict.dtheta_Mu = [d2theta_MuMu; d2theta_MuconjMu; d2theta_Muvnormv;
    d2theta_Muu; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace w by f
wSet1 = [w, dw_Mu, dw_conjMu, dw_vnormv, dw_u];
f =u*(1+u*conj(u))^3*conj(daV_u)+(1+u*conj(u))^4/2*conj(d2aV_uu);

dfVec = df_main4_CR_funWv2(f,CVarW7,derivativeDict,gamma);
df_Mu = dfVec(1);
df_conjMu = dfVec(2);
df_vnormv = dfVec(3);
df_u = dfVec(4);

% % break point % %
wSet2 = [d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_uu,...
    d2w_MuMu, d2w_MuconjMu, d2w_Muvnormv, d2w_conjMuconjMu,...
    d2w_conjMuvnormv, d2w_vnormvvnormv];
    
d2f_uVec = df_main4_CR_funWv2(df_u,CVarW7,derivativeDict,gamma);
d2f_uMu = d2f_uVec(1);
d2f_uconjMu = d2f_uVec(2);
d2f_uvnormv = d2f_uVec(3);
d2f_uu = d2f_uVec(4);
d2f_MuVec = df_main4_CR_funWv2(df_Mu,CVarW7,derivativeDict,gamma);
d2f_conjMuVec = df_main4_CR_funWv2(df_conjMu,CVarW7,derivativeDict,gamma);
d2f_vnormvVec = df_main4_CR_funWv2(df_vnormv,CVarW7,derivativeDict,gamma);
d2f_MuMu = d2f_MuVec(1);
d2f_MuconjMu = d2f_MuVec(2);
d2f_Muvnormv = d2f_MuVec(3);
d2f_conjMuconjMu = d2f_conjMuVec(2);
d2f_conjMuvnormv = d2f_conjMuVec(3);
d2f_vnormvvnormv = d2f_vnormvVec(3);

% % break point % %
wSet3 = [d3w_uMuMu, d3w_uMuconjMu, d3w_uMuvnormv, d3w_uconjMuconjMu,...
    d3w_uconjMuvnormv, d3w_uuMu, d3w_uuconjMu, d3w_uuvnormv, d3w_uuu];

d3f_uMuVec = df_main4_CR_funWv2(d2f_uMu,CVarW7,derivativeDict,gamma);
d3f_uconjMuVec = df_main4_CR_funWv2(d2f_uconjMu,CVarW7,derivativeDict,gamma);
d3f_uMuMu = d3f_uMuVec(1);
d3f_uMuconjMu = d3f_uMuVec(2);
d3f_uMuvnormv = d3f_uMuVec(3);
d3f_uconjMuconjMu = d3f_uconjMuVec(2);
d3f_uconjMuvnormv = d3f_uconjMuVec(3);

d3f_uuVec = df_main4_CR_funWv2(d2f_uu,CVarW7,derivativeDict,gamma);
d3f_uuMu = d3f_uuVec(1);
d3f_uuconjMu = d3f_uuVec(2);
d3f_uuvnormv = d3f_uuVec(3);
d3f_uuu = d3f_uuVec(4);

% % break point % %
wSet4 = [d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuconjMuconjMu, d4w_uuuMu,...
    d4w_uuuconjMu, d4w_uuuu]; 
d4f_uuMuVec = df_main4_CR_funWv2(d3f_uuMu,CVarW7,derivativeDict,gamma);
d4f_uuconjMuVec = df_main4_CR_funWv2(d3f_uuconjMu,CVarW7,derivativeDict,gamma);
d4f_uuuVec = df_main4_CR_funWv2(d3f_uuu,CVarW7,derivativeDict,gamma);

d4f_uuMuMu = d4f_uuMuVec(1);
d4f_uuMuconjMu = d4f_uuMuVec(2);
d4f_uuconjMuconjMu = d4f_uuconjMuVec(2);
d4f_uuuMu = d4f_uuuVec(1);
d4f_uuuconjMu = d4f_uuuVec(2);
d4f_uuuu = d4f_uuuVec(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ric12Sub = d2aV_conjuconjMu...
%     -(1+u*conj(u))^2*(daV_u*daV_conju)...
%     +(1+u*conj(u))^2*conj(daV_u)*conj(daV_conju)...
%     -2*daV_vnormv +2*conj(daV_vnormv) +4*i*dtheta_vnormv...
%     +4*i*theta*(aV+conj(aV)) +2*conj(aV)^2 - 2*aV^2;

for j=1:120
    temp = Weylf(j,5);
    for k=1:length(wSet4)
        wChar = char(wSet4(k));
        fChar = strrep(wChar,'w','f');
        eval(['temp=subs(temp,', wChar, ',', fChar, ');']);
    end
    
    for k=1:length(wSet3)
        wChar = char(wSet3(k));
        fChar = strrep(wChar,'w','f');
        eval(['temp=subs(temp,', wChar, ',', fChar, ');']);
    end
    
    for k=1:length(wSet2)
        wChar = char(wSet2(k));
        fChar = strrep(wChar,'w','f');
        eval(['temp=subs(temp,', wChar, ',', fChar, ');']);
    end
    
    for k=1:length(wSet1)
        wChar = char(wSet1(k));
        fChar = strrep(wChar,'w','f');
        eval(['temp=subs(temp,', wChar, ',', fChar, ');']);
    end 
%     temp = subs(temp, conj(d2aV_conjuconjMu),Ric12Sub);
    Weylf(j,5) = temp;
end

clearvars wChat wSet1 wSet2 wSet3 wSet4 wSet fChar variableSet1 subSet1
clearvars temp j k 
save('DataWeyl7new_NatTwoR_CR_rmW.mat');