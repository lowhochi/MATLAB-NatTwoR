% Weyl2_NatTwoR_CR_rmW.m
load('DataWeyl_NatTwoR_CR_rmW.mat');
symSetWeylTwo = symvar(WeylTwo);

wSetTwo = [d2w_MuMu, d2w_MuconjMu, d2w_Muvnormv, d2w_conjMuconjMu,...
    d2w_conjMuvnormv, d2w_uMu, d2w_uconjMu, d2w_uu, d2w_uvnormv,...
    d2w_vnormvvnormv, d3w_uMuMu, d3w_uMuconjMu, d3w_uMuvnormv,...
    d3w_uconjMuconjMu, d3w_uconjMuvnormv, d3w_uuMu, d3w_uuconjMu,...
    d3w_uuu, d3w_uuvnormv, d4w_uuMuMu, d4w_uuMuconjMu, d4w_uuconjMuconjMu,...
    d4w_uuuMu, d4w_uuuconjMu, d4w_uuuu, dw_Mu, dw_conjMu, dw_u, dw_vnormv,...
    w, u];

variableSet5  = [d2aV_conjuconju, d2aV_uconju, d3aV_uuu, d3aV_uuconju];

variableSet6 = [d3aV_conjuconjuMu, d3aV_conjuconjuconjMu, d3aV_conjuconjuvnormv,...
    d3aV_conjuconjuconju, d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv,...
    d3aV_uconjuconju, d4aV_uuconjuMu, d4aV_uuconjuconjMu, d4aV_uuconjuconju,...
    d4aV_uuuconju, d4aV_uuconjuvnormv];

variableSet7 = [d4aV_uconjuMuMu, d4aV_uconjuMuconjMu, d4aV_uconjuMuvnormv,...
    d4aV_uconjuconjMuconjMu, d4aV_uconjuconjMuvnormv, d4aV_uconjuconjuMu,...
    d4aV_uconjuconjuconjMu, d4aV_uconjuconjuconju, d4aV_uconjuconjuvnormv,...
    d4aV_uconjuvnormvvnormv];

MVarWeylTwo =[aV, daV_Mu, daV_conjMu, daV_vnormv, daV_u, daV_conju,...
    d2aV_MuMu, d2aV_MuconjMu, d2aV_Muvnormv, d2aV_conjMuconjMu,...
    d2aV_conjMuvnormv, d2aV_conjuMu, d2aV_conjuconjMu,...
    d2aV_conjuvnormv, d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv,...
    d2aV_vnormvvnormv, d2aV_uu, d3aV_conjuMuMu, d3aV_conjuMuconjMu,...
    d3aV_conjuMuvnormv, d3aV_conjuconjMuconjMu, d3aV_conjuconjMuvnormv,...
    d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv, d3aV_uconjMuvnormv,...
    d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv, ...
    d2theta_MuMu, d2theta_MuconjMu, dtheta_Mu, ...
    dtheta_vnormv, theta];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet5: %ok
d2aV_conjuconjuR = -(2*u/(1+u*conj(u)))*(daV_conju +conj(daV_u))-conj(d2aV_uu);
d2aV_uconjuR = (2/(1+u*conj(u))^2)*(2*i*theta+conj(aV)-2*aV);
d3aV_uuuR = (-6*conj(u)^2/(1+u*conj(u))^2)*daV_u...
    - (6*conj(u)/(1+u*conj(u)))*d2aV_uu;
d3aV_uuconjuR = -4*conj(u)/(1+u*conj(u))^3*(2*i*theta + conj(aV) - 2*aV)...
    + 2/(1+u*conj(u))^2*(-2*daV_u + conj(daV_conju));

% subSet5 = [d2aV_conjuconjuR, d2aV_uconjuR, d3aV_uuuR, d3aV_uuconjuR];
% symvar([d2aV_uconjuR, d2aV_conjuconjuR, d3aV_uuuR, d3aV_uuconjuR]) = 
%   [ aV, d2aV_uu, daV_conju, daV_u, theta, u]

% aMSubTwo = -(1+u*conj(u))^2*daV_u;
% T4SubTwo = ((1+u*conj(u))^2/2)*daV_conju;
% bVSubTwo = u*(1+u*conj(u))*conj(daV_u) + (1+u*conj(u))^2/2*conj(d2aV_uu);
h11SubTwo = i*(1+u*conj(u))^2*aV + (i/2)*(1+u*conj(u))^4*d2aV_uconjuR;
% variableSet1 = [aM, bV, T4, h11, rho]; %last
subSet1Two =[aMSub, bVSub, T4Sub, h11SubTwo, rhoSub]; %last

CVarWeylTwo =[u, theta, aV, daV_u, daV_conju, d2aV_uu, daV_Mu, daV_conjMu,...
    daV_vnormv, dtheta_Mu, dtheta_vnormv];

myNumber2 = length(CVarWeylTwo);
dCVarWeylTwo = sym('dCVar',[5,myNumber2]);
dCVarWeylTwo(:,1) = [0; 0; 0; 1; 0];%u
dCVarWeylTwo(:,2) = [dtheta_Mu; dtheta_conjMu; dtheta_vnormv; 0; 0]; %theta
dCVarWeylTwo(:,3) = [daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVarWeylTwo(:,4) = [d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv; d2aV_uu; d2aV_uconjuR]; %daV_u
dCVarWeylTwo(:,5) = [d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;
    d2aV_uconjuR; d2aV_conjuconjuR]; %daV_conju
dCVarWeylTwo(:,6) = [d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv; d3aV_uuuR;
    d3aV_uuconjuR]; %d2aV_uu 

dCVarWeylTwo(:,7) = [d2aV_MuMu; d2aV_MuconjMu; d2aV_Muvnormv; d2aV_Muu; 
    d2aV_Muconju]; %daV_Mu
dCVarWeylTwo(:,8) = [d2aV_conjMuMu; d2aV_conjMuconjMu; d2aV_conjMuvnormv;
    d2aV_conjMuu; d2aV_conjMuconju]; %daV_conjMu
dCVarWeylTwo(:,9) = [d2aV_vnormvMu; d2aV_vnormvconjMu; d2aV_vnormvvnormv; 
    d2aV_vnormvu; d2aV_vnormvconju]; %daV_vnormv
dCVarWeylTwo(:,10)=[d2theta_MuMu; d2theta_MuconjMu; d2theta_Muvnormv;
    d2theta_Muu; 0]; %dtheta_Mu
dCVarWeylTwo(:,11)=[d2theta_vnormvMu; d2theta_vnormvconjMu; d2theta_vnormvvnormv;
    d2theta_vnormvu; d2theta_vnormvconju]; %dtheta_vnormv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet6 = [d3aV_conjuconjuMu, d3aV_conjuconjuconjMu, d3aV_conjuconjuvnormv,...
%     d3aV_conjuconjuconju, d3aV_uconjuMu, d3aV_uconjuconjMu, d3aV_uconjuvnormv,...
%     d3aV_uconjuconju, d4aV_uuconjuMu, d4aV_uuconjuconjMu, d4aV_uuconjuconju,...
%     d4aV_uuuconju, d4aV_uuconjuvnormv];

d3aV_conjuconjuVec = df_NatTwo_MuSet_CR_rmW(d2aV_conjuconjuR, CVarWeylTwo, dCVarWeylTwo); 
%ok
d3aV_conjuconjuMuR = d3aV_conjuconjuVec(1) ;
d3aV_conjuconjuconjMuR = d3aV_conjuconjuVec(2);
d3aV_conjuconjuvnormvR = d3aV_conjuconjuVec(3);
d3aV_conjuconjuconjuR = d3aV_conjuconjuVec(5);

d3aV_uconjuVec = df_NatTwo_MuSet_CR_rmW(d2aV_uconjuR, CVarWeylTwo, dCVarWeylTwo); 
%ok
d3aV_uconjuMuR = d3aV_uconjuVec(1);
d3aV_uconjuconjMuR = d3aV_uconjuVec(2);
d3aV_uconjuvnormvR = d3aV_uconjuVec(3);
d3aV_uconjuconjuR = d3aV_uconjuVec(5);

d4aV_uuconjuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uuconjuR, CVarWeylTwo, dCVarWeylTwo);
% ok
d4aV_uuconjuMuR = d4aV_uuconjuVec(1);
d4aV_uuconjuconjMuR = d4aV_uuconjuVec(2);
d4aV_uuconjuvnormvR = d4aV_uuconjuVec(3);
d4aV_uuconjuconjuR = d4aV_uuconjuVec(5);

d4aV_uuuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uuuR, CVarWeylTwo, dCVarWeylTwo);
% ok
d4aV_uuuconjuR = d4aV_uuuVec(5);

% diffVariable6 = symvar([d3aV_uconjuMuR, d3aV_uconjuconjMuR, d3aV_uconjuvnormvR,...
%     d3aV_uconjuconjuR]);
% 
% subSet6 = [d3aV_conjuconjuMuR, d3aV_conjuconjuconjMuR, d3aV_conjuconjuvnormvR,...
%      d3aV_conjuconjuconjuR, d3aV_uconjuMuR, d3aV_uconjuconjMuR, d3aV_uconjuvnormvR,...
%      d3aV_uconjuconjuR, d4aV_uuconjuMuR, d4aV_uuconjuconjMuR, d4aV_uuconjuconjuR,...
%      d4aV_uuuconjuR, d4aV_uuconjuvnormvR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet7 = [d4aV_uconjuMuMu, d4aV_uconjuMuconjMu, d4aV_uconjuMuvnormv,...
%     d4aV_uconjuconjMuconjMu, d4aV_uconjuconjMuvnormv, d4aV_uconjuconjuMu,...
%     d4aV_uconjuconjuconjMu, d4aV_uconjuconjuconju, d4aV_uconjuconjuvnormv,...
%     d4aV_uconjuvnormvvnormv];

d4aV_uconjuMuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuMuR, CVarWeylTwo, dCVarWeylTwo);
% ok
d4aV_uconjuMuMuR = d4aV_uconjuMuVec(1);
d4aV_uconjuMuconjMuR = d4aV_uconjuMuVec(2);
d4aV_uconjuMuvnormvR =d4aV_uconjuMuVec(3);

d4aV_uconjuconjMuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuconjMuR, CVarWeylTwo, dCVarWeylTwo);
% ok
d4aV_uconjuconjMuconjMuR = d4aV_uconjuconjMuVec(2);
d4aV_uconjuconjMuvnormvR = d4aV_uconjuconjMuVec(3);

d4aV_uconjuvnormvVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuvnormvR, CVarWeylTwo, dCVarWeylTwo);
% ok
d4aV_uconjuvnormvvnormvR = d4aV_uconjuvnormvVec(3);

d4aV_uconjuconjuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuconjuR, CVarWeylTwo, dCVarWeylTwo);
% ok
d4aV_uconjuconjuMuR = d4aV_uconjuconjuVec(1);
d4aV_uconjuconjuconjMuR = d4aV_uconjuconjuVec(2);
d4aV_uconjuconjuvnormvR = d4aV_uconjuconjuVec(3);
d4aV_uconjuconjuconjuR = d4aV_uconjuconjuVec(5);

% subSet7 = [d4aV_uconjuMuMuR, d4aV_uconjuMuconjMuR, d4aV_uconjuMuvnormvR,...
%      d4aV_uconjuconjMuconjMuR, d4aV_uconjuconjMuvnormvR, d4aV_uconjuconjuMuR,...
%      d4aV_uconjuconjuconjMuR, d4aV_uconjuconjuconjuR, d4aV_uconjuconjuvnormvR,...
%      d4aV_uconjuvnormvvnormvR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MVarWeylTwo = [MVarWeylTwo, wSetTwo];

for j=1:120
    m = WeylTwo(j,1);
    n = WeylTwo(j,2);
    k = WeylTwo(j,3);
    ll =WeylTwo(j,4);
    temp = WeylTwo(j,5);
    for k=1:length(variableSet7)
        myChar = char(variableSet7(k));
        myCharR = [myChar,'R'];
        eval(['temp=subs(temp,', myChar,',',myCharR,');']);
    end
    
    for k=1:length(variableSet6)
        myChar = char(variableSet6(k));
        myCharR = [myChar,'R'];
        eval(['temp=subs(temp,', myChar,',',myCharR,');']);
    end
    
    for k=1:length(variableSet5)
        myChar = char(variableSet5(k));
        myCharR = [myChar,'R'];
        eval(['temp=subs(temp,', myChar,',',myCharR,');']);
    end
    temp = subs(temp, variableSet1, subSet1Two);
    temp = subs(temp, conj(dtheta_vnormv), dtheta_vnormv);
    temp = subs(temp, conj(theta), theta);
    % Y=1+u*conj(u);
%     temp = subs(temp, u*conj(u), Y-1);
%     temp = complex_simple3(temp, MVarWeylTwo);
%     temp = subs(temp, u*conj(u), Y-1);
    temp = complex_simple3(temp, MVarWeylTwo);      
    WeylTwo(j,5)=temp;
end

clearvars temp j m n k ll
clearvars myNumber myNumber2 subSet1 subSet1Two subSet2 subSet3 subSet4
clearvars diffVariable2 diffVariable3 
clearvars variableSet1 variableSet2 variableSet3 variableSet4
clearvars variableSet5 variableSet6 variableSet7

clearvars drhoVec dbVVec dh11Vec daMVec dT4Vec d4aV_uuuVec d4aV_uuconjuVec
clearvars d4aV_uconjuvnormvVec d4aV_uconjuconjuVec d4aV_uconjuconjMuVec
clearvars d4aV_uconjuMuVec d3w_uuRow d3aV_uconjuRow d3T4_uconjMuVec d3T4_uconjuVec
clearvars d2rho_conjuVec d3T4_uuVec d3T4_uMuVec d2rho_uVec d2h11_conjMuVec
clearvars d2h11_conjuVec d2h11_uVec d2h11_vnormvVec d2rho_conjMuVec d2aM_conjMuVec
clearvars d2aM_conjuVec d2T4_uVec d2T4_conjMuVec d2T4_MuVec d2aV_uRow d2aV_conjuRow
clearvars d2aM_uVec d2rho_MuVec d3aV_uconjuVec dthetaRow

save('DataWeyl2_NatTwoR_CR_rmW.mat');