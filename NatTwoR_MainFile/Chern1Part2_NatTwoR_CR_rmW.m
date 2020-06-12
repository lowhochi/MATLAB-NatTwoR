% chern_NatTwoR_CR_rmW_file1Part2.m
load('DataChern1_Part1_NatTwoR_CR_rmW.mat');
% (2) Replacement of higher derivatives of aV in R_in_aV.
syms theta Y0 real % theta = G12_3 + G23_1 + G31_2, Y0 = 1+u*conj(u);

% Replace these variables in MVarChern1:
% V4  [d2aV_conjuconju, d2aV_uconju, d3aV_uuu,...
% V5  d3aV_uconjuconju, d3aV_uuconju, d3aV_uconjuMu,...
%     d3aV_uconjuconjMu, d3aV_uconjuvnormv,...
% V6  d4aV_uconjuMuMu, d4aV_uconjuMuconjMu, d4aV_uconjuMuvnormv,...
%     d4aV_uconjuconjMuconjMu, d4aV_uconjuconjMuvnormv, d4aV_uconjuvnormvvnormv,...
%     d4aV_uuconjuMu, d4aV_uuconjuconjMu, d4aV_uuconjuvnormv,...
%     d4aV_uconjuconjuMu, d4aV_uconjuconjuconjMu, d4aV_uconjuconjuvnormv,...
% V7  d4aV_uuuconju, d4aV_uuconjuconju, d4aV_uconjuconjuconju];

d2aV_conjuconjuR = -(2*u/(1+u*conj(u)))*(daV_conju +conj(daV_u))-conj(d2aV_uu);
d2aV_uconjuR = (2/(1+u*conj(u))^2)*(2*i*theta+conj(aV)-2*aV);
d3aV_uuuR = (-6*conj(u)^2/(1+u*conj(u))^2)*daV_u...
    - (6*conj(u)/(1+u*conj(u)))*d2aV_uu;
d3aV_uuconjuR = -4*conj(u)/(1+u*conj(u))^3*(2*i*theta + conj(aV) - 2*aV)...
    + 2/(1+u*conj(u))^2*(-2*daV_u + conj(daV_conju));

variable_NatTwoR_Chern1_Part2_CR_rmW

variableSet4 = [d2aV_conjuconju, d2aV_uconju, d3aV_uuu]; 
subSet4 = [d2aV_conjuconjuR, d2aV_uconjuR, d3aV_uuuR];

% %
d3aV_uconjuVec = df_NatTwo_MuSet_CR_rmW(d2aV_uconjuR, CVarChern1S, dCVarChern1S);
% d3aV_uuconjuR = d3aV_uconjuVec(4);
d3aV_uconjuconjuR = d3aV_uconjuVec(5);
d3aV_uconjuMuR = d3aV_uconjuVec(1);
d3aV_uconjuconjMuR = d3aV_uconjuVec(2);
d3aV_uconjuvnormvR = d3aV_uconjuVec(3);

variableSet5 = [d3aV_uuconju, d3aV_uconjuconju, d3aV_uconjuMu,...
   d3aV_uconjuconjMu, d3aV_uconjuvnormv];
subSet5 = [d3aV_uuconjuR, d3aV_uconjuconjuR, d3aV_uconjuMuR,...
   d3aV_uconjuconjMuR, d3aV_uconjuvnormvR];

% Differentiate further on [d3aV_uconjuMuR, d3aV_uconjuconjMuR, 
%   d3aV_uconjuvnormvR, d3aV_uuconjuR, d3aV_uconjuconjuR];

% d3aV_uconjuMuR = d3aV_uconjuVec(1);
d4aV_uconjuMuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuMuR, CVarChern1S, dCVarChern1S);
d4aV_uconjuMuMuR = d4aV_uconjuMuVec(1);
d4aV_uconjuMuconjMuR = d4aV_uconjuMuVec(2);
d4aV_uconjuMuvnormvR = d4aV_uconjuMuVec(3);
% d3aV_uconjuconjMuR = d3aV_uconjuVec(2);
d4aV_uconjuconjMuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuconjMuR, CVarChern1S, dCVarChern1S);
d4aV_uconjuconjMuconjMuR = d4aV_uconjuconjMuVec(2);
d4aV_uconjuconjMuvnormvR = d4aV_uconjuconjMuVec(3);
% d3aV_uconjuvnormv = d3aV_uconjuVec(3);
d4aV_uconjuvnormvVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuvnormvR, CVarChern1S, dCVarChern1S);
d4aV_uconjuvnormvvnormvR = d4aV_uconjuvnormvVec(3);
% d3aV_uuconjuR = d3aV_uconjuVec(4);
d4aV_uuconjuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uuconjuR, CVarChern1S, dCVarChern1S);
d4aV_uuconjuMuR = d4aV_uuconjuVec(1);
d4aV_uuconjuconjMuR = d4aV_uuconjuVec(2);
d4aV_uuconjuvnormvR = d4aV_uuconjuVec(3);
d4aV_uuuconjuR = d4aV_uuconjuVec(4);
d4aV_uuconjuconjuR = d4aV_uuconjuVec(5);

% d3aV_uconjuconjuR = d3aV_uconjuVec(5);
d4aV_uconjuconjuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuconjuR, CVarChern1S, dCVarChern1S);
d4aV_uconjuconjuMuR = d4aV_uconjuconjuVec(1);
d4aV_uconjuconjuconjMuR = d4aV_uconjuconjuVec(2);
d4aV_uconjuconjuvnormvR = d4aV_uconjuconjuVec(3);
d4aV_uconjuconjuconjuR = d4aV_uconjuconjuVec(5);

variableSet6 = [d4aV_uconjuMuMu, d4aV_uconjuMuconjMu, d4aV_uconjuMuvnormv,...
    d4aV_uconjuconjMuconjMu, d4aV_uconjuconjMuvnormv, d4aV_uconjuvnormvvnormv,...
    d4aV_uuconjuMu, d4aV_uuconjuconjMu, d4aV_uuconjuvnormv, d4aV_uconjuconjuMu,...
    d4aV_uconjuconjuconjMu, d4aV_uconjuconjuvnormv];

subSet6 = [d4aV_uconjuMuMuR, d4aV_uconjuMuconjMuR, d4aV_uconjuMuvnormvR,...
    d4aV_uconjuconjMuconjMuR, d4aV_uconjuconjMuvnormvR, d4aV_uconjuvnormvvnormvR,...
    d4aV_uuconjuMuR, d4aV_uuconjuconjMuR, d4aV_uuconjuvnormvR, d4aV_uconjuconjuMuR,...
    d4aV_uconjuconjuconjMuR, d4aV_uconjuconjuvnormvR];

% %
% d4aVTemp0 = -6*conj(u)^2/(1+u*conj(u))^2;
% d4aVTemp1 = 6*conj(u)/(1+u*conj(u));
% d4aVTemp2 = 2/(1+u*conj(u))^2;
% d4aVTemp3 = complexdiff3(d4aVTemp2,u,0);
% d4aVTemp4 = complexdiff3(d4aVTemp3,u,1);
% d4aVTemp5 = complexdiff3(d4aVTemp3,u,0);

% d3aV_uuconjuR = d4aVTemp3*(2*i*theta +conj(aV) -2*aV)...
%    +d4aVTemp2*(conj(daV_conju)-2*daV_u);
% d4aV_uuuconjuR = complexdiff3(d4aVTemp0,u,1)*daV_u...
%     - complexdiff3(d4aVTemp1,u,1)*d2aV_uu...
%     + d4aVTemp0*d2aV_uconjuR - d4aVTemp1*d3aV_uuconjuR;

% d4aV_uuconjuconjuR = d4aVTemp4*(2*i*theta+conj(aV)-2*aV)...
%     + d4aVTemp3*(conj(daV_u)-2*daV_conju)...
%     + conj(d4aVTemp3)*(conj(daV_conju)-2*daV_u)...
%     + d4aVTemp2*(conj(d2aV_uconjuR)-2*d2aV_uconjuR);

% d4aV_uconjuconjuconjuR = conj(d4aVTemp5)*(2*i*theta+conj(aV)-2*aV)...
%     + 2*conj(d4aVTemp3)*(conj(daV_u)-2*daV_conju)...
%     + d4aVTemp2*(conj(d2aV_uu)-2*d2aV_conjuconjuR);

variableSet7 = [d4aV_uuuconju, d4aV_uuconjuconju, d4aV_uconjuconjuconju];
subSet7 = [d4aV_uuuconjuR, d4aV_uuconjuconjuR, d4aV_uconjuconjuconjuR];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Substitution

for j=1:16
    m = indexChern(j,1);
    n = indexChern(j,2);
    k = indexChern(j,3);
    ll = indexChern(j,4);
    temp = Chern_in_aV(m,n,k,ll);
    temp = subs(temp, variableSet7, subSet7);
    temp = subs(temp, variableSet6, subSet6);
    temp = subs(temp, variableSet5, subSet5);
    temp = subs(temp, variableSet4, subSet4);
    % Y0=1+u*conj(u);
    temp = subs(temp, u*conj(u), Y0-1);
    temp = subs(temp, variableSet1, subSet1Two);
    temp = complex_simple3(temp, MVarChern1S);
    Chern_in_aV(m,n,k,ll) = temp;
end

clearvars j m n k ll temp 
clearvars d4aVTemp0 d4aVTemp1 d4aVTemp2 d4aVTemp3 d4aVTemp4 d4aVTemp5
clearvars variableSet1 variableSet2 variableSet3 variableSet4
clearvars variableSet5 variableSet6 variableSet7
clearvars subSet1 subSet2 subSet3 subSet4 subSet5 subSet6 subSet7   

clearvars subSet1Two daVRow d2aV_uRow d3aV_uconjuRow
clearvars dthetaRow

list_of_variables_Chern1 = who;

save('DataChern1_Part2_NatTwoR_CR_rmW.mat');