% RicCurv2_NatTwoR_CR_rmW.m
load('DataRicCurv1_Oct7_NatTwoR_CR_rmW.mat');
clearvars MVarRicCurv1 MVarRicCurv2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symSetRicCurv2 = symvar(RicCurv_in_aV);
syms Y % Y= 1+u*conj(u);

% [u, w, theta,...
%   dtheta_Mu, dtheta_conjMu, dtheta_vnormv,...
%   dw_Mu, dw_conjMu, dw_u, dw_vnormv, d2w_uMu, d2w_uconjMu, d2w_uu,...
%   d2w_uvnormv, d3w_uuMu, d3w_uuconjMu, d3w_uuu, d3w_uuvnormv,...
%   d4w_uuuMu, d4w_uuuu];

% [aV, daV_Mu, daV_conjMu, daV_conju, daV_u, daV_vnormv,...
%   d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuconju,...
%   d2aV_conjuvnormv, d2aV_uMu, d2aV_uconjMu, d2aV_uconju, d2aV_uu,...
%   d2aV_uvnormv,...
%   d3aV_conjuconjuconjMu, d3aV_conjuconjuconju, d3aV_uconjuMu,...
%   d3aV_uconjuconjMu, d3aV_uconjuconju, d3aV_uconjuvnormv,... 
%   d3aV_uuMu, d3aV_uuconju, d3aV_uuu,...
%   d4aV_uconjuconjuMu, d4aV_uconjuconjuconjMu, d4aV_uconjuconjuconju,...
%   d4aV_uconjuconjuvnormv, d4aV_uuconjuMu, d4aV_uuconjuconjMu,...
%   d4aV_uuconjuconju, d4aV_uuconjuvnormv, d4aV_uuuconju];

d2aV_conjuconjuR = -(2*u/(1+u*conj(u)))*(daV_conju +conj(daV_u))-conj(d2aV_uu);
d2aV_uconjuR = (2/(1+u*conj(u))^2)*(2*i*theta+conj(aV)-2*aV);
d3aV_uuuR = (-6*conj(u)^2/(1+u*conj(u))^2)*daV_u...
    - (6*conj(u)/(1+u*conj(u)))*d2aV_uu;
d3aV_uuconjuR = -4*conj(u)/(1+u*conj(u))^3*(2*i*theta + conj(aV) - 2*aV)...
    + 2/(1+u*conj(u))^2*(-2*daV_u + conj(daV_conju));

% symvar([d2aV_uconjuR, d2aV_conjuconjuR, d3aV_uuuR, d3aV_uuconjuR]) = 
%   [ aV, d2aV_uu, daV_conju, daV_u, theta, u]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aMSubTwo = -(1+u*conj(u))^2*daV_u;
T4SubTwo = ((1+u*conj(u))^2/2)*daV_conju;
bVSubTwo = u*(1+u*conj(u))*conj(daV_u) + (1+u*conj(u))^2/2*conj(d2aV_uu);
h11SubTwo = i*(1+u*conj(u))^2*aV + (i/2)*(1+u*conj(u))^4*d2aV_uconjuR;
% variableSet1 = [aM, bV, T4, h11, rho]; %last
subSet1Two =[aMSubTwo, bVSubTwo, T4SubTwo, h11SubTwo, rhoSub]; %last

CVarRicCurv2 =[u, theta, aV, daV_u, daV_conju, d2aV_uu];

myNumber2 = length(CVarRicCurv2);
dCVarRicCurv2 = sym('dCVar',[5,myNumber2]);
dCVarRicCurv2(:,1) = [0; 0; 0; 1; 0];%u
dCVarRicCurv2(:,2) = [dtheta_Mu; dtheta_conjMu; dtheta_vnormv; 0; 0]; %theta
dCVarRicCurv2(:,3) = [daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVarRicCurv2(:,4) = [d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv; d2aV_uu; d2aV_uconjuR]; %daV_u
dCVarRicCurv2(:,5) = [d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;
    d2aV_uconjuR; d2aV_conjuconjuR]; %daV_conju
dCVarRicCurv2(:,6) = [d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv; d3aV_uuuR;
    d3aV_uuconjuR]; %d2aV_uu 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d3aV_uconjuVec = df_NatTwo_MuSet_CR_rmW(d2aV_uconjuR, CVarRicCurv2, dCVarRicCurv2);
d3aV_uconjuMuR = d3aV_uconjuVec(1);
d3aV_uconjuconjMuR = d3aV_uconjuVec(2);
d3aV_uconjuvnormvR = d3aV_uconjuVec(3);
d3aV_uconjuconjuR = d3aV_uconjuVec(5); %Need to differentiate this term

d3aV_conjuconjuVec =  df_NatTwo_MuSet_CR_rmW(d2aV_conjuconjuR, CVarRicCurv2, dCVarRicCurv2);
d3aV_conjuconjuconjMuR = d3aV_conjuconjuVec(2);
d3aV_conjuconjuconjuR = d3aV_conjuconjuVec(5);

d4aV_uuconjuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uuconjuR, CVarRicCurv2, dCVarRicCurv2);
d4aV_uuconjuMuR = d4aV_uuconjuVec(1);
d4aV_uuconjuconjMuR =  d4aV_uuconjuVec(2);
d4aV_uuconjuvnormvR = d4aV_uuconjuVec(3);
d4aV_uuuconjuR = d4aV_uuconjuVec(4);
d4aV_uuconjuconjuR = d4aV_uuconjuVec(5);

d4aV_uconjuconjuVec = df_NatTwo_MuSet_CR_rmW(d3aV_uconjuconjuR, CVarRicCurv2, dCVarRicCurv2);
d4aV_uconjuconjuMuR = d4aV_uconjuconjuVec(1);
d4aV_uconjuconjuconjMuR = d4aV_uconjuconjuVec(2);
d4aV_uconjuconjuvnormvR = d4aV_uconjuconjuVec(3);
d4aV_uconjuconjuconjuR = d4aV_uconjuconjuVec(5);

% variableSet6 = symvar(RicCurv_in_aV);
variableSet6 = []; %starting 4 terms
variableSet7 = []; %other terms
subSet6 = [];
subSet7 = [];
MVarRicCurv2 = []; %simplify

myNumber2 = length(symSetRicCurv2);
for j=1:myNumber2
    myVar = symSetRicCurv2(j);
    myChar = [char(myVar), 'R'];
    if exist(myChar)==0
        eval(['MVarRicCurv2 = [MVarRicCurv2,', char(myVar), '];']);
        continue
    end 
    if strfind(myChar,'d2')==1 
        eval(['variableSet6 = [variableSet6,', char(myVar), '];']);
        eval(['subSet6 = [subSet6,', myChar, '];']);
    elseif strfind(myChar,'d3aV_uu')==1
        eval(['variableSet6 = [variableSet6,', char(myVar), '];']);
        eval(['subSet6 = [subSet6,', myChar, '];']);
    else
        eval(['variableSet7 = [variableSet7,', char(myVar), '];']);
        eval(['subSet7 = [subSet7,', myChar, '];']);
    end
end

MVarRicCurv2 = [MVarRicCurv2, Y];

for m=1:6
    for n=1:6
        temp = RicCurv_in_aV(m,n);
        temp = subs(temp, variableSet7, subSet7);
        temp = subs(temp, variableSet6, subSet6);
        temp = subs(temp, variableSet1, subSet1Two);
        % theta is real
        temp = subs(temp, conj(dtheta_vnormv), dtheta_vnormv);
        temp = subs(temp, dtheta_conjMu, conj(dtheta_Mu));
        temp = subs(temp, conj(theta), theta);
        temp = subs(temp, u*conj(u), Y-1);
        temp = complex_simple3(temp, MVarRicCurv2);
        temp = subs(temp, u*conj(u), Y-1);
        temp = complex_simple3(temp, MVarRicCurv2);      
        RicCurv_in_aV(m,n) = temp;
    end
end
clearvars m n temp

save('DataTempRicCurv2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars d2T4_uVec d2T4_conjuVec d2aM_uVec d2rho_conjuVec d3T4_uconjuVec 
clearvars d3T4_uuVec d3aV_uconjuRow  d4aV_uconjuconjuVec d4aV_uuconjuVec   
clearvars daMVec j subSet1 subSet1Two subSet2 subSet3 subSet4 subSet5
clearvars subSet6 subSet7 variableSet1 variableSet2 variableSet3
clearvars variableSet4 variableSet5 variableSet6 variableSet7

clearvars d2aV_uRow d2h11_uVec d2rho_uVec d3aV_conjuconjuVec
clearvars d3aV_uconjuVec dh11Vec dT4Vec drhoVec diffVariable2 diffVariable4

list_of_variables_RicCurv2 = who;
save('DataRicCurv2_Oct7_NatTwoR_CR_rmW.mat');
