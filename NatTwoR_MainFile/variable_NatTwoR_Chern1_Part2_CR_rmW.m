% variable_NatTwoR_Chern1_Part2_CR_rmW.m
% Replacement of derivatives is from the highest to the lowest.
% We may create Lie brackets between Mu, conjMu and vnormv in the process,
% so we redefine  aM T4 and bV for final substitution.

aMSubTwo = -(1+u*conj(u))^2*daV_u;
T4SubTwo = ((1+u*conj(u))^2/2)*daV_conju;
bVSubTwo = u*(1+u*conj(u))*conj(daV_u) + (1+u*conj(u))^2/2*conj(d2aV_uu);
h11SubTwo = i*(1+u*conj(u))^2*aV + (i/2)*(1+u*conj(u))^4*d2aV_uconjuR;
% variableSet1 = [aMu, bV, T4, h11]; %last
subSet1Two =[aMSubTwo, bVSubTwo, T4SubTwo, h11SubTwo]; %last

% Derivatives of theta
syms dtheta_Mu 
syms dtheta_vnormv real
dtheta_conjMu = conj(dtheta_Mu); %theta is real
% The second derivatives
syms d2theta_MuMu d2theta_MuconjMu d2theta_Muvnormv
syms d2theta_vnormvvnormv real

d2theta_conjMuconjMu = conj(d2theta_MuMu); %theta is real
d2theta_conjMuMu = conj(d2theta_MuconjMu); %theta is real
d2theta_conjMuvnormv = conj(d2theta_Muvnormv); %theta is real

dthetaRow = [dtheta_Mu, dtheta_conjMu, dtheta_vnormv];
d2theta_vnormvMu = d2theta_Muvnormv + dthetaRow*lieMu{1,3}; %it has T4 aV bV h11
d2theta_vnormvconjMu = conj(d2theta_vnormvMu);  %theta is real
d2theta_Muu = dthetaRow*lieMu{4,1};
d2theta_conjMuconju = conj(d2theta_Muu); %theta is real
d2theta_vnormvu = dthetaRow*lieMu{4,3};
d2theta_vnormvconju = conj(d2theta_vnormvu);  %theta is real

CVarChern1S = sym('CVarChern1S',[1,11]);
CVarChern1S(1) = u;
CVarChern1S(2) = theta;
CVarChern1S(3) = aV;
CVarChern1S(4) = daV_u;
CVarChern1S(5) = daV_conju;
CVarChern1S(6) = d2aV_uu;
CVarChern1S(7) = daV_Mu;
CVarChern1S(8) = daV_conjMu;
CVarChern1S(9) = daV_vnormv;
CVarChern1S(10) = dtheta_Mu;
CVarChern1S(11) = dtheta_vnormv;

dCVarChern1S = sym('dCVarChern1S',[5,11]);
% We use the -R version of a derivative of aV if it can be simplified:
% d2aV_uconjuR, d2aV_conjuconjuR, d3aV_uuuR, d3aV_uuconjuR

dCVarChern1S(:,1) = [0; 0; 0; 1; 0];%u
dCVarChern1S(:,2) = [dtheta_Mu; dtheta_conjMu; dtheta_vnormv; 0; 0]; %theta
dCVarChern1S(:,3) = [daV_Mu; daV_conjMu; daV_vnormv; daV_u; daV_conju]; %aV
dCVarChern1S(:,4) = [d2aV_uMu; d2aV_uconjMu; d2aV_uvnormv; d2aV_uu; d2aV_uconjuR]; %daV_u
dCVarChern1S(:,5) = [d2aV_conjuMu; d2aV_conjuconjMu; d2aV_conjuvnormv;
    d2aV_uconjuR; d2aV_conjuconjuR]; %daV_conju
dCVarChern1S(:,6) = [d3aV_uuMu; d3aV_uuconjMu; d3aV_uuvnormv; d3aV_uuuR;
    d3aV_uuconjuR]; %d2aV_uu 
dCVarChern1S(:,7) = [d2aV_MuMu; d2aV_MuconjMu; d2aV_Muvnormv; d2aV_Muu; d2aV_Muconju];%daV_Mu;
dCVarChern1S(:,8) = [d2aV_conjMuMu; d2aV_conjMuconjMu; d2aV_conjMuvnormv;
    d2aV_conjMuu; d2aV_conjMuconju]; %daV_conjMu;
dCVarChern1S(:,9) = [d2aV_vnormvMu; d2aV_vnormvconjMu; d2aV_vnormvvnormv;...
    d2aV_vnormvu; d2aV_vnormvconju]; %daV_vnormv;
dCVarChern1S(:,10) = [d2theta_MuMu; d2theta_MuconjMu; d2theta_Muvnormv;
     d2theta_Muu; 0];%dtheta_Mu;
dCVarChern1S(:,11) = [d2theta_vnormvMu; d2theta_vnormvconjMu; d2theta_vnormvvnormv;
    d2theta_vnormvu; d2theta_vnormvconju];%dtheta_vnormv;

MVarChern1S = [u, aV, daV_Mu, daV_conjMu, daV_vnormv, daV_u, daV_conju,...
    d2aV_MuMu, d2aV_MuconjMu, d2aV_Muvnormv, d2aV_conjMuconjMu,...
    d2aV_conjMuvnormv, d2aV_vnormvvnormv, d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv,...
    d2aV_uu, d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv,... 
    d3aV_uMuMu, d3aV_uMuconjMu, d3aV_uMuvnormv, d3aV_uconjMuconjMu,... 
    d3aV_uconjMuvnormv, d3aV_uvnormvvnormv, d3aV_uuMu, d3aV_uuconjMu, d3aV_uuvnormv,...
    dtheta_Mu, d2theta_MuMu, d2theta_MuconjMu, d2theta_Muvnormv,...
    d2w_uu, d2w_uconjMu, d2w_conjMuconjMu, dw_u, dw_conjMu, dw_vnormv, w];