load('Data_WeylChern_Part1_May30.mat');
clearvars WeylX WeylX0 MVarX0
assumeAlso([x y z u1 u2 gamma], 'real');
% Model w = f(x,u). So rho = 0.
variableFive_Nat_CR_rmW
% df = df_NatTwo_WeylPart5_CR_rmW(f, CVar5, dCVar5)

f = -i/2*mu1*(mu3*G11_2 + mu1*G12_3 - mu2*G11_3)...
    -i/2*mu2*(-mu3*G22_1 + mu1*G22_3 + mu2*G23_1)...
    -i/2*mu3*(mu3*G31_2 - mu1*G33_2 + mu2*G33_1);
dfdu = complexdiff3(f,u,0);
d2fdu2 = complexdiff3(dfdu,u,0);
phiF = d2fdu2 - 6*conj(u)/(1+u*conj(u))*dfdu + 12*conj(u)^2/(1+u*conj(u))^2*f;
rhoSub = i*(phiW-conj(phiW)) - i*(phiF-conj(phiF));

h11Sub = -mu1*conj(mu1)*G12_3 - mu2*conj(mu2)*G23_1 - mu3*conj(mu3)*G31_2...
    + 1/(2*norm_of_v)*(G11_2*v1*v3 - G11_3*v1*v2 - G22_1*v3*v2...
    + G22_3*v2*v1 + G33_1*v2*v3 - G33_2*v1*v3);
T4Sub = -i/(2*norm_of_v)*v1*(mu3*G11_2 - mu2*G11_3 + mu1*G12_3)...
    -i/(2*norm_of_v)*v2*(-mu3*G22_1 + mu2*G23_1 + mu1*G22_3)...
    -i/(2*norm_of_v)*v3*(mu3*G31_2 + mu2*G33_1 -mu1*G33_2);
aMuSub0 = -i*(conj(u)^2-1)*(v1*(G23_1+G31_2) + v2*G11_3 - v3*G11_2)...
    -2*i*conj(u)*(-v1*G22_3 + v2*(G12_3+G31_2) + v3*G22_1)...
    -(conj(u)^2+1)*(v1*G33_2 - v2*G33_1 + v3*(G12_3+G23_1));
aMuSub = 1/(2*(1+u*conj(u))^2)*aMuSub0;
aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);
bVSub = 1/(2*(1+u*conj(u))^2)*(-i*mu2^2*G23_1-i*mu3^2*G31_2-i*mu1^2*G12_3...
    -i*mu1*mu3*G11_2 + i*mu1*mu2*G11_3 + i*mu2*mu3*G22_1...
    -i*mu2*mu1*G22_3 + i*mu3*mu1*G33_2 - i*mu3*mu2*G33_1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dh11N = [dh11_Mu, dh11_conjMu, dh11_vnormv, dh11_u, dh11_conju];
dh11Vec = df_NatTwo_WeylPart5_CR_rmW(h11Sub, CVar5, dCVar5);

d2h11N = [d2h11_uMu, d2h11_uconjMu, d2h11_uvnormv, d2h11_uu, d2h11_uconju, ...
    d2h11_conjuMu, d2h11_conjuconjMu, d2h11_conjuvnormv, d2h11_conjuconju, ...
    d2h11_MuconjMu, d2h11_Muvnormv, d2h11_conjMuvnormv, d2h11_vnormvvnormv];

dh11_uSub = dh11Vec(4);
dh11_conjuSub = dh11Vec(5);
dh11_MuSub = dh11Vec(1);
dh11_conjMuSub = dh11Vec(2);
dh11_vnormvSub = dh11Vec(3);

d2h11_uVec = df_NatTwo_WeylPart5_CR_rmW(dh11_uSub, CVar5, dCVar5);
d2h11_conjuVec = df_NatTwo_WeylPart5_CR_rmW(dh11_conjuSub, CVar5, dCVar5);
d2h11_MuVec =  df_NatTwo_WeylPart5_CR_rmW(dh11_MuSub, CVar5, dCVar5);
d2h11_conjMuVec = df_NatTwo_WeylPart5_CR_rmW(dh11_conjMuSub, CVar5, dCVar5);
d2h11_vnormvVec = df_NatTwo_WeylPart5_CR_rmW(dh11_vnormvSub, CVar5, dCVar5);

d2h11Vec = [d2h11_uVec(1),d2h11_uVec(2),d2h11_uVec(3),d2h11_uVec(4),d2h11_uVec(5),...
    d2h11_conjuVec(1),d2h11_conjuVec(2),d2h11_conjuVec(3),d2h11_conjuVec(5),...
    d2h11_MuVec(2),d2h11_MuVec(3),d2h11_conjMuVec(3),d2h11_vnormvVec(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dT4N = [dT4_Mu, dT4_conjMu, dT4_vnormv, dT4_u, dT4_conju];
dT4Vec = df_NatTwo_WeylPart5_CR_rmW(T4Sub, CVar5, dCVar5);
%
d2T4N = [d2T4_uMu, d2T4_uconjMu, d2T4_uvnormv, d2T4_uu, d2T4_uconju, ...
    d2T4_conjuMu, d2T4_conjuconjMu, d2T4_conjuvnormv, d2T4_conjuconju, ...
    d2T4_MuMu, d2T4_MuconjMu, d2T4_Muvnormv, d2T4_conjMuconjMu, ....
    d2T4_conjMuvnormv];
    
dT4_MuSub = dT4Vec(1);
dT4_conjMuSub = dT4Vec(2);
dT4_vnormvSub = dT4Vec(3);
dT4_uSub = dT4Vec(4);
dT4_conjuSub = dT4Vec(5);
d2T4_uVec = df_NatTwo_WeylPart5_CR_rmW(dT4_uSub, CVar5, dCVar5);
d2T4_conjuVec = df_NatTwo_WeylPart5_CR_rmW(dT4_conjuSub, CVar5, dCVar5);
d2T4_MuVec = df_NatTwo_WeylPart5_CR_rmW(dT4_MuSub, CVar5, dCVar5);
d2T4_conjMuVec = df_NatTwo_WeylPart5_CR_rmW(dT4_conjMuSub, CVar5, dCVar5);
d2T4_vnormvVec = df_NatTwo_WeylPart5_CR_rmW(dT4_vnormvSub, CVar5, dCVar5);

d2T4Vec = [d2T4_uVec(1),d2T4_uVec(2),d2T4_uVec(3),d2T4_uVec(4),d2T4_uVec(5),...
    d2T4_conjuVec(1),d2T4_conjuVec(2),d2T4_conjuVec(3),d2T4_conjuVec(5),...
    d2T4_MuVec(1),d2T4_MuVec(2),d2T4_MuVec(3),d2T4_conjMuVec(2),d2T4_conjMuVec(3)];
%
d3T4N = [d3T4_uuMu, d3T4_uuconjMu, d3T4_uuvnormv, d3T4_uuu, d3T4_uuconju, ...
    d3T4_uconjuMu, d3T4_uconjuconjMu, d3T4_uconjuvnormv, d3T4_uconjuconju, ...
    d3T4_uMuMu, d3T4_uMuconjMu, d3T4_uconjMuconjMu];

d3T4_uMuVec = df_NatTwo_WeylPart5_CR_rmW(d2T4_uVec(1), CVar5, dCVar5);
d3T4_uconjMuVec = df_NatTwo_WeylPart5_CR_rmW(d2T4_uVec(2), CVar5, dCVar5);
d3T4_uuVec = df_NatTwo_WeylPart5_CR_rmW(d2T4_uVec(4), CVar5, dCVar5);
d3T4_uconjuVec = df_NatTwo_WeylPart5_CR_rmW(d2T4_uVec(5), CVar5, dCVar5);

d3T4Vec = [d3T4_uuVec(1),d3T4_uuVec(2),d3T4_uuVec(3),d3T4_uuVec(4),d3T4_uuVec(5),...
    d3T4_uconjuVec(1),d3T4_uconjuVec(2),d3T4_uconjuVec(3),d3T4_uconjuVec(5),...
    d3T4_uMuVec(1), d3T4_uMuVec(2), d3T4_uconjMuVec(2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daMuN = [daMu_Mu, daMu_conjMu, daMu_vnormv, daMu_u, daMu_conju];
daMuVec = df_NatTwo_WeylPart5_CR_rmW(aMuSub, CVar5, dCVar5);
%
d2aMuN = [d2aMu_uMu, d2aMu_uconjMu, d2aMu_uvnormv, d2aMu_uu, d2aMu_uconju, ...
    d2aMu_conjuMu, d2aMu_conjuconjMu, d2aMu_conjuvnormv, d2aMu_conjuconju, ...
    d2aMu_MuMu, d2aMu_MuconjMu, d2aMu_Muvnormv, d2aMu_conjMuvnormv];

d2aMu_uVec = df_NatTwo_WeylPart5_CR_rmW(daMuVec(4), CVar5, dCVar5);
d2aMu_conjuVec = df_NatTwo_WeylPart5_CR_rmW(daMuVec(5), CVar5, dCVar5);
d2aMu_MuVec = df_NatTwo_WeylPart5_CR_rmW(daMuVec(1), CVar5, dCVar5);
d2aMu_conjMuVec = df_NatTwo_WeylPart5_CR_rmW(daMuVec(2), CVar5, dCVar5);

d2aMuVec = [d2aMu_uVec(1),d2aMu_uVec(2),d2aMu_uVec(3),d2aMu_uVec(4),d2aMu_uVec(5),...
    d2aMu_conjuVec(1),d2aMu_conjuVec(2),d2aMu_conjuVec(3),d2aMu_conjuVec(5),...
    d2aMu_MuVec(1), d2aMu_MuVec(2), d2aMu_MuVec(3), d2aMu_conjMuVec(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daVN = [daV_Mu, daV_conjMu, daV_vnormv, daV_u, daV_conju];
daVVec = df_NatTwo_WeylPart5_CR_rmW(aVSub, CVar5, dCVar5);
%
d2aVN = [d2aV_uMu, d2aV_uconjMu, d2aV_uvnormv, d2aV_uu, d2aV_uconju, ...
    d2aV_conjuMu, d2aV_conjuconjMu, d2aV_conjuvnormv, d2aV_conjuconju, ...
    d2aV_MuMu, d2aV_MuconjMu, d2aV_conjMuconjMu];

d2aV_uVec = df_NatTwo_WeylPart5_CR_rmW(daVVec(4), CVar5, dCVar5);
d2aV_conjuVec = df_NatTwo_WeylPart5_CR_rmW(daVVec(5), CVar5, dCVar5);
d2aV_MuVec = df_NatTwo_WeylPart5_CR_rmW(daVVec(1), CVar5, dCVar5);
d2aV_conjMuVec = df_NatTwo_WeylPart5_CR_rmW(daVVec(2), CVar5, dCVar5);

d2aVVec = [d2aV_uVec(1),d2aV_uVec(2),d2aV_uVec(3),d2aV_uVec(4),d2aV_uVec(5),...
    d2aV_conjuVec(1),d2aV_conjuVec(2),d2aV_conjuVec(3),d2aV_conjuVec(5),...
    d2aV_MuVec(1), d2aV_MuVec(2), d2aV_conjMuVec(2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbVN = [dbV_Mu, dbV_conjMu, dbV_vnormv, dbV_u, dbV_conju];
dbVVec = df_NatTwo_WeylPart5_CR_rmW(bVSub, CVar5, dCVar5);
%
drhoN = [drho_Mu, drho_conjMu, drho_vnormv, drho_u, drho_conju];
drhoVec =  df_NatTwo_WeylPart5_CR_rmW(rhoSub, CVar5, dCVar5);
%
d2rhoN = [d2rho_uMu,d2rho_uconjMu,d2rho_uvnormv,d2rho_uu,d2rho_uconju,...
    d2rho_conjuMu, d2rho_conjuconjMu, d2rho_conjuvnormv, d2rho_conjuconju,...
    d2rho_MuMu, d2rho_MuconjMu, d2rho_conjMuconjMu];
d2rho_uVec = df_NatTwo_WeylPart5_CR_rmW(drhoVec(4), CVar5, dCVar5);
d2rho_conjuVec = df_NatTwo_WeylPart5_CR_rmW(drhoVec(5), CVar5, dCVar5);
d2rho_MuVec = df_NatTwo_WeylPart5_CR_rmW(drhoVec(1), CVar5, dCVar5);
d2rho_conjMuVec = df_NatTwo_WeylPart5_CR_rmW(drhoVec(2), CVar5, dCVar5);

d2rhoVec = [d2rho_uVec(1),d2rho_uVec(2),d2rho_uVec(3),d2rho_uVec(4),...
    d2rho_uVec(5),...
    d2rho_conjuVec(1),d2rho_conjuVec(2),d2rho_conjuVec(3),d2rho_conjuVec(5),...
    d2rho_MuVec(1),d2rho_MuVec(2),d2rho_conjMuVec(2)];

MVarFive = [G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3, G33_1, G33_2,...
    dG12_3_Mu, dG12_3_conjMu, dG12_3_vnormv,...
    dG23_1_Mu, dG23_1_conjMu, dG23_1_vnormv,...
    dG31_2_Mu, dG31_2_conjMu, dG31_2_vnormv,...
    dG11_2_Mu, dG11_2_conjMu, dG11_2_vnormv,...
    dG11_3_Mu, dG11_3_conjMu, dG11_3_vnormv,...
    dG22_1_Mu, dG22_1_conjMu, dG22_1_vnormv,...
    dG22_3_Mu, dG22_3_conjMu, dG22_3_vnormv,...
    dG33_1_Mu, dG33_1_conjMu, dG33_1_vnormv,...
    dG33_2_Mu, dG33_2_conjMu, dG33_2_vnormv,...
    d2G12_3_MuMu, d2G12_3_MuconjMu, d2G12_3_Muvnormv,...
    d2G12_3_conjMuconjMu, d2G12_3_conjMuvnormv, d2G12_3_vnormvvnormv,...
    d2G23_1_MuMu, d2G23_1_MuconjMu, d2G23_1_Muvnormv,...
    d2G23_1_conjMuconjMu, d2G23_1_conjMuvnormv, d2G23_1_vnormvvnormv,...
    d2G31_2_MuMu, d2G31_2_MuconjMu, d2G31_2_Muvnormv,...
    d2G31_2_conjMuconjMu, d2G31_2_conjMuvnormv, d2G31_2_vnormvvnormv,...
    d2G11_2_MuMu, d2G11_2_MuconjMu, d2G11_2_Muvnormv,...
    d2G11_2_conjMuconjMu, d2G11_2_conjMuvnormv, d2G11_2_vnormvvnormv,...
    d2G11_3_MuMu, d2G11_3_MuconjMu, d2G11_3_Muvnormv,...
    d2G11_3_conjMuconjMu, d2G11_3_conjMuvnormv, d2G11_3_vnormvvnormv,...
    d2G22_1_MuMu, d2G22_1_MuconjMu, d2G22_1_Muvnormv,...
    d2G22_1_conjMuconjMu, d2G22_1_conjMuvnormv, d2G22_1_vnormvvnormv,...
    d2G22_3_MuMu, d2G22_3_MuconjMu, d2G22_3_Muvnormv,...
    d2G22_3_conjMuconjMu, d2G22_3_conjMuvnormv, d2G22_3_vnormvvnormv,...
    d2G33_1_MuMu, d2G33_1_MuconjMu, d2G33_1_Muvnormv,...
    d2G33_1_conjMuconjMu, d2G33_1_conjMuvnormv, d2G33_1_vnormvvnormv,...
    d2G33_2_MuMu, d2G33_2_MuconjMu, d2G33_2_Muvnormv,...
    d2G33_2_conjMuconjMu, d2G33_2_conjMuvnormv, d2G33_2_vnormvvnormv,...
    w, dw_Mu, dw_conjMu, dw_vnormv, dw_u, ...
    d2w_MuMu, d2w_MuconjMu, d2w_Muvnormv, d2w_conjMuconjMu, d2w_conjMuvnormv,...
    d2w_vnormvvnormv, d2w_uMu, d2w_uconjMu, d2w_uvnormv, d2w_uu];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countSet225 = [];
for m=1:6
    for n=1:6
        for k=1:6
            for ll=1:6
                rowTemp = [m,n,k,ll];
                if (m<n)&&(k<ll)
                    countSet225 = [countSet225; rowTemp];
                end
            end
        end
    end
end
numberCount = zeros(225,2);
countSetTwo = []; % length(countSetTwo)=120
for j=1:225
    numberCount(j,1) = countSet225(j,1)*10 + countSet225(j,2);
    numberCount(j,2) = countSet225(j,3)*10 + countSet225(j,4);
    if numberCount(j,1)<=numberCount(j,2)
        countSetTwo = [countSetTwo; countSet225(j,:)];
    end      
end
clearvars countSet225 numberCount j m n k ll
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WeylF = sym('WeylF',[120 5]);
for j=1:120
    m = countSetTwo(j,1);
    n = countSetTwo(j,2);
    k = countSetTwo(j,3);
    ll = countSetTwo(j,4);
    temp = Weyl(m,n,k,ll);
    WeylF(j,:) = [m, n, k, ll, temp];
end
% symvarWeyl

clearvars m n k ll temp
for j=1:120
    temp = WeylF(j,5);
    temp = subs(temp, dh11N, dh11Vec);
    temp = subs(temp, d2h11N, d2h11Vec);
    temp = subs(temp, dT4N, dT4Vec); 
    temp = subs(temp, d2T4N, d2T4Vec);
    temp = subs(temp, d3T4N, d3T4Vec);
    temp = subs(temp, daMuN, daMuVec);
    temp = subs(temp, d2aMuN, d2aMuVec);
    temp = subs(temp, daVN, daVVec);
    temp = subs(temp, d2aVN, d2aVVec);
    temp = subs(temp, dbVN, dbVVec);
    temp = subs(temp, drhoN, drhoVec);
    temp = subs(temp, d2rhoN, d2rhoVec);
    temp = subs(temp, [h11, T4, aMu, aV, bV, rho],...
        [h11Sub, T4Sub, aMuSub, aVSub, bVSub, rhoSub]);
    % temp = complex_simple3(temp, MVarFive);
    WeylF(j,5) = temp;
end
clearvars temp j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('DataZero_WeylChern_Part5_Jun10.mat');


%%
load('DataZero_WeylChern_Part5_Jun10.mat');
assumeAlso([x y z u1 u2 gamma], 'real');

remove_variable_Jun10

save('Data_WeylChern_Part5_Jun10.mat');
