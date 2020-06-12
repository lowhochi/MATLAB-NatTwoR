load('DataMain4_NatTwoR_CR_rmW.mat');
syms theta % theta = G12_3 + G23_1 + G31_2;
% symSetWeyl = symbolic variables of Weyl(m,n,k,ll). length = 129
% Replace Weyl curvature tensor by aV and its derivatives.
variableSet1 = [aM, bV, T4, h11, rho];

wSet = [u, w, dw_u, dw_Mu, dw_conjMu, dw_vnormv,...
    d2w_uu, d2w_uMu, d3w_uuu, d3w_uuMu, d2w_uconjMu, d2w_uvnormv,...
    d2w_MuconjMu, d3w_uuconjMu, d2w_Muvnormv, d3w_uuvnormv, d3w_uMuconjMu,...
    d3w_uMuvnormv, d2w_conjMuconjMu, d2w_conjMuvnormv, d3w_uconjMuconjMu,...
    d2w_vnormvvnormv, d3w_uconjMuvnormv];
  
aVSet = [aV, d2aV_uu, d2aV_uMu, d2aV_MuMu, d2aV_uconju, d2aV_conjuMu,...
    d2aV_uconjMu, d2aV_MuconjMu, d2aV_conjuconju, d2aV_conjuconjMu, d2aV_conjMuconjMu,...
    daV_u, daV_Mu, daV_conju, daV_conjMu, daV_vnormv];

variable_NatTwoR_Weyl1_CR_rmW

WeylTwo = sym('Weyl',[120,5]);
for j=1:120
    m = countIndex120(j,1);
    n = countIndex120(j,2);
    k = countIndex120(j,3);
    ll = countIndex120(j,4);
    temp = Weyl(m,n,k,ll);
    WeylTwo(j,:) = [m,n,k,ll,temp];
end

% comment the below session in actual running.
% checkSet1 = union(variableSet1, wSet);
% checkSet1 = union(checkSet1, aVSet);
% variableSet2 =[]; %first derivative
% variableSet3 =[]; %second derivative
% variableSet4 = []; %third derivative
%   
% for j=1:129
%     myVar = symSetWeyl(j);
%     myStr = string(myVar);
%     if (ismember(myVar,checkSet1)==1)
%         continue
%     end   
%     if (strfind(myStr,'d2')==1)
%         variableSet3 = [variableSet3, myVar];
%     elseif (strfind(myStr,'d3')==1)
%         variableSet4 = [variableSet4, myVar];
%     else
%         variableSet2 = [variableSet2, myVar];
%     end    
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variableSet1 = [aM, bV, T4, h11, rho];
aMSub = -(1+u*conj(u))^2*daV_u;
T4Sub = ((1+u*conj(u))^2/2)*daV_conju;
bVSub = u*(1+u*conj(u))*conj(daV_u) + (1+u*conj(u))^2/2*conj(d2aV_uu);
h11Sub = i*(1+u*conj(u))^2*aV + (i/2)*(1+u*conj(u))^4*d2aV_uconju;

phiW = d2w_uu - 6*conj(u)/(1+u*conj(u))*dw_u + 12*conj(u)^2/(1+u*conj(u))^2*w;
rhoSub = i*(phiW-conj(phiW)) + 16*theta + 12*i*(aV-conj(aV));

subSet1 = [aMSub,bVSub, T4Sub, h11Sub, rhoSub];
% symvar(subSet1) = [aV, d2aV_uconju, d2aV_uu, d2w_uu, daV_conju,...
%   daV_u, dw_u, theta, u, w]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSet2 =[dT4_u, dT4_Mu, dT4_conju, dT4_conjMu, dT4_vnormv,...
    daM_u, daM_Mu, daM_conju, daM_conjMu, daM_vnormv,...
    dbV_u, dbV_Mu, dbV_conju, dbV_conjMu, dbV_vnormv,...
    dh11_u, dh11_Mu, dh11_conju, dh11_conjMu, dh11_vnormv,...
    drho_u, drho_Mu, drho_conju, drho_conjMu, drho_vnormv];

daMVec = df_NatTwo_MuSet_CR_rmW(aMSub, CVarWeyl, dCVarWeyl); %ok
daM_MuSub = daMVec(1);
daM_conjMuSub = daMVec(2);
daM_vnormvSub = daMVec(3);
daM_uSub = daMVec(4);
daM_conjuSub = daMVec(5);

dT4Vec = df_NatTwo_MuSet_CR_rmW(T4Sub, CVarWeyl, dCVarWeyl); %ok
dT4_MuSub = dT4Vec(1);
dT4_conjMuSub = dT4Vec(2);
dT4_vnormvSub = dT4Vec(3);
dT4_uSub = dT4Vec(4);
dT4_conjuSub = dT4Vec(5);

dbVVec = df_NatTwo_MuSet_CR_rmW(bVSub, CVarWeyl, dCVarWeyl); %ok
dbV_MuSub = dbVVec(1);
dbV_conjMuSub = dbVVec(2);
dbV_vnormvSub = dbVVec(3);
dbV_uSub = dbVVec(4);
dbV_conjuSub = dbVVec(5);

dh11Vec = df_NatTwo_MuSet_CR_rmW(h11Sub, CVarWeyl, dCVarWeyl); %ok
dh11_MuSub = dh11Vec(1);
dh11_conjMuSub = dh11Vec(2);
dh11_vnormvSub = dh11Vec(3);
dh11_uSub = dh11Vec(4);
dh11_conjuSub = dh11Vec(5);

drhoVec = df_NatTwo_MuSet_CR_rmW(rhoSub, CVarWeyl, dCVarWeyl);%ok
drho_MuSub = drhoVec(1);
drho_conjMuSub = drhoVec(2);
drho_vnormvSub = drhoVec(3);
drho_uSub = drhoVec(4);
drho_conjuSub = drhoVec(5);

subSet2 = [dT4_uSub, dT4_MuSub, dT4_conjuSub, dT4_conjMuSub, dT4_vnormvSub,...
     daM_uSub, daM_MuSub, daM_conjuSub, daM_conjMuSub, daM_vnormvSub,...
     dbV_uSub, dbV_MuSub, dbV_conjuSub, dbV_conjMuSub, dbV_vnormvSub,...
     dh11_uSub, dh11_MuSub, dh11_conjuSub, dh11_conjMuSub, dh11_vnormvSub,...
     drho_uSub, drho_MuSub, drho_conjuSub, drho_conjMuSub, drho_vnormvSub];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffVariable2 = symvar([dT4_uSub, dT4_conjuSub, dT4_MuSub, dT4_conjMuSub,...
    daM_uSub, daM_conjuSub, daM_MuSub, daM_conjMuSub,...
    dh11_uSub, dh11_conjuSub, dh11_MuSub, dh11_conjMuSub, dh11_vnormvSub,...
    drho_uSub, drho_conjuSub, drho_MuSub, drho_conjMuSub]);

variableSet3 = [d2T4_uu, d2T4_uMu, d2T4_MuMu, d2T4_uconju, d2T4_conjuMu,...
    d2T4_uconjMu, d2T4_uvnormv, d2T4_MuconjMu, d2T4_Muvnormv, d2T4_conjuconju,...
    d2T4_conjuconjMu, d2T4_conjuvnormv, d2T4_conjMuconjMu, d2T4_conjMuvnormv,...
    d2aM_uu, d2aM_uMu, d2aM_MuMu, d2aM_uconju, d2aM_conjuMu, d2aM_uconjMu,...
    d2aM_uvnormv, d2aM_MuconjMu, d2aM_Muvnormv, d2aM_conjuconju, d2aM_conjuconjMu,...
    d2aM_conjuvnormv, d2aM_conjMuvnormv, ...
    d2h11_uu, d2h11_uMu, d2h11_uconju, d2h11_conjuMu, d2h11_uconjMu, d2h11_uvnormv,...
    d2h11_MuconjMu, d2h11_Muvnormv, d2h11_conjuconju, d2h11_conjuconjMu,...
    d2h11_conjuvnormv, d2h11_conjMuvnormv, d2h11_vnormvvnormv,...
    d2rho_uu, d2rho_uMu, d2rho_MuMu, d2rho_uconju,...
    d2rho_conjuMu, d2rho_uconjMu, d2rho_MuconjMu, d2rho_conjuconju,...
    d2rho_conjuconjMu, d2rho_conjMuconjMu]; 
%   length=50;
 
d2T4_uVec = df_NatTwo_MuSet_CR_rmW(dT4_uSub, CVarWeyl, dCVarWeyl); %ok
d2T4_uMuSub = d2T4_uVec(1);
d2T4_uconjMuSub = d2T4_uVec(2);
d2T4_uvnormvSub = d2T4_uVec(3);
d2T4_uuSub = d2T4_uVec(4);
d2T4_uconjuSub = d2T4_uVec(5);
d2T4_conjuVec = df_NatTwo_MuSet_CR_rmW(dT4_conjuSub, CVarWeyl, dCVarWeyl); %ok
d2T4_conjuMuSub = d2T4_conjuVec(1);
d2T4_conjuconjMuSub = d2T4_conjuVec(2);
d2T4_conjuvnormvSub = d2T4_conjuVec(3);
d2T4_conjuconjuSub = d2T4_conjuVec(5);
d2T4_MuVec = df_NatTwo_MuSet_CR_rmW(dT4_MuSub, CVarWeyl, dCVarWeyl); %ok
d2T4_MuMuSub = d2T4_MuVec(1);
d2T4_MuconjMuSub = d2T4_MuVec(2);
d2T4_MuvnormvSub = d2T4_MuVec(3);
d2T4_conjMuVec = df_NatTwo_MuSet_CR_rmW(dT4_conjMuSub, CVarWeyl, dCVarWeyl); %ok
d2T4_conjMuconjMuSub = d2T4_conjMuVec(2);
d2T4_conjMuvnormvSub = d2T4_conjMuVec(3);

d2aM_uVec = df_NatTwo_MuSet_CR_rmW(daM_uSub, CVarWeyl, dCVarWeyl); %ok
d2aM_uMuSub = d2aM_uVec(1);
d2aM_uconjMuSub = d2aM_uVec(2); 
d2aM_uvnormvSub = d2aM_uVec(3);
d2aM_uuSub = d2aM_uVec(4);
d2aM_uconjuSub = d2aM_uVec(5);
d2aM_conjuVec = df_NatTwo_MuSet_CR_rmW(daM_conjuSub, CVarWeyl, dCVarWeyl); %ok
d2aM_conjuMuSub = d2aM_conjuVec(1);
d2aM_conjuconjMuSub = d2aM_conjuVec(2);
d2aM_conjuvnormvSub = d2aM_conjuVec(3);
d2aM_conjuconjuSub = d2aM_conjuVec(5);
d2aM_MuVec = df_NatTwo_MuSet_CR_rmW(daM_MuSub, CVarWeyl, dCVarWeyl);
d2aM_conjMuVec = df_NatTwo_MuSet_CR_rmW(daM_conjMuSub, CVarWeyl, dCVarWeyl); %ok
d2aM_MuMuSub = d2aM_MuVec(1);
d2aM_MuconjMuSub = d2aM_MuVec(2);
d2aM_MuvnormvSub = d2aM_MuVec(3);
d2aM_conjMuvnormvSub = d2aM_conjMuVec(3);

d2h11_uVec = df_NatTwo_MuSet_CR_rmW(dh11_uSub, CVarWeyl, dCVarWeyl); %ok
d2h11_uMuSub = d2h11_uVec(1);
d2h11_uconjMuSub = d2h11_uVec(2);
d2h11_uvnormvSub = d2h11_uVec(3);
d2h11_uuSub = d2h11_uVec(4);
d2h11_uconjuSub = d2h11_uVec(5);
d2h11_conjuVec = df_NatTwo_MuSet_CR_rmW(dh11_conjuSub, CVarWeyl, dCVarWeyl); %ok
d2h11_conjuMuSub = d2h11_conjuVec(1) ;
d2h11_conjuconjMuSub = d2h11_conjuVec(2);
d2h11_conjuvnormvSub = d2h11_conjuVec(3);
d2h11_conjuconjuSub = d2h11_conjuVec(5);
d2h11_MuVec = df_NatTwo_MuSet_CR_rmW(dh11_MuSub, CVarWeyl, dCVarWeyl); %ok
d2h11_conjMuVec = df_NatTwo_MuSet_CR_rmW(dh11_conjMuSub, CVarWeyl, dCVarWeyl); %ok
d2h11_vnormvVec = df_NatTwo_MuSet_CR_rmW(dh11_vnormvSub, CVarWeyl, dCVarWeyl); %ok
d2h11_MuconjMuSub = d2h11_MuVec(2);
d2h11_MuvnormvSub = d2h11_MuVec(3);
d2h11_conjMuvnormvSub = d2h11_conjMuVec(3);
d2h11_vnormvvnormvSub = d2h11_vnormvVec(3);

d2rho_uVec = df_NatTwo_MuSet_CR_rmW(drho_uSub, CVarWeyl, dCVarWeyl); %ok
d2rho_uMuSub = d2rho_uVec(1);
d2rho_uconjMuSub = d2rho_uVec(2);
d2rho_uuSub = d2rho_uVec(4);
d2rho_uconjuSub = d2rho_uVec(5);

d2rho_conjuVec = df_NatTwo_MuSet_CR_rmW(drho_conjuSub, CVarWeyl, dCVarWeyl); %ok
d2rho_conjuMuSub = d2rho_conjuVec(1);
d2rho_conjuconjMuSub = d2rho_conjuVec(2);
d2rho_conjuconjuSub = d2rho_conjuVec(5);

d2rho_MuVec = df_NatTwo_MuSet_CR_rmW(drho_MuSub, CVarWeyl, dCVarWeyl); %ok
d2rho_conjMuVec = df_NatTwo_MuSet_CR_rmW(drho_conjMuSub, CVarWeyl, dCVarWeyl); %ok
d2rho_MuMuSub = d2rho_MuVec(1);
d2rho_MuconjMuSub = d2rho_MuVec(2);
d2rho_conjMuconjMuSub = d2rho_conjMuVec(2);

subSet3 = [d2T4_uuSub, d2T4_uMuSub, d2T4_MuMuSub, d2T4_uconjuSub, d2T4_conjuMuSub,...
    d2T4_uconjMuSub, d2T4_uvnormvSub, d2T4_MuconjMuSub, d2T4_MuvnormvSub, d2T4_conjuconjuSub,...
    d2T4_conjuconjMuSub, d2T4_conjuvnormvSub, d2T4_conjMuconjMuSub, d2T4_conjMuvnormvSub,...
    d2aM_uuSub, d2aM_uMuSub, d2aM_MuMuSub, d2aM_uconjuSub, d2aM_conjuMuSub, d2aM_uconjMuSub,...
    d2aM_uvnormvSub, d2aM_MuconjMuSub, d2aM_MuvnormvSub, d2aM_conjuconjuSub, d2aM_conjuconjMuSub,...
    d2aM_conjuvnormvSub, d2aM_conjMuvnormvSub, ...
    d2h11_uuSub, d2h11_uMuSub, d2h11_uconjuSub, d2h11_conjuMuSub, d2h11_uconjMuSub, d2h11_uvnormvSub,...
    d2h11_MuconjMuSub, d2h11_MuvnormvSub, d2h11_conjuconjuSub, d2h11_conjuconjMuSub,...
    d2h11_conjuvnormvSub, d2h11_conjMuvnormvSub, d2h11_vnormvvnormvSub,...
    d2rho_uuSub, d2rho_uMuSub, d2rho_MuMuSub, d2rho_uconjuSub,...
    d2rho_conjuMuSub, d2rho_uconjMuSub, d2rho_MuconjMuSub, d2rho_conjuconjuSub,...
    d2rho_conjuconjMuSub, d2rho_conjMuconjMuSub]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diffVariable3 = symvar([d2T4_uuSub, d2T4_uconjuSub, d2T4_uMuSub, d2T4_uconjMuSub]);
 
variableSet4 = [d3T4_uuu, d3T4_uuMu, d3T4_uMuMu, d3T4_uuconju,...
    d3T4_uconjuMu, d3T4_uuconjMu, d3T4_uMuconjMu, d3T4_uconjuconju,...
    d3T4_uconjuconjMu, d3T4_uconjMuconjMu];

d3T4_uuVec = df_NatTwo_MuSet_CR_rmW(d2T4_uuSub, CVarWeyl, dCVarWeyl); %ok
d3T4_uuMuSub = d3T4_uuVec(1);
d3T4_uuconjMuSub = d3T4_uuVec(2);
d3T4_uuuSub = d3T4_uuVec(4);
d3T4_uuconjuSub = d3T4_uuVec(5);
d3T4_uconjuVec = df_NatTwo_MuSet_CR_rmW(d2T4_uconjuSub, CVarWeyl, dCVarWeyl); %ok
d3T4_uconjuMuSub = d3T4_uconjuVec(1);
d3T4_uconjuconjMuSub = d3T4_uconjuVec(2);
d3T4_uconjuconjuSub = d3T4_uconjuVec(5);

d3T4_uMuVec = df_NatTwo_MuSet_CR_rmW(d2T4_uMuSub, CVarWeyl, dCVarWeyl); %ok
d3T4_uconjMuVec = df_NatTwo_MuSet_CR_rmW(d2T4_uconjMuSub, CVarWeyl, dCVarWeyl); %ok
d3T4_uMuMuSub = d3T4_uMuVec(1);
d3T4_uMuconjMuSub = d3T4_uMuVec(2);
d3T4_uconjMuconjMuSub = d3T4_uconjMuVec(2);

subSet4 = [d3T4_uuuSub, d3T4_uuMuSub, d3T4_uMuMuSub, d3T4_uuconjuSub,...
    d3T4_uconjuMuSub, d3T4_uuconjMuSub, d3T4_uMuconjMuSub, d3T4_uconjuconjuSub,...
    d3T4_uconjuconjMuSub, d3T4_uconjMuconjMuSub];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:120
    m = WeylTwo(j,1);
    n = WeylTwo(j,2);
    k = WeylTwo(j,3);
    ll =WeylTwo(j,4);
    temp = WeylTwo(j,5);
    
    for k=1:length(variableSet4)
        myChar = char(variableSet4(k));
        myCharSub = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar,',',myCharSub,');']);
    end
    
    for k=1:length(variableSet3)
        myChar = char(variableSet3(k));
        myCharSub = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar,',',myCharSub,');']);
    end
    
    for k=1:length(variableSet2)
        myChar = char(variableSet2(k));
        myCharSub = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar,',',myCharSub,');']);
    end
    
    for k=1:length(variableSet1)
        myChar = char(variableSet1(k));
        myCharSub = [myChar,'Sub'];
        eval(['temp=subs(temp,', myChar,',',myCharSub,');']);
    end
    
%   temp = complex_simple3(temp, MVarWeylOne);
    WeylTwo(j,5)=temp;
end

clearvars j m n k ll temp 
clearvars myChar myCharSub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('DataWeyl_NatTwoR_CR_rmW.mat');
