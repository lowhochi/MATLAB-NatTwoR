load('DataMain4_NatTwoR_CR_rmW.mat');
load('ricciData.mat');
% ricci = ricci_curv{1};

% RicCurv(m,n) = Ric(u_m,u_n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSetRic = symvar(RicCurv);
subSetRic = variableSetRic;
RicFlat = sym('Ric',[6 6]);
for j=1:86
    myVar = variableSetRic(j);
    myChar = char(myVar);
    if strfind(myChar,'dw')==1
       subSetRic(j)=myVar;
   elseif strfind(myChar,'d2w')==1
       subSetRic(j)=myVar;
   elseif strfind(myChar,'d3w')==1
       subSetRic(j)=myVar;
   elseif strfind(myChar,'d4w')==1
       subSetRic(j)=myVar;
   elseif strfind(myChar, 'd2rho')==1
       subSetRic(j)=myVar;
   elseif strfind(myChar, 'drho')==1
       subSetRic(j)=myVar;
   elseif myVar==u
       subSetRic(j)=myVar;
   elseif myVar==w
       subSetRic(j)=myVar;
   elseif myVar==rho
       subSetRic(j)=myVar;
   else
       subSetRic(j)=0;
       continue
   end
end
for m=1:6
    for n=1:6
        temp=RicCurv(m,n);
        temp=subs(temp, variableSetRic, subSetRic);
        RicFlat(m,n)=temp;
    end
end
clearvars j m n temp myChar myVar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms d4w_uuuMu d4w_uuuu
phiW = d2w_uu-6*conj(u)*dw_u/(1+u*conj(u))+12*conj(u)^2*w/(1+u*conj(u))^2;
dPhiW_Mu = d3w_uuMu -6*conj(u)*d2w_uMu/(1+u*conj(u))...
    +12*conj(u)^2*dw_Mu/(1+u*conj(u))^2;
dPhiW_conjMu = d3w_uuconjMu -6*conj(u)*d2w_uconjMu/(1+u*conj(u))...
    +12*conj(u)^2*dw_conjMu/(1+u*conj(u))^2;
dPhiW_vnormv = d3w_uuvnormv -6*conj(u)*d2w_uvnormv/(1+u*conj(u))...
    +12*conj(u)^2*dw_vnormv/(1+u*conj(u))^2;

rhoSub = i*phiW - i*conj(phiW);
drho_MuSub = i*dPhiW_Mu -i*conj(dPhiW_conjMu);
drho_conjMuSub = i*dPhiW_conjMu -i*conj(dPhiW_Mu);
drho_vnormvSub = i*dPhiW_vnormv -i*conj(dPhiW_vnormv);

drho_uSub = i*d3w_uuu -6*i*conj(u)/(1+u*conj(u))*d2w_uu...
    +18*i*conj(u)^2/(1+u*conj(u))^2*dw_u...
    -24*i*conj(u)^3/(1+u*conj(u))^3*w...
    +6*i/(1+u*conj(u))^2*conj(dw_u) -24*i*u/(1+u*conj(u))^3*conj(w);
drho_conjuSub = conj(drho_uSub);
d2rho_uMuSub = i*d4w_uuuMu -6*i*conj(u)/(1+u*conj(u))*d3w_uuMu...
    +18*i*conj(u)^2/(1+u*conj(u))^2*d2w_uMu...
    -24*i*conj(u)^3/(1+u*conj(u))^3*dw_Mu...
    +6*i/(1+u*conj(u))^2*conj(d2w_uconjMu)...
    -24*i*u/(1+u*conj(u))^3*conj(dw_conjMu);
d2rho_conjuconjMuSub = conj(d2rho_uMuSub);

d2rho_uuSub= i*d4w_uuuu -6*i*conj(u)/(1+u*conj(u))*d3w_uuu...
    +24*i*conj(u)^2/(1+u*conj(u))^2*d2w_uu -60*i*conj(u)^3/(1+u*conj(u))^3*dw_u...
    +72*i*conj(u)^4/(1+u*conj(u))^4*w...
    -12*i*conj(u)/(1+u*conj(u))^3*conj(dw_u)...
    +72*i*u*conj(u)/(1+u*conj(u))^4*conj(w)...
    -24*i/(1+u*conj(u))^3*conj(w);
d2rho_conjuconjuSub = conj(d2rho_uuSub);

rhoSet = [rho, drho_Mu, drho_conjMu, drho_vnormv, drho_u, drho_conju,...
    d2rho_uMu, d2rho_conjuconjMu, d2rho_uu, d2rho_conjuconju];
rhoSubSet = [rhoSub, drho_MuSub, drho_conjMuSub, drho_vnormvSub, drho_uSub,...
    drho_conjuSub, d2rho_uMuSub, d2rho_conjuconjMuSub, d2rho_uuSub, d2rho_conjuconjuSub];
for m=1:6
    for n=1:6
        temp=RicFlat(m,n);
        temp=subs(temp, rhoSet, rhoSubSet);
        RicFlat(m,n)=temp;
    end
end
MVarRicF =[d2w_uMu, d2w_uconjMu, d2w_uu, d2w_uvnormv, d3w_uuMu,...
    d3w_uuconjMu, d3w_uuu, d3w_uuvnormv, d4w_uuuMu, d4w_uuuu,...
    dw_Mu, dw_conjMu, dw_u, dw_vnormv, u, w];
for m=1:6
    for n=1:6
        temp=RicFlat(m,n);
        RicFlat(m,n) = complex_simple3(temp, MVarRicF);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSet2 = [T4, aM, aV, bV, d2h11_uMu, d2h11_uconju,...
    d2h11_uu, d2w_uconjMu, d2w_uu, dT4_conjMu,...
    dT4_conju, dT4_u, daM_Mu, daM_u, dh11_Mu,...
    dh11_conju, dh11_u, dw_conjMu, dw_u, h11, u, w];
subSet2 = [0, 0, 0, 0, 0, 0,...
    0, d2w_uconjMu, d2w_uu, 0,...
    0, 0, 0, 0, 0,...
    0, 0, dw_conjMu, dw_u, 0, u, w];
MVarSet2 = [u, w, dw_u, d2w_uu, dw_conjMu, d2w_uconjMu];
for m=1:2
    for n=1:2
        temp = ricci(m,n);
        temp = subs(temp, variableSet2, subSet2);
        ricci(m,n) = complex_simple3(temp,MVarSet2);
    end
end
variableSet3 = [u, w, dw_u, dw_conjMu, ...
    T4, aM, aV, bV, h11, dT4_u, dT4_conju, dT4_conjMu,...
    dh11_u, dh11_conju, dh11_conjMu];
subSet3 = [u, w, dw_u, dw_conjMu, ...
    0, 0, 0, 0, 0, 0, 0, 0,...
    0, 0, 0];
for n=1:2
   for k=1:2
       tempT = Gamma.T(n,k);
       tempT = subs(tempT, variableSet3, subSet3);
       for m=1:2
           temp1 = Gamma.holo(m,n,k);
           temp1 = subs(temp1, variableSet3, subSet3);
           temp2 = Gamma.antiholo(m,n,k);
           temp2 = subs(temp2, variableSet3, subSet3);
           Gamma.holo(m,n,k) = complex_simple3(temp1, MVarSet2);
           Gamma.antiholo(m,n,k) = complex_simple3(temp2, MVarSet2);
       end
       Gamma.T(n,k) = complex_simple3(tempT, MVarSet2);
   end
end
clearvars m n k temp1 temp2 tempT temp
clearvars variableSet2 variableSet3 subSet2 subSet3
save('Data_Ricci_Check_Jan30.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RicZ11 = Ric(Z1,conj(Z1));
load('Data_Ricci_Check_Jan30.mat');
hTwo = [0, -i;
    i, 0];
term1 = i/2*(Gamma.holo(1,1,1)+Gamma.holo(1,2,2)...
    -conj(Gamma.antiholo(1,1,1)) -conj(Gamma.antiholo(1,2,2)));
term1bar = i/2*(Gamma.antiholo(1,1,1)+Gamma.antiholo(1,2,2)...
    -conj(Gamma.holo(1,1,1)) -conj(Gamma.holo(1,2,2)));
term2 = i/2*(Gamma.holo(2,1,1)+Gamma.holo(2,2,2)...
    -conj(Gamma.antiholo(2,1,1)) -conj(Gamma.antiholo(2,2,2)));
term2bar = i/2*(Gamma.antiholo(2,1,1)+Gamma.antiholo(2,2,2)...
    -conj(Gamma.holo(2,1,1)) -conj(Gamma.holo(2,2,2)));

RicZ11 = RicFlat(1,2)-term1*RicFlat(6,2)-term1bar*RicFlat(1,6)...
    +term1*term1bar*RicFlat(6,6);
RicZ12 = RicFlat(1,4)-term1*RicFlat(6,4)-term2bar*RicFlat(1,6)...
    +term1*term2bar*RicFlat(6,6);
RicZ21 = RicFlat(3,2)-term2*RicFlat(6,2)-term1bar*RicFlat(3,6)...
    +term2*term1bar*RicFlat(6,6);
RicZ22 = RicFlat(3,4)-term2*RicFlat(6,4)-term2bar*RicFlat(3,6)...
    +term2*term2bar*RicFlat(6,6);

test11 = 1/2*ricci(1,1)+1/12*rhoSub*hTwo(1,1);
test12 = 1/2*ricci(1,2)+1/12*rhoSub*hTwo(1,2);
test21 = 1/2*ricci(2,1)+1/12*rhoSub*hTwo(2,1);
test22 = 1/2*ricci(2,2)+1/12*rhoSub*hTwo(2,2);
diff11 = complex_simple3(RicZ11-test11, MVarSet2);
diff12 = complex_simple3(RicZ12-test12, MVarSet2);
diff21 = complex_simple3(RicZ21-test21, MVarSet2);
diff22 = complex_simple3(RicZ22-test22, MVarSet2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%