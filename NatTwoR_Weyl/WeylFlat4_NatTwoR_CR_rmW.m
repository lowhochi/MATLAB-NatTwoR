% WeylFlat4_NatTwoR_CR_rmW.m
% Verify that system(F) in (7.16) is correct.
load('DataWeylFlat3_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableSet = [d2phi_pconjp, d2phi_pp, d2phi_py, d2phi_yy,...
    d3phi_pconjpy, d3phi_ppconjp, d3phi_ppp, d3phi_ppy,...
    d3phi_pyy, d3phi_yyy];
syms p
syms Psi dPsi_p d2Psi_pp
syms Zeta1 dZeta1_p d2Zeta1_pp
syms Zeta2 dZeta2_p d2Zeta2_pp
syms K0 dK0_y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2phi_ppSub = y*dPsi_p+ i*dZeta1_p*conj(p) + i*Zeta2;
d2phi_pconjpSub = i*(Zeta1-conj(Zeta1));
d2phi_pySub = Psi;
d2phi_yySub = -2*K0;
d3phi_pconjpySub = 0;
d3phi_ppconjpSub = i*dZeta1_p;
d3phi_pppSub = y*d2Psi_pp +i*d2Zeta1_pp*conj(p) +i*dZeta2_p;
d3phi_ppySub = dPsi_p;
d3phi_yyySub = -2*dK0_y;
d3phi_pyySub = 0;

MVarFlat3 = [u, y, p, Psi, dPsi_p, d2Psi_pp, K0, dK0_y,...
    Zeta1, dZeta1_p, d2Zeta1_pp, Zeta2, dZeta2_p, d2Zeta2_pp];
subSet = [d2phi_pconjpSub, d2phi_ppSub, d2phi_pySub, d2phi_yySub,...
    d3phi_pconjpySub, d3phi_ppconjpSub, d3phi_pppSub, d3phi_ppySub,...
    d3phi_pyySub, d3phi_yyySub];

WflatSet03 = sym('Wflat',[1,4]);
for number=1:4
    temp = WflatSet02(number);
    temp = subs(temp, variableSet, subSet);
    temp = subs(temp, conj(y), y);
    temp = subs(temp, conj(K0), K0);
    temp = subs(temp, conj(dK0_y), dK0_y);
    WflatSet03(number) = complex_simple3(temp, MVarFlat3);
end
clearvars variableSet subSet temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
termF1 = Psi*conj(dZeta1_p) +dZeta1_p*conj(Psi)...
    +1/2*dPsi_p*(conj(dZeta1_p)*p+conj(Zeta2))...
    +1/2*conj(dPsi_p)*(dZeta1_p*conj(p)+Zeta2);
% % % % %
termF21 = i*dPsi_p*conj(d2Psi_pp);
termF22 = dPsi_p*(conj(d2Zeta1_pp)*p+conj(dZeta2_p))...
    -conj(d2Psi_pp)*(dZeta1_p*conj(p)+Zeta2)... 
    -dZeta1_p*conj(dPsi_p);
termF23 = -2*conj(dZeta1_p)*(K0+i*(Zeta1-conj(Zeta1)))...
    +i*(conj(d2Zeta1_pp)*p+conj(dZeta2_p))*(dZeta1_p*conj(p)+Zeta2)...
    +i*dZeta1_p*(conj(dZeta1_p)*p+conj(Zeta2));
termF2 = termF21*y^2 + termF22*y + termF23;
% % % % %
termF31 = i*(Psi*conj(d2Psi_pp) -dK0_y*conj(dPsi_p));
termF32 = -dK0_y*(conj(dZeta1_p)*p+conj(Zeta2))...
    -i*conj(dPsi_p)*(K0+i*(Zeta1-conj(Zeta1)))...
    +Psi*(conj(d2Zeta1_pp)*p+conj(dZeta2_p))...
    -conj(Psi)*conj(dZeta1_p);
termF3 = termF31*y +termF32;
% % % % %
Wflat1212Test = (1+u*conj(u))^2*...
    (-termF1*(u^2*conj(u)^2-4*u*conj(u)+1)...
    -termF2*(u^2*conj(u)-u)-conj(termF2)*(conj(u)^2*u-conj(u))...
    -termF3*u^2 - conj(termF3)*conj(u)^2);
% % % % %
Wflat1215Test = (1+u*conj(u))*...
    (-3/2*termF1*(u*conj(u)^2-conj(u))...
    -1/4*termF2*(3*u*conj(u)-1)-1/4*conj(termF2)*(3*conj(u)^2-u*conj(u)^3)...
    -1/2*termF3*u +1/2*conj(termF3)*conj(u)^3);
% % % % %
Wflat1515Test = -3/2*termF1*conj(u)^2 ...
    -1/2*termF2*conj(u) +1/2*conj(termF2)*conj(u)^3 ...
    -1/4*termF3 -1/4*conj(termF3)*conj(u)^4;
% % % % %
Wflat1525Test = -1/4*termF1*(u^2*conj(u)^2-4*u*conj(u)+1)...
    -1/4*termF2*(u^2*conj(u)-u)...
    -1/4*conj(termF2)*(u*conj(u)^2-conj(u))...
    -1/4*termF3*u^2 -1/4*conj(termF3)*conj(u)^2;
% % % % %
WflatTestSet = [Wflat1212Test, Wflat1215Test, Wflat1515Test, Wflat1525Test];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for number=1:4
    myDiff = WflatSet03(number) -WflatTestSet(number);
    myDiff = subs(myDiff, conj(y), y);
    myDiff = subs(myDiff, conj(K0), K0);
    myDiff = subs(myDiff, conj(dK0_y), dK0_y);
    myDiff = complex_simple3(myDiff, MVarFlat3);
    disp('The difference is');
    disp(myDiff);
end
clearvars number myDiff temp
save('DataWeylFlat4_NatTwoR_CR_rmW.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%