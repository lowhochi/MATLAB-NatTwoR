% Condition dPsi_p = 0 and d2Psi_pp = 0.
load('DataTemp_WeylCurv_Part45_Aug31.mat');
assumeAlso([x y z u1 u2 gamma K0 dK0_y d2K0_yy], 'real');
syms conju
MVarPsi3 = [u, p, Psi, dPsi_p, d2Psi_pp, d3Psi_ppp,...
    Zeta1, dZeta1_p, d2Zeta1_pp, d3Zeta1_ppp,...
    Zeta2, dZeta2_p, d2Zeta2_pp, d3Zeta2_ppp];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WN1212 = WeylSN(1,5)/(1+u*conj(u))^2;
WN1215 = WeylSN(2,5)/(1+u*conj(u));
WN1225 = WeylSN(3,5)/(1+u*conj(u));
WN1515 = WeylSN(4,5);
WN1525 = WeylSN(5,5);
WN2525 = WeylSN(6,5);
WNvector = [WN1212, WN1215, WN1225, WN1515, WN1525, WN2525];
for j=1:6
    temp = WNvector(j);
    temp = complex_simple3(temp, MVarPsi3);
    temp = subs(temp, conj(u), conju);
    WNvector(j) = temp;
    clear temp
end
indexNonLinear = {'u1212', 'u1215', 'u1225', 'u1515', 'u1525', 'u2525';
    'term1212', 'term1215', 'term1225', 'term1515', 'term1525',...
    'term2525'};

[CoeffSN.term1212, CoeffSN.u1212] = coeffs(WNvector(1), [u,conju,y]);
[CoeffSN.term1215, CoeffSN.u1215] = coeffs(WNvector(2), [u,conju,y]);
[CoeffSN.term1225, CoeffSN.u1225] = coeffs(WNvector(3), [u,conju,y]);
[CoeffSN.term1515, CoeffSN.u1515] = coeffs(WNvector(4), [u,conju,y]);
[CoeffSN.term1525, CoeffSN.u1525] = coeffs(WNvector(5), [u,conju,y]);
[CoeffSN.term2525, CoeffSN.u2525] = coeffs(WNvector(6), [u,conju,y]);

save('DataTemp2_WeylCurv_Part45_Aug31.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
load('DataTemp2_WeylCurv_Part45_Aug31.mat');
assumeAlso([x y z u1 u2 gamma K0 dK0_y d2K0_yy], 'real');
syms Psi0 b0 b1
syms r0 real
% dPsi_p = 0 and d2Psi_pp = 0.
% WN1212:
term1212v0 = sym('term1212v0',[1,length(CoeffSN.u1212)]);
Zeta1Sub = [(-i*r0/conj(Psi0))*p+b0, -i*r0/conj(Psi0), 0];
Zeta2Sub = [Zeta2, -i*r0*Psi0/conj(Psi0)^2*p+b1, -i*r0*Psi0/conj(Psi0)^2];
b1Sub = 2*Psi0*conj(Psi0)*(b0-conj(b0)-i*K0)/(conj(Psi0)^2-Psi0^2);
PsiSub = [Psi0, 0, 0];
MVarPsi3 = union(MVarPsi3,[Psi0, b0, b1]);
for jj=1:length(CoeffSN.u1212)
    termTemp = CoeffSN.term1212(jj);
    termTemp = subs(termTemp,[Psi,dPsi_p,d2Psi_pp],PsiSub);
    termTemp = subs(termTemp, [dK0_y, d2K0_yy], [0, 0]);
    termTemp = subs(termTemp,[Zeta1,dZeta1_p,d2Zeta1_pp],Zeta1Sub);
    termTemp = subs(termTemp,[Zeta2,dZeta2_p,d2Zeta2_pp],Zeta2Sub);
%   termTemp = subs(termTemp, b1, b1Sub);
    termTemp = complex_simple3(termTemp, MVarPsi3);
    term1212v0(jj) = termTemp;
    clear termTemp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%