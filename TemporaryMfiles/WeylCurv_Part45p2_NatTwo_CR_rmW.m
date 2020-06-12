% Combine and organize 
% (1) WeylCurv_temp1Part45_NatTwo_CR_rmW.m
% (2) WeylCurv_temp2Part45_NatTwo_CR_rmW.m
% 
% b/f WeylCurv_Part45_NatTwo_CR_rmW.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataTemp_WeylCurv_Part45_Aug31.mat');
assumeAlso([x y z u1 u2 gamma K0 dK0_y d2K0_yy], 'real');
MVarPsi3 = [u, p, Psi, dPsi_p, d2Psi_pp, d3Psi_ppp,...
    Zeta1, dZeta1_p, d2Zeta1_pp, d3Zeta1_ppp,...
    Zeta2, dZeta2_p, d2Zeta2_pp, d3Zeta2_ppp];
syms conju
WN1212 = WeylSN(1,5)/(1+u*conj(u))^2;
WN1215 = WeylSN(2,5)/(1+u*conj(u));
WN1225 = WeylSN(3,5)/(1+u*conj(u));
WN1515 = WeylSN(4,5);
WN1525 = WeylSN(5,5);
WN2525 = WeylSN(6,5);
WNvector = [WN1212, WN1215, WN1225, WN1515, WN1525, WN2525];

PsiN = [Psi, dPsi_p, d2Psi_pp, d3Psi_ppp];
Zeta1N = [Zeta1, dZeta1_p, d2Zeta1_pp, d3Zeta1_ppp];
Zeta2N = [Zeta2, dZeta2_p, d2Zeta2_pp, d3Zeta2_ppp];
K0N = [K0, dK0_y, d2K0_yy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Condition (1): dZeta1_p = 0.
for j=1:6
    temp = WNvector(j);
    temp = subs(temp, Zeta1N, [Zeta1, 0, 0, 0]);
    temp = complex_simple3(temp, MVarPsi3);
    temp = subs(temp, conj(u), conju);
    WNvector(j) = temp;
    clear temp
end

indexNonLinear = {'u1212', 'u1215', 'u1225', 'u1515', 'u1525', 'u2525';
    'term1212', 'term1215', 'term1225', 'term1515', 'term1525',...
    'term2525'};

[CoeffSN.term1212, CoeffSN.u1212] = coeffs(WNvector(1), [u,conju]);
[CoeffSN.term1215, CoeffSN.u1215] = coeffs(WNvector(2), [u,conju]);
[CoeffSN.term1225, CoeffSN.u1225] = coeffs(WNvector(3), [u,conju]);
[CoeffSN.term1515, CoeffSN.u1515] = coeffs(WNvector(4), [u,conju]);
[CoeffSN.term1525, CoeffSN.u1525] = coeffs(WNvector(5), [u,conju]);
[CoeffSN.term2525, CoeffSN.u2525] = coeffs(WNvector(6), [u,conju]);

latex_Part45p2 = fopen('latex__Part45p2_c1.txt','w');
fprintf(latex_Part45p2, 'The term W1212 when dZeta1_p = 0\n');
N = length(CoeffSN.u1212);
for j=1:N
    uTemp = CoeffSN.u1212(j);
    termTemp = CoeffSN.term1212(j);
    fprintf(latex_Part45p2,'$%s $: $%s $ \n',latex(uTemp),latex(termTemp));
    fprintf(latex_Part45p2, '%s\n', ' ');
    clearvars uTemp termTemp
end
fclose(latex_Part45p2);


