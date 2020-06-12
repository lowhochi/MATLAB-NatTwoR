% ScalarCurv_NatTwo_CR_rmW.m
% load('DataMain4_Nat_CR_rmW_Ap21.mat');
load('DataMain4S_Nat_CR_rmW_Apr21.mat'); %with STwo defined as below
assumeAlso([x y z u1 u2 gamma], 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STwo = complex_simple3(S, symSetRCurv0);
% % Verify that S = 5/3*rho.
% save('DataMain4S_Nat_CR_rmW_Apr21.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(1)
symSetS = symvar(STwo);
phiW = d2w_uu - 6*conj(u)*dw_u/(1+u*conj(u)) + 12*conj(u)^2*w/(1+u*conj(u))^2;
S_h11Term = 2*d2h11_uconju - 4*u/(1+u*conj(u))*dh11_u...
    -4*conj(u)/(1+u*conj(u))*dh11_conju + 8*h11/(1+u*conj(u)) -14*h11/(1+u*conj(u))^2;
S_T4Term = -3*i*dT4_u + 5*i*conj(dT4_u) + 6*i*conj(u)/(1+u*conj(u))*T4...
    - 10*i*u/(1+u*conj(u))*conj(T4);
S_aMuaVTerm = 2*i*aV - 3*i/2*daMu_conju + 3*i/2*conj(daMu_conju)...
    + 3*i*u/(1+u*conj(u))*aMu - 3*i*conj(u)/(1+u*conj(u))*conj(aMu);

Stest = 1/6*rho + 3*i/2*(phiW - conj(phiW)) + S_h11Term + S_T4Term + S_aMuaVTerm;
Sdifference = complex_simple3(STwo-Stest, symSetS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(2)
syms conjaMu conjdaMu_conju
SThree = subs(Stest, conj(aMu), conjaMu);
SThree = subs(SThree, conj(daMu_conju), conjdaMu_conju);
%
aMu_Sub = 2*i*conj(u)/(1+u*conj(u))*h11 - i*dh11_u -2*conj(T4);
conjaMu_Sub = -2*i*u/(1+u*conj(u))*h11 + i*dh11_conju - 2*T4;
aV_Sub = -conj(aV) - dT4_u - conj(dT4_u) + 2*conj(u)/(1+u*conj(u))*T4...
    + 2*u/(1+u*conj(u))*conj(T4);
daMu_conju_Sub = 2*i*(1+2*u*conj(u))/(1+u*conj(u))^2*h11...
    - 2*i*u/(1+u*conj(u))*dh11_u + 4*conj(u)/(1+u*conj(u))*T4...
    - 2*conj(aV) - 2*dT4_u - 2*conj(dT4_u);
conjdaMu_conju_Sub = -2*i*(1+2*u*conj(u))/(1+u*conj(u))^2*h11...
    + 2*i*conj(u)/(1+u*conj(u))*dh11_conju - 4*conj(u)/(1+u*conj(u))*T4 + 2*conj(aV);
%
Stest3 = subs(SThree, aMu, aMu_Sub);
Stest3 = subs(Stest3, conjaMu, conjaMu_Sub);
Stest3 = subs(Stest3, aV, aV_Sub);
Stest3 = subs(Stest3, daMu_conju, daMu_conju_Sub);
Stest3 = subs(Stest3, conjdaMu_conju, conjdaMu_conju_Sub);
symSetS3 = symvar(Stest3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%