load('DataTemp2_WeylCurv_Part45_Aug31.mat');
assumeAlso([x y z u1 u2 gamma K0 dK0_y d2K0_yy], 'real');
syms Psi0 b0 b1
syms r0 real
% dPsi_p = 0 and d2Psi_pp = 0.
Zeta1Sub = [(-i*r0/conj(Psi0))*p+b0, -i*r0/conj(Psi0), 0];
Zeta2Sub = [Zeta2, -i*r0*Psi0/conj(Psi0)^2*p+b1, -i*r0*Psi0/conj(Psi0)^2];
b1Sub = 2*Psi0*conj(Psi0)*(b0-conj(b0)-i*K0)/(conj(Psi0)^2-Psi0^2);
PsiSub = [Psi0, 0, 0];
MVarPsi3 = union(MVarPsi3,[Psi0, b0, b1]);
% WNvector=[WN1212, WN1215, WN1225, WN1515, WN1525, WN2525];
WNvectorTwo = sym('WNvectorTwo',[1,6]);
for kk=1:6
    temp = WNvector(kk);
    temp = subs(temp, [Psi,dPsi_p,d2Psi_pp], PsiSub);
    temp = subs(temp, [dK0_y,d2K0_yy], [0,0]);
    temp = subs(temp,[Zeta1,dZeta1_p,d2Zeta1_pp],Zeta1Sub);
    temp = subs(temp,[Zeta2,dZeta2_p,d2Zeta2_pp],Zeta2Sub);
    temp = subs(temp, b1, b1Sub);
    temp = complex_simple3(temp, MVarPsi3);
    WNvectorTwo(kk) = temp;
    clear temp
end

%%
load('DataTemp2_WeylCurv_Part45_Aug31.mat');