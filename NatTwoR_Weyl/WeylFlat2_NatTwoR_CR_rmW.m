% WeylFlat2_NatTwoR_CR_rmW.m
load('DataWeylFlat1_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WflatSet01 = WflatSet(1:4);
Y0 = 1+u*conj(u);
Aterm = dlambda1_z + conj(dlambda1_z) -i*dlambda1_x +i*conj(dlambda1_x);
Bterm = dK_z +i*dK_x +2*i*dlambda0_x -2*dlambda0_z +i*dlambda1_y;
Cterm = dlambda1_z +i*dlambda1_x +4*i*dlambda0_y;
test1214 = Aterm/2*(1 -4*u*conj(u) +u^2*conj(u)^2)...
    + (Bterm*conj(u)+conj(Bterm)*u)*(1-u*conj(u))...
    - (Cterm*conj(u)^2 +conj(Cterm)*u^2);
test1235 = 1/Y0*(conj(u)^2/4*(u*conj(u)-3)*Bterm...
    +1/4*(1-3*u*conj(u))*conj(Bterm) +3/4*conj(u)*(u*conj(u)-1)*Aterm...
    +conj(u)^3/2*Cterm -u/2*conj(Cterm));
test1535 = 1/Y0^2*(3*conj(u)^2/4*Aterm +conj(u)^3/2*Bterm...
    -conj(u)/2*conj(Bterm) -conj(u)^4/4*Cterm -conj(Cterm)/4);
test1545 = 1/(4*Y0^2)*test1214;
testSet = [test1214, test1235, test1535, test1545];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diffSet = sym('diffSet',[1,4]);
for number=1:4 
    temp = WflatSet01(number)+testSet(number);
    temp = subs(temp, conjKset, Kset);
    temp = complex_simple3(temp, MVarFlat1);
    diffSet(number) = temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%