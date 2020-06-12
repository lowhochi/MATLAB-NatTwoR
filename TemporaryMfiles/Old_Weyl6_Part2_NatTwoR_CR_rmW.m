% Weyl6_Part2_NatTwoR_CR_rmW.m
load('DataWeyl6_NatTwoR_CR_rmW.mat');
assumeAlso(gamma,'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume rho=0 when M is a flat space.
% w(x,u) = lamba0 +lambda1*u +K*u^2 -conj(lambda1)*u^3
%   + conj(lambda0)*u^4;
Y0 = 1+u*conj(u);
% Use Aterm, Bterm, Cterm to express all of them.
W1214 = WeylFlat(3,5);
W1223 = WeylFlat(6,5);
W1235 = WeylFlat(11,5);
W1245 = WeylFlat(13,5);
W1415 = WeylFlat(31,5);
W1425 = WeylFlat(35,5);
W1523 = WeylFlat(45,5);
W1535 = WeylFlat(50,5);
W1545 = WeylFlat(52,5);
W2325 = WeylFlat(68,5);
W2535 = WeylFlat(88,5);
W2545 = WeylFlat(90,5);

wLinearVec = [W1214, W1223, W1235, W1245, W1415, W1425, W1523, W1535,...
    W1545, W2325, W2535, W2545];

Aterm = dlambda1_z + conj(dlambda1_z) -i*dlambda1_x +i*conj(dlambda1_x);
Bterm = dK_z +i*dK_x +2*i*dlambda0_x -2*dlambda0_z +i*dlambda1_y;
Cterm = dlambda1_z +i*dlambda1_x +4*i*dlambda0_y;
% Convention difference: with a minus sign on every term.
test1214 = Aterm/2*(1 -4*u*conj(u) +u^2*conj(u)^2)...
    + (Bterm*conj(u)+conj(Bterm)*u)*(1-u*conj(u))...
    - (Cterm*conj(u)^2 +conj(Cterm)*u^2);
% diff1214 = W1214 + test1214;
% diff1214 = complex_simple3(diff1214, MVarW6);

test1235 = 1/Y0*(conj(u)^2/4*(u*conj(u)-3)*Bterm...
    +1/4*(1-3*u*conj(u))*conj(Bterm) +3/4*conj(u)*(u*conj(u)-1)*Aterm...
    +conj(u)^3/2*Cterm -u/2*conj(Cterm));
diff1235 = W1235 + test1235;
diff1235 = subs(diff1235, conjKset, Kset);
diff1235 = complex_simple3(diff1235, MVarW6);

test1415 = 1/Y0*(3*conj(u)/4*(u*conj(u)-1)*Aterm...
    +conj(u)^2/4*(u*conj(u)-3)*Bterm + 1/4*(1-3*u*conj(u))*conj(Bterm)...
    + 1/2*(conj(u)^3*Cterm-u*conj(Cterm)));
diff1415 = W1415 + test1415;
diff1415 = subs(diff1415, conjKset, Kset);
diff1415 = complex_simple3(diff1415, MVarW6);

test1425 = -conj(test1235);
diff1425 = W1425 + test1425;
diff1425 = subs(diff1425, conjKset, Kset);
diff1425 = complex_simple3(diff1425, MVarW6);

test1535 = 1/Y0^2*(3*conj(u)^2/4*Aterm +conj(u)^3/2*Bterm...
    -conj(u)/2*conj(Bterm) -conj(u)^4/4*Cterm -conj(Cterm)/4);
diff1535 = W1535 + test1535;
diff1535 = subs(diff1535, conjKset, Kset);
diff1535 = complex_simple3(diff1535, MVarW6);

test1545 = 1/(4*Y0^2)*test1214;
diff1545 = W1545 + test1545;
diff1545 = subs(diff1545, conjKset, Kset);
diff1545 = complex_simple3(diff1545, MVarW6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%