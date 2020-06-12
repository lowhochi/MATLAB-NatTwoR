% WeylCurv_Part42_NatTwo_CR_rmW.m
load('Data_WeylChern_Part4_May30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms C0
WeylSxz = sym('WeylSxz',[120 5]);
subSetSxz = sym('subSetSxz', [length(symSetWeylS),1]);
for j=1:length(symSetWeylS)
    variable = string(symSetWeylS(j));
    switch variable
        case string(u)
            subSetSxz(j) = u;
        case string(lambda0)
            subSetSxz(j) = -y/2;
        case string(lambda1)
            subSetSxz(j) = x+i*z;
        case string(dlambda0_Mu)
            subSetSxz(j) = -u;
        case string(dlambda0_conjMu)
            subSetSxz(j) = -conj(u);
        case string(dlambda0_vnormv)
            subSetSxz(j) = (u*conj(u)-1)/(2*(1+u*conj(u)));
        case string(dlambda1_Mu)
            subSetSxz(j) = -2;
        case string(dlambda1_conjMu)
            subSetSxz(j) = 2*conj(u)^2;
        case string(dlambda1_vnormv)
            subSetSxz(j) = T1 + i*T3;
        otherwise 
            subSetSxz(j) = 0;
    end
end

for j=1:120
    m = countSetTwo(j,1);
    n = countSetTwo(j,2);
    k = countSetTwo(j,3);
    ll = countSetTwo(j,4);
    temp = WeylSS(j,5);
    temp = subs(temp, symSetWeylS, subSetSxz);
    temp = subs(temp, conj(K), K);
    temp = complex_simple3(temp, [u]);
    WeylSxz(j,:) = [m,n,k,ll,temp]; % WeylSxz(j,5)=0 for any j.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%