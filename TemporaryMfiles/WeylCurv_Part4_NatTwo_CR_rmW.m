load('Data_WeylChern_Part3_May30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MVarSTwo = [u, lambda0, lambda1, K, ...
    dlambda0_Mu, dlambda0_conjMu, dlambda0_vnormv, ...
    d2lambda0_MuMu, d2lambda0_MuconjMu, d2lambda0_Muvnormv, ...
    d2lambda0_conjMuconjMu, d2lambda0_conjMuvnormv, d2lambda0_vnormvvnormv, ...
    dlambda1_Mu, dlambda1_conjMu, dlambda1_vnormv, ...
    d2lambda1_MuMu, d2lambda1_MuconjMu, d2lambda1_Muvnormv, ...
    d2lambda1_conjMuconjMu, d2lambda1_conjMuvnormv, d2lambda1_vnormvvnormv, ...
    dK_Mu, dK_conjMu, dK_vnormv, ...
    d2K_MuMu, d2K_MuconjMu, d2K_Muvnormv, ...
    d2K_conjMuconjMu, d2K_conjMuvnormv, d2K_vnormvvnormv]; 
symSetWeylS = union(MVarSTwo, [aMu;h11]); % symvar of WeylS
countSet225 = [];
for m=1:6
    for n=1:6
        for k=1:6
            for ll=1:6
                rowTemp = [m,n,k,ll];
                if (m<n)&&(k<ll)
                    countSet225 = [countSet225; rowTemp];
                end
            end
        end
    end
end
numberCount = zeros(225,2);
countSetTwo = []; % length(countSetTwo)=120
for j=1:225
    numberCount(j,1) = countSet225(j,1)*10 + countSet225(j,2);
    numberCount(j,2) = countSet225(j,3)*10 + countSet225(j,4);
    if numberCount(j,1)<=numberCount(j,2)
        countSetTwo = [countSetTwo; countSet225(j,:)];
    end      
end

WeylSS = sym('WeylSS',[120,5]);
for j=1:120
    m = countSetTwo(j,1);
    n = countSetTwo(j,2);
    k = countSetTwo(j,3);
    ll = countSetTwo(j,4);
    temp = WeylS(m,n,k,ll);
    WeylSS(j,:)=[m,n,k,ll, temp];
end
clearvars j numberCount countSet225 m n k ll rowTemp temp
save('Data_WeylChern_Part4_May30.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section2/ w(x,u)= K(y)*u^2
load('Data_WeylChern_Part4_May30.mat');
syms K2 K22 real % K2 = dK/dy, K22 = d2K/dy2
WeylSK = sym('WeylSK',[120 5]);
subSetSK = sym('subSetSK', [length(symSetWeylS),1]);

for j=1:length(symSetWeylS)
    variable = string(symSetWeylS(j));
    switch variable
        case string(u)
            subSetSK(j) = u;
        case string(K)
            subSetSK(j) = K;
        case string(dK_Mu)
            subSetSK(j) = 2*u*K2;
        case string(dK_conjMu)
            subSetSK(j) = 2*conj(u)*K2;
        case string(dK_vnormv)
            subSetSK(j) = (1-u*conj(u))/(1+u*conj(u))*K2;
        case string(d2K_MuMu)
            subSetSK(j) = 4*u^2*K22;
        case string(d2K_conjMuconjMu)
            subSetSK(j) = 4*conj(u)^2*K22;
        case string(d2K_MuconjMu)
            subSetSK(j) = 4*u*conj(u)*K22;
        case string(d2K_Muvnormv)
            subSetSK(j) = 2*u*(1-u*conj(u))/(1+u*conj(u))*K22;
        case string(d2K_conjMuvnormv)
            subSetSK(j) = 2*conj(u)*(1-u*conj(u))/(1+u*conj(u))*K22;
        case string(d2K_vnormvvnormv)
            subSetSK(j) = (1-u*conj(u))^2/(1+u*conj(u))^2*K22;
        otherwise 
            subSetSK(j) = 0;
    end
end

for j=1:120
    m = countSetTwo(j,1);
    n = countSetTwo(j,2);
    k = countSetTwo(j,3);
    ll = countSetTwo(j,4);
    temp = WeylSS(j,5);
    temp = subs(temp, symSetWeylS, subSetSK);
    temp = subs(temp, conj(K), K);
    temp = complex_simple3(temp, [u,K]);
    WeylSK(j,:) = [m,n,k,ll,temp]; % WeylSK(k,5)=0 for any j.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%