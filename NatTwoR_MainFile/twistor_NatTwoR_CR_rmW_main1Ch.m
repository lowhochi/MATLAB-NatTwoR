% twistor_NatTwoR_CR_rmW_main1Ch.m
% Ch = Chern curvature tensor C(m,n,k,ll);
load('DataMain15_NatTwoR_CR_rmW.mat');

countIndexChern = [1,1,1,1; 1,1,1,2; 1,1,2,1; 1,1,2,2;
    1,2,1,1; 1,2,1,2; 1,2,2,1; 1,2,2,2;
    2,1,1,1; 2,1,1,2; 2,1,2,1; 2,1,2,2;
    2,2,1,1; 2,2,1,2; 2,2,2,1; 2,2,2,2];

MVarMain1Ch = [];
for j=1:16
    m = countIndexChern(j,1);
    n = countIndexChern(j,2);
    k = countIndexChern(j,3);
    ll= countIndexChern(j,4);
    tempSet = symvar(R(m,n,k,ll));
    MVarMain1Ch = union(MVarMain1Ch, tempSet);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find out the chern tensor coefficients.
% Chern(m,n,k,l) = C_{m\bar{n}k\bar{l}}.
% Main15: R(m,n,k,l) = R_m^n_k_ol where R(X_k, conj_X_l)X_m = R_m^n_k_ol X_n. 
% Main1Ch:RThree(m,n,k,ll) = R(m,p,k,ll)*h_{p\bar n} = R_{m\bar nk\bar ll}

Chern = sym('Chern',[2,2,2,2]);
RThree = sym('RThree',[2 2 2 2]);
symSetChern = []; % specify the symbolic variables in Chern tensor

for j=1:16
    m = countIndexChern(j,1);
    n = countIndexChern(j,2);
    k = countIndexChern(j,3);
    ll= countIndexChern(j,4);
    
   temp = 0;
   for p=1:2
       temp = temp+ h(p,n)*R(m,p,k,ll);
   end
   RThree(m,n,k,ll) = temp;
   tempPart1 = RThree(m,n,k,ll);
   tempPart2 = (-1/4)*(h(k,ll)*ric(m,n) + h(m,ll)*ric(k,n)...
                    + ric(k,ll)*h(m,n) + ric(m,ll)*h(k,n));
   tempPart3 = (rho/12)*(h(k,ll)*h(m,n) + h(m,ll)*h(k,n));
   tempChern = tempPart1 + tempPart2 + tempPart3;                
   Chern(m,n,k,ll) = complex_simple3(tempChern,MVarMain1Ch);
   tempSet = symvar(Chern(m,n,k,ll));
   symSetChern = union(symSetChern, tempSet);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars j m n k ll p
clearvars temp tempChern tempPart1 tempPart2 tempPart3 tempSet
list_of_variables_Main1Ch = who;
save('DataMain1Ch_NatTwoR_CR_rmW.mat');


