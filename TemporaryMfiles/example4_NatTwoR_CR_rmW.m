% example4_NatTwoR_CR_rmW.m
variable_example4_CR_rmW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metric tensor g = d2phi
% basis = {d/dp, d/dconjp, d/dy}
d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, d2phi_conjpconjp, d2phi_conjpy;
    d2phi_py, d2phi_conjpy, d2phi_yy];
d2phiInv = sym('d2phiInv',[3 3]);
for ii=1:3
    for j=1:3
        temp = myCofactor(d2phi,ii,j);
        d2phiInv(ii,j) = temp/detPhi;
    end
end
Gamma = sym('Gamma',[3,3,3]);
for ii=1:3
    for j=1:3
        for k=1:3
            temp=0;
            dgij = df_example_pconjpy(d2phi(ii,j),CVar,derivativeDict); 
            for ll=1:3
                dgil = df_example_pconjpy(d2phi(ii,ll),CVar,derivativeDict); 
                dgjl = df_example_pconjpy(d2phi(j,ll),CVar,derivativeDict); 
                temp = temp + 1/2*d2phiInv(k,ll)*(dgil(j)...
                    + dgjl(ii) -dgij(ll));
            end
            temp = myRealVariableFun(temp, realVariable);
            Gamma(ii,j,k) = complex_simple3(temp, MVar);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars ii j k ll temp dgij dgil dgjl
% load('DataTemp_example4_v0.mat');
% (1,3) Riemann tensor
% Riem(i,j,k,l) = R_{ijk}^l
Riem = sym('Riem',[3 3 3 3]);
for ii=1:3
    for j=1:3
        for k=1:3
            for ll=1:3
                dGammajkl = df_example_pconjpy(Gamma(j,k,ll),...
                    CVar, derivativeDict); 
                dGammaikl = df_example_pconjpy(Gamma(ii,k,ll),...
                    CVar, derivativeDict);
                partOne = dGammajkl(ii) - dGammaikl(j);
                partTwo = 0;
                for m=1:3
                    partTwo = partTwo + Gamma(j,k,m)*Gamma(ii,m,ll)...
                        - Gamma(ii,k,m)*Gamma(j,m,ll);
                end
                temp = partOne + partTwo;
                temp = myRealVariableFun(temp, realVariable);
                Riem(ii,j,k,ll) = complex_simple3(temp, MVar);
            end
        end
    end
end
clearvars ii j k ll m temp partOne partTwo dGammajkl dGammaikl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Data_example4_v1.mat');