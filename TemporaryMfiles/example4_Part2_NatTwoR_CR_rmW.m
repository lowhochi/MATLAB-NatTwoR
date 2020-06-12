% example4_Part2_NatTwoR_CR_rmW.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataTemp_example4_v1.mat');
Ric = sym('Ric',[3 3]);
for ii=1:3
    for j=1:3
        temp = 0;
        for k=1:3
            temp = temp + Riem(k,ii,j,k);
        end
        temp = myRealVariableFun(temp, realVariable);
        Ric(ii,j) = complex_simple3(temp, MVar);
    end
end
S = 0;
for ii=1:3
    for j=1:3
        S = S + Ric(ii,j)*d2phiInv(j,ii);
    end
end
S = myRealVariableFun(S, realVariable);
S = complex_simple3(S, MVar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schouten tensor: Pij
P = sym('P',[3 3]);
Cotton = sym('Cotton',[3 3 3]);
for ii=1:3
    for j=1:3
        P(ii,j) = Ric(ii,j)-S/4*d2phi(ii,j);
        P(ii,j) = complex_simple3(P(ii,j),MVar);
    end
end
% % % % %
for ii=1:3
    for j=1:3
        for k=1:3
            dPij = df_example_pconjpy(P(ii,j),CVar,derivativeDict); 
            dPik = df_example_pconjpy(P(ii,k),CVar,derivativeDict); 
            partOne = dPij(k) - dPik(j); 
            partTwo = 0;
            for ll=1:3
                partTwo = Gamma(j,ii,ll)*P(ll,k)+Gamma(j,k,ll)*P(ii,ll)...
                    -Gamma(k,ii,ll)*P(ll,j)-Gamma(k,j,ll)*P(ii,ll);
            end
            temp = partOne + partTwo;
            temp = myRealVariableFun(temp, realVariable);
            Cotton(ii,j,k) = complex_simple3(temp, MVar);
        end
    end
end
clearvars ii j k partOne partTwo temp dPij dPik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('DataTemp_example4_v2.mat');