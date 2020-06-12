% Weyl4_NatTwoR_CR_rmW.m
load('DataWeyl2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part (1) Weyl=0 on the following (m,n,k,ll)-index.
% They are from {X1,conj(X1),X2,conj(X2)} but excluded from 'W=C' part.
% [1,2,1,3; 1,2,2,4; 1,3,1,3; 1,3,1,4;
%    1,3,2,3; 1,3,2,4; 1,3,3,4; 1,4,2,4; 2,3,2,4; 2,4,2,4; 2,4,3,4];
% 11 out of 21 combinations within {X1,conj(X1),X2,conj(X2)}.

indexPart1 = [1,2,1,3; 1,2,2,4; 1,3,1,3; 1,3,1,4;
    1,3,2,3; 1,3,2,4; 1,3,3,4; 1,4,2,4; 2,3,2,4; 2,4,2,4; 2,4,3,4];
for j=1:11
    m = indexPart1(j,1);
    n = indexPart1(j,2);
    k = indexPart1(j,3);
    ll =indexPart1(j,4);
    myNumber = 1;
    while (myNumber<=120)
        rowTemp = WeylTwo(myNumber,1:4);
        if (rowTemp==[m,n,k,ll])
            myString = sprintf('row number: %d ',myNumber);
            disp(myString);
            disp(WeylTwo(myNumber,1:4));
            disp(WeylTwo(myNumber,5));
            break
        end
        myNumber = myNumber+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part (2)
indexWithFive = []; %size = [34,4];
for j=1:120
    m = countIndex120(j,1);
    n = countIndex120(j,2);
    k = countIndex120(j,3);
    ll = countIndex120(j,4);
    mySet = [m,n,k,ll];
    if ismember(6,mySet)==1
        continue
    end
    if ismember(5,mySet)==1
        indexWithFive = [indexWithFive; mySet];
    end
end

indexPart2 = []; %size = [15,4];
for j=1:34
    m = indexWithFive(j,1);
    n = indexWithFive(j,2);
    k = indexWithFive(j,3);
    ll = indexWithFive(j,4);
    mySet = [m,n,k,ll];
    myNumber = sum(ismember(mySet,3)) + sum(ismember(mySet,4));
    if myNumber>1
        indexPart2 = [indexPart2; mySet]; 
    end
end

WeylPart2 = sym('WeylPart4',[15 5]);
% WeylPart2(j,:) = [indexPart2(j,:), W(indexPart2(j,:))];
for j=1:15
    m = indexPart2(j,1);
    n = indexPart2(j,2);
    k = indexPart2(j,3);
    ll = indexPart2(j,4);
    WeylPart2(j,1:4)=[m,n,k,ll];
    
    myNumber = 1;
    while (myNumber<=120)
        rowTemp = WeylTwo(myNumber,1:4);
        if (rowTemp==[m,n,k,ll])
            WeylPart2(j,5) = WeylTwo(myNumber,5);
            break
        end
        myNumber = myNumber+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% latexFile = fopen('latex_Weyl4Part2.txt','w');
% for j=1:15
%     m = WeylPart2(j,1);
%     n = WeylPart2(j,2);
%     k = WeylPart2(j,3);
%     ll = WeylPart2(j,4);
%     temp = WeylPart2(j,5);
%     fprintf(latexFile,'$Weyl(u_%d,u_%d,u_%d,u_%d)$ is: \\\\ \n',[m,n,k,ll]);
%     fprintf(latexFile, '%s\n', ' ');
%     fprintf(latexFile, '$\\ds %s $ \\\\[0.1in] \n', latex(temp));
%     fprintf(latexFile, '%s\n', ' ');   
% end
% fclose(latexFile);

% Check the values of Weyl curvature coefficients in Part(2).
test1435 = WeylPart2(3,5) +1/24*drho_uSub;
test1435 = complex_simple3(test1435, symSetWeylTwo);

test1445 = WeylPart2(4,5) -1/24*drho_conjuSub; 
test1445 = complex_simple3(test1445, symSetWeylTwo);

test1534 = WeylPart2(5,5) +1/24*drho_uSub; 
test1534 = complex_simple3(test1534, symSetWeylTwo);

test2335 = WeylPart2(6,5) -1/24*drho_uSub;
test2335 = complex_simple3(test2335, symSetWeylTwo);

test2345 = WeylPart2(7,5) +1/24*drho_conjuSub;
test2345 = complex_simple3(test2345, symSetWeylTwo);

test2534 = WeylPart2(10,5)-1/24*drho_conjuSub;
test2534 = complex_simple3(test2534, symSetWeylTwo);

d2rho_uuSubR = subs(d2rho_uuSub,d2aV_conjuconju,d2aV_conjuconjuR);
test3535 = WeylPart2(13,5) +1/48*d2rho_uuSubR +conj(u)/(24*(1+u*conj(u)))*drho_uSub;
test3535 = complex_simple3(test3535, symSetWeylTwo);

test3545 = WeylPart2(14,5) +1/(24*(1+u*conj(u))^2)*rhoSub;
test3545 = complex_simple3(test3545, symSetWeylTwo);

d2rho_conjuconjuSubR = subs(d2rho_conjuconjuSub,d2aV_conjuconju,d2aV_conjuconjuR);
test4545 = WeylPart2(15,5) +1/48*d2rho_conjuconjuSubR +u/(24*(1+u*conj(u)))*drho_conjuSub;
test4545 = complex_simple3(test4545, symSetWeylTwo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part(3) Check the following terms.
W1315 = WeylTwo(18,5);
W1325 = WeylTwo(22,5);
W1524 = WeylTwo(46,5);
W2425 = WeylTwo(77,5);

WeylPart3 = [W1325, W1524]; 

% 'Ric12'-substitution
KtermSub = d2aV_conjuconjMu...
    -(1+u*conj(u))^2*(daV_u*daV_conju)...
    +(1+u*conj(u))^2*conj(daV_u)*conj(daV_conju)...
    -2*daV_vnormv +2*conj(daV_vnormv) +4*i*dtheta_vnormv...
    +4*i*theta*(aV+conj(aV)) +2*conj(aV)^2 - 2*aV^2;

WeylPart3(1) = subs(WeylPart3(1),conj(d2aV_conjuconjMu),KtermSub);
WeylPart3(1) = complex_simple3(WeylPart3(1),symSetWeylTwo);

WeylPart3(2) = subs(WeylPart3(2),conj(d2aV_conjuconjMu),KtermSub);
WeylPart3(2) = complex_simple3(WeylPart3(2),symSetWeylTwo);

variable_NatTwoR_Chern2_CR_rmW % Use the Chern2 file

variableSet0 = [aM, bV, h11, T4];
subSet0 = [aMSub, bVSub, T4Sub, h11SubTwo];
% substitute variable0 BEFORE [theta, aV, daV_u, daV_conju, d2aV_uu].

variableSetG1 = [theta, aV, daV_u, daV_conju, d2aV_uu];
variableSetG2 = [daV_conjMu, daV_vnormv, dtheta_vnormv];
variableSetG3 = [d2aV_uMu, d2aV_uvnormv, d2aV_conjuconjMu, d3aV_uuMu];
thetaSub = G12_3 + G23_1 + G31_2;
dtheta_vnormvSub = dG12_3_vnormv + dG23_1_vnormv + dG31_2_vnormv;

aVSub = 1/(2*(1+u*conj(u))^2)*(i*(mu1*conj(mu1)+mu3*conj(mu3))*G23_1...
    +i*(mu1*conj(mu1)+mu2*conj(mu2))*G31_2...
    +i*(mu2*conj(mu2)+mu3*conj(mu3))*G12_3...
    -i*conj(mu1)*mu3*G11_2 + i*conj(mu1)*mu2*G11_3...
    +i*conj(mu2)*mu3*G22_1 - i*conj(mu2)*mu1*G22_3...
    +i*conj(mu3)*mu1*G33_2 -i*conj(mu3)*mu2*G33_1);

daVVec = df_NatTwo_MuSet_CR_rmW(aVSub, CVarG1, dCVarG1);
daV_conjMuSub = daVVec(2);
daV_vnormvSub = daVVec(3);
daV_uSub = daVVec(4);
daV_conjuSub = daVVec(5);

d2aV_uVec = df_NatTwo_MuSet_CR_rmW(daV_uSub, CVarG1, dCVarG1);
d2aV_uMuSub = d2aV_uVec(1);
d2aV_uvnormvSub = d2aV_uVec(3);
d2aV_uuSub = d2aV_uVec(4);
d2aV_conjuVec = df_NatTwo_MuSet_CR_rmW(daV_conjuSub, CVarG1, dCVarG1);
d2aV_conjuconjMuSub = d2aV_conjuVec(2);

d3aV_uuVec = df_NatTwo_MuSet_CR_rmW(d2aV_uuSub, CVarG1, dCVarG1);
d3aV_uuMuSub = d3aV_uuVec(1);

% variableSetG1 = [theta, aV, daV_u, daV_conju, d2aV_uu];
% variableSetG2 = [daV_conjMu, daV_vnormv, dtheta_vnormv];
% variableSetG3 = [d2aV_uMu, d2aV_uvnormv, d2aV_conjuconjMu, d3aV_uuMu];
subSetG1 = [thetaSub, aVSub, daV_uSub, daV_conjuSub, d2aV_uuSub];
subSetG2 = [daV_conjMuSub, daV_vnormvSub, dtheta_vnormvSub];
subSetG3 = [d2aV_uMuSub, d2aV_uvnormvSub, d2aV_conjuconjMuSub, d3aV_uuMuSub];

for j=1:2
    temp = WeylPart3(j);
    temp = subs(temp, variableSetG3, subSetG3);
    temp = subs(temp, variableSetG2, subSetG2);
    temp = subs(temp, variableSet0, subSet0);
    temp = subs(temp, variableSetG1, subSetG1);
    % Second substitution
    temp = subs(temp, variableSetMu, subSetMu);
    temp = subs(temp, variableSetConjMu, subSetConjMu);
    temp = subs(temp, variableSetVnormv, subSetVnormv);
    % Bianchi Identity
    temp = subs(temp, variableSetBianchi1, subSetBianchi1);
    temp = complex_simple3(temp, [u]);
    WeylPart3(j)=temp;
end