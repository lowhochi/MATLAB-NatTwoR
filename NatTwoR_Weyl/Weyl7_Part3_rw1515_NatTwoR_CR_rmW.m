% Weyl7_Part3_rw1515_NatTwoR_CR_rmW.m
load('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat');

variable_NatTwoR_W725_CR_rmW

temp = group3Vec(3);
temp = subs(temp, d2G12_3N, d2G12_3Vec);
temp = subs(temp, d2G23_1N, d2G23_1Vec);
temp = subs(temp, d2G31_2N, d2G31_2Vec);
temp = subs(temp, d2G11_2N, d2G11_2Vec);
temp = subs(temp, d2G11_3N, d2G11_3Vec);
temp = subs(temp, d2G22_1N, d2G22_1Vec);
temp = subs(temp, d2G22_3N, d2G22_3Vec);
temp = subs(temp, d2G33_1N, d2G33_1Vec);
temp = subs(temp, d2G33_2N, d2G33_2Vec);
temp = subs(temp, variableSetMu, subSetMu);
temp = subs(temp, variableSetConjMu, subSetConjMu);
temp = subs(temp, variableSetVnormv, subSetVnormv);
temp = subs(temp, variableSetBianchi2, subSetBianchi2);
temp = subs(temp, variableSetBianchi1, subSetBianchi1);
Wf1515t = temp;

% % NC substitution
% Wf1515t = subs(Wf1515t, cottonVariableSet3, cottonSubSet3);
% Wf1515t = subs(Wf1515t, cottonVariableSet2, cottonSubSet2);
% Wf1515t = subs(Wf1515t, cottonVariableSet1, cottonSubSet1);
% % Wf1515t = subs(WF1515t, d2G12_3_by23, d2G12_3_by23Sub);
% % Wf1515t = subs(WF1515t, d2G33_2_by11, d2G33_2_by11Sub);
% Wf1515t = subs(Wf1515t, d2G23_1_by12, cotton223 -d2G11_3_by22 - term223(3));
% 
% for j=1:length(Gset)
%     try 
%         Wf1515t = subs(Wf1515t, Gset(j), 0);
%     catch
%         disp(Gset(j));
%     end
% end
% Wf1515t = complex_simple3(Wf1515t, u);
% [termNC, gVecNC] = coeffs(Wf1515t, totalSet);
% for j=1:length(gVecNC)
%     termNC(j) = complex_simple3(termNC(j),u);
%     disp(gVecNC(j));
%     disp(termNC(j));
% end

testPart1 = i*(1-conj(u)^4)*(cotton(1,1,2)-cotton(2,1,1));
testPart2 = (conj(u)-conj(u)^3)...
    *(cotton(1,1,3)-cotton(3,1,1)+cotton(3,2,2)-cotton(2,2,3));
testPart3 =i*(conj(u)^3+conj(u))...
    *(cotton(1,3,3)-cotton(3,1,3)+cotton(2,1,2)-cotton(1,2,2));
testPart4 = -(3*conj(u)^2+conj(u)^4/2+1/2)*cotton(1,2,3)...
    +(conj(u)^4+1)*cotton(2,1,3)...
    +(3*conj(u)^2-conj(u)^4/2-1/2)*cotton(3,1,2);

diff1515 = Wf1515t-testPart1-testPart2-testPart3-testPart4;
diff1515 = subs(diff1515, variableSetBianchi2, subSetBianchi2);
diff1515 = subs(diff1515, variableSetBianchi1, subSetBianchi1);

[termVec, gVec] = coeffs(diff1515, d2Gset); %length of gVec = 46
dGterm1515 = termVec(end);
for j=1:length(gVec)
    disp(gVec(j));
    if (j~=length(gVec))
        termVec(j) = complex_simple3(termVec(j),u);
        disp(termVec(j));
    end
end

try 
    [termVec1515, gVec1515] = coeffs(dGterm1515, Gset);
catch
    disp('Error');
end
disp(length(gVec1515)); %length = 381
% extraTerm1515 = 0;
for j= 1:length(gVec1515)
    disp(j);
    termVec1515(j) = complex_simple3(termVec1515(j),u);
    disp(termVec1515(j));
%     if (termVec1515(j)~=0)
%         extraTerm1515 = extraTerm1515 + termVec1515(j)*gVec1515(j);
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'termVec1515', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'gVec1515', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'diff1515', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'Wf1515t', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'dGterm1515', '-append');




