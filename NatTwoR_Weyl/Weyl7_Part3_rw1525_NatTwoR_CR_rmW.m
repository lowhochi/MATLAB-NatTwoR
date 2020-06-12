% Weyl7_Part3_rw15 25_NatTwoR_CR_rmW.m
load('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat');
variable_NatTwoR_W725_CR_rmW

temp = group3Vec(4);
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
Wf1525t = temp;

testPart1 = i*(u^2-conj(u)^2)*(cotton(1,1,2)-cotton(2,1,1));
testPart2 = (u*conj(u)^2/2 +u^2*conj(u)/2 -conj(u)/2 -u/2)...
    *(cotton(1,1,3) -cotton(3,1,1) +cotton(3,2,2) -cotton(2,2,3));
testPart3 = i/2*(u-conj(u)+u*conj(u)^2-u^2*conj(u))...
    *(cotton(1,2,2) -cotton(2,1,2) +cotton(3,1,3) -cotton(1,3,3));
testPart4 = 1/2*(4*u*conj(u)-u^2*conj(u)^2-u^2-conj(u)^2-1)*cotton(1,2,3)...
    +(u^2 +conj(u)^2)*cotton(2,1,3)...
    +1/2*(u^2*conj(u)^2-conj(u)^2-u^2-4*u*conj(u)+1)*cotton(3,1,2);

diff1525 = Wf1525t-testPart1-testPart2-testPart3-testPart4;
diff1525 = subs(diff1525, variableSetBianchi2, subSetBianchi2);
diff1525 = subs(diff1525, variableSetBianchi1, subSetBianchi1);

[termVec, gVec] = coeffs(diff1525, d2Gset); %length of gVec = 46
dGterm1525 = termVec(end);
for j=1:length(gVec)
    disp(gVec(j));
    if (j~=length(gVec))
        termVec(j) = complex_simple3(termVec(j),u);
        disp(termVec(j));
    end
end

try 
    [termVec1525, gVec1525] = coeffs(dGterm1525, Gset);
catch
    disp('Error');
end
disp(length(gVec1525)); %length = 381
for j= 1:length(gVec1525)
    disp(j);
    termVec1525(j) = complex_simple3(termVec1525(j),u);
    disp(termVec1525(j));
end

save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'termVec1525', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'gVec1525', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'diff1525', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'Wf1525t', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'dGterm1525', '-append');