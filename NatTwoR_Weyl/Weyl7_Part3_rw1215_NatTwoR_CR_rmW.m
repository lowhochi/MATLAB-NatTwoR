% Weyl7_Part3_rw1215_NatTwoR_CR_rmW.m
load('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat');
temp = group3Vec(2);
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
Wf1215t = temp;

testPart1 = 2*i*y0*(u+conj(u)^3)*(cotton(1,1,2)-cotton(2,1,1));
testPart2 = i/2*y0*(u*conj(u)^3-3*conj(u)^2+3*u*conj(u)-1)...
    *(cotton(2,1,2)-cotton(1,2,2)+cotton(1,3,3)-cotton(3,1,3));
testPart3 = 1/2*y0*(u*conj(u)^3-3*conj(u)^2-3*u*conj(u)+1)...
    *(cotton(3,1,1)-cotton(1,1,3)+cotton(2,2,3)-cotton(3,2,2));
testPart4 = -y0*(u-3*conj(u)+3*u*conj(u)^2-conj(u)^3)*cotton(1,2,3)...
    +2*y0*(u-conj(u)^3)*cotton(2,1,3)...
    -y0*(u+3*conj(u)-3*u*conj(u)^2-conj(u)^3)*cotton(3,1,2);
 
diff1215 = Wf1215t-testPart1-testPart2-testPart3-testPart4;
diff1215 = subs(diff1215, variableSetBianchi2, subSetBianchi2);
diff1215 = subs(diff1215, variableSetBianchi1, subSetBianchi1);

[termVec1, gVec1] = coeffs(diff1215, d2Gset); %length of gVec = 46
for j=1:length(gVec1)
    disp(gVec1(j));
    if (j~=length(gVec1))
        termVec1(j) = complex_simple3(termVec1(j),u);
        disp(termVec1(j));
    end
end

dGterm1215 = termVec1(end);
try
    [termVec1215,gVec1215] = coeffs(dGterm1215, Gset);
catch
    disp('error');
end
myNum = length(gVec1215);
disp(myNum);
for j=51:myNum
    disp(j);
    termVec1215(j) = complex_simple3(termVec1215(j), u);
    disp(gVec1215(j));
    disp(termVec1215(j));
end

save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'termVec1215', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'gVec1215', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'diff1215', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'Wf1215t', '-append');
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat', 'dGterm1215', '-append');
