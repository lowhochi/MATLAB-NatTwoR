load('DataWeyl7_Part25_NatTwoR_CR_rmW.mat');
variable_NatTwoR_W725_CR_rmW
y0 = 1+u*conj(u);
% Wf1212 = Weylf(1,5);
% Wf1215 = Weylf(4,5);
% Wf1515 = Weylf(43,5);
% Wf1525 = Weylf(47,5);
% group3Vec = [Wf1212, Wf1215, Wf1515, Wf1525];
% Second Substitution (for Wf1212);
temp = group3Vec(1);
    
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
Wf1212t = temp;

% NC substitution
% Wf1212t = subs(Wf1212t, covPVariableSet3, covPSubSet3);
% Wf1212t = subs(Wf1212t, covPVariableSet2, covPSubSet2);
% Wf1212t = subs(Wf1212t, covPVariableSet1, covPSubSet1);
% % d2G23_1_by12 = covP223 -d2G11_3_by22 - term223(3);
% Wf1212t = subs(Wf1212t, d2G23_1_by12, covP223 -d2G11_3_by22 - term223(3)); 
% for j=1:length(totalSet)
%   Wf1212t = subs(Wf1212t,conj(totalSet(j)),totalSet(j)); 
% end
% Wf1212t = subs(Wf1212t, Gset, zeros(1,length(Gset)));
% [term1212, gVec1212] = coeffs(Wf1212t, totalSet);
% for j=1:length(gVec1212)
%     term1212(j) = complex_simple3(term1212(j),u);
%     disp(gVec1212(j));
%     disp(term1212(j));
% end

testPart1 = 4*i*y0^2*(conj(u)^2-u^2)*(covP(2,1,1)-covP(1,1,2));
testPart2 = 2*y0^2*(u+conj(u))*(1-u*conj(u))...
    *(covP(3,1,1)+covP(2,2,3)-covP(1,1,3)-covP(3,2,2));
testPart3 = 2*i*y0^2*(u-conj(u))*(1-u*conj(u))...
    *(covP(1,2,2) -covP(2,1,2) +covP(3,1,3) -covP(1,3,3));
testPart4 = 4*y0^2*(u^2+conj(u)^2)*covP(2,1,3)...
    -2*y0^2*((u-conj(u))^2+(1-u*conj(u))^2)*covP(1,2,3)...
    -2*y0^2*((u+conj(u))^2-(1-u*conj(u))^2)*covP(3,1,2);

diff1212 = Wf1212t -testPart1 - testPart2 - testPart3 - testPart4;
diff1212 = subs(diff1212, variableSetBianchi2, subSetBianchi2);
diff1212 = subs(diff1212, variableSetBianchi1, subSetBianchi1);
% diff1212 = complex_simple3(diff1212, u);
% save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat');
[termVec1, gVec1] = coeffs(diff1212, d2Gset);
for j=1:length(gVec1)
    if j~=length(gVec1);
        termVec1(j) = complex_simple3(termVec1(j), u);
        disp(gVec1(j));
        disp(termVec1(j));
    end
end

dGterm1212 = termVec1(end);
try
    [termVec2,gVec2] = coeffs(dGterm1212, Gset);
catch
    disp('error');
end
myNum = length(gVec2);
disp(myNum);
for j=1:myNum
    disp(j);
    termVec2(j) = complex_simple3(termVec2(j), u);
    disp(gVec2(j));
    disp(termVec2(j));
end

termVec1212 = termVec2;
gVec1212 = gVec2;

clearvars wChar testPart1 testPart2 testPart3 testPart4
clearvars temp temp01 subChar part1 part2 n m ll j 
clearvars list_of_variables_Main2 list_of_variables_Main3 list_of_variables_Main4 
clearvars myChar myCharR myNum k gVec1 termVec1 dricTemp
clearvars termVec2 gVec2
save('DataWeyl7new_Part25_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%