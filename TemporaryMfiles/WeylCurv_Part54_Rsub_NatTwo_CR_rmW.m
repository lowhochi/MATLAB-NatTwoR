load('DataTemp_July15.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
clearvars timer1 timer2 timer3 timer4 timer5 timer6
clearvars timer7 timer8 timer9 timer10 timer11 timer12
clearvars WF4180 WF81120
clearvars WFp53_4180 WFp53_81120
clearvars MVarFive MVarFive2
% Check out symbolic variables in WeylFf2.
% Detect the dG-terms and list them out. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms rm1212 rm1213 rm1223 rm1313 rm1323 rm2323
dGijkN = [dG22_1_by1, dG11_3_by2, dG22_3_by1, dG33_1_by1, dG33_2_by1,...
    dG33_2_by2, dG11_2_by3, dG22_1_by3, dG33_1_by2];

dG22_1_by1Sub = -rm1212 - dG11_2_by2 + G23_1*G12_3 + G11_3*G22_3 ...
    - G23_1*G31_2 + G22_1*G22_1 + G11_2*G11_2 - G12_3*G31_2;

dG11_3_by2Sub = -rm1213 -dG23_1_by1 - G22_1*G12_3 - G11_2*G22_3 ...
    + G22_1*G23_1 + G23_1*G33_1 + G11_2*G11_3 + G12_3*G33_1;

dG22_3_by1Sub = rm1223 + dG12_3_by2 - G22_1*G11_3 + G11_2*G23_1 ...
    + G22_1*G22_3 - G23_1*G33_2 - G11_2*G12_3 - G12_3*G33_2; %updated

dG33_1_by1Sub = rm1313 + dG11_3_by3 - G31_2*G12_3 - G11_2*G33_2 ...
    + G31_2*G23_1 - G33_1*G33_1 - G11_3*G11_3 + G12_3*G23_1; 

dG33_2_by1Sub = -rm1323 - dG12_3_by3 - G31_2*G11_3 - G11_2*G33_1 ...
    + G31_2*G22_3 + G33_1*G33_2 + G11_3*G12_3 + G12_3*G22_3;

dG33_2_by2Sub = -rm2323 - dG22_3_by3 + G31_2*G23_1 + G22_1*G33_1 ...
    - G31_2*G12_3 + G33_2*G33_2 - G23_1*G12_3 + G22_3*G22_3;

dG11_2_by3Sub = -rm1213 + dG31_2_by1 + G33_1*G12_3 - G11_3*G33_2 ...
    - G31_2*G22_1 - G33_1*G31_2 + G11_3*G11_2 - G12_3*G22_1;

dG22_1_by3Sub = rm1223 - dG31_2_by2 - G33_1*G22_3 - G23_1*G33_2 ...
    + G31_2*G11_2 + G33_2*G31_2 + G23_1*G11_2 + G22_3*G22_1; %updated

dG33_1_by2Sub = -rm1323 + dG23_1_by3 + G31_2*G22_3 - G22_1*G33_2 ...
    - G31_2*G11_3 + G33_2*G33_1 - G23_1*G11_3 - G22_3*G23_1;

dGijkSub = [dG22_1_by1Sub, dG11_3_by2Sub, dG22_3_by1Sub, dG33_1_by1Sub,...
    dG33_2_by1Sub, dG33_2_by2Sub, dG11_2_by3Sub, dG22_1_by3Sub,...
    dG33_1_by2Sub];

MVarFive4 = setdiff(MVarFive3, d2GijkE);
MVarFive4 = setdiff(MVarFive4, dGijkN);
MVarFive4 = union(MVarFive4, [rm1212,rm1213,rm1223,rm1313,rm1323,rm2323]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rewrite here. 
% tempVecSix = [Weyl(1,2,5,6), Weyl(1,5,2,6), Weyl(1,5,5,6), Weyl(1,6,2,5),
%   Weyl(2,5,5,6)];
tempVecSix = [WeylFf2(15,5), WeylFf2(48,5), WeylFf2(54,5), ...
    WeylFf2(58,5), WeylFf2(92,5)];
for j=1:5
    temp = tempVecSix(j);
    temp = subs(temp, dGijkN, dGijkSub);
    temp = complex_simple3(temp, MVarFive4);
    tempVecSix(j) = temp;
end
clearvars j temp
clearvars dGijkN dG22_1_by1Sub dG11_3_by2Sub dG22_3_by1Sub
clearvars dG33_1_by1Sub dG33_2_by1Sub dG33_2_by2Sub dG11_2_by3Sub
clearvars dG22_1_by3Sub dG33_1_by2Sub dGijkSub
save('DataTemp2_VecSix_Aug9.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
load('DataTemp2_VecSix_Aug9.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
W1556 = tempVecSix(3);
W1556 = expand(W1556);
varSet1556 = symvar(W1556);

Gset1556 = [G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3, G33_1, G33_2];
dGset1556 = [dG12_3_by1, dG12_3_by2, dG12_3_by3, dG23_1_by1, dG23_1_by2, dG23_1_by3, ...
    dG31_2_by1, dG31_2_by2, dG31_2_by3, dG11_2_by1, dG11_2_by2, dG11_3_by1, dG11_3_by3, ...
    dG22_1_by2, dG22_3_by2, dG22_3_by3, dG33_1_by3, dG33_2_by3];
length_of_dG1556 = length(dGset1556);
rmSet1556 = [rm1212, rm1213, rm1223, rm1313, rm1323, rm2323];

% W1556part1 = the dGij_k-component of W1556
W1556part1 = subs(W1556, Gset1556, zeros(1,9));
W1556part1 = subs(W1556part1, rmSet1556, zeros(1,6));
% W1556part2 = the Gijk-Gijk-component of W1556
W1556part2 = subs(W1556, dGset1556, zeros(1,length_of_dG1556));
W1556part2 = subs(W1556part2, rmSet1556, zeros(1,6));
% W1556part3 = the rm-component of W1556
W1556part3 = subs(W1556, Gset1556, zeros(1,9)); 
W1556part3 = subs(W1556part3, dGset1556, zeros(1,length_of_dG1556));

W1556v0 = W1556part1 + W1556part2 + W1556part3; 
% test1 = W1556 - W1556v0
test1 = W1556 - W1556v0;
test1 = complex_simple3(test1, MVarFive4);

W1556part1 = complex_simple3(W1556part1, MVarFive4); % = 0
W1556part2 = complex_simple3(W1556part2, MVarFive4);
W1556part3 = complex_simple3(W1556part3, MVarFive4); % = 0

[term1556, u1556] = coeffs(W1556part2, Gset1556);
term1556v2 = sym('term1556v2',[1,45]);
for j = 1:45
    term1556v2(j) = complex_simple3(term1556(j), [u]);
end
save('DataTemp3_VecSix_Aug19.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
