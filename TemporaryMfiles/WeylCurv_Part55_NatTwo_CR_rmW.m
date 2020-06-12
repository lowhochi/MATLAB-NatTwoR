% b/f from WeylCurv_Part53_NatTwo_CR_rmW.m
% Investigate the terms Weyl(1,3,2,5) and Weyl(1,5,2,4). 
load('DataTemp_July15.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
clearvars timer1 timer2 timer3 timer4 timer5 timer6
clearvars timer7 timer8 timer9 timer10 timer11 timer12
clearvars WF4180 WF81120 countSet
clearvars WFp53_4180 WFp53_81120
clearvars MVarFive MVarFive2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W1325 and W1524 are in u, Gij_k's and dGij_k_byll. 
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

Gset1325 = [G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3, G33_1, G33_2];
dGset1325 = [dG12_3_by1, dG12_3_by2, dG12_3_by3, dG23_1_by1, dG23_1_by2, dG23_1_by3, ...
    dG31_2_by1, dG31_2_by2, dG31_2_by3, dG11_2_by1, dG11_2_by2, dG11_3_by1, dG11_3_by3, ...
    dG22_1_by2, dG22_3_by2, dG22_3_by3, dG33_1_by3, dG33_2_by3];
rmSet1325 = [rm1212, rm1213, rm1223, rm1313, rm1323, rm2323];
totalSet1325 = [Gset1325, dGset1325, rmSet1325];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W1325 = subs(WeylFf2(22,5), dGijkN, dGijkSub);
W1325 = expand(W1325);
W1325 = complex_simple3(W1325, MVarFive4);
[term1325, G1325] = coeffs(W1325, totalSet1325);

W1524 = subs(WeylFf2(46,5), dGijkN, dGijkSub);
W1524 = expand(W1524);
W1524 = complex_simple3(W1524, MVarFive4);
[term1524, G1524] = coeffs(W1524, totalSet1325);

for j=1:length(G1325)
   temp1325 = term1325(j);
   temp1325 = complex_simple3(temp1325, [u]);
   term1325(j) = temp1325;
   clear temp1325
end

for j=1:length(G1524)
    temp1524 = term1524(j);
    temp1524 = complex_simple3(temp1524, [u]);
    term1524(j) = temp1524;
    clear temp1524
end

% Check that W1325 and W1524
% W1325part1 = dGij_k_byll component
% W1325part2 = Gij_k*Gij_k component
% W1325part3 = rm(i,j,k,ll) component
W1325part1 = subs(W1325, Gset1325, zeros(1,9));
W1325part1 = subs(W1325part1, rmSet1325, zeros(1,6));
W1325part2 = subs(W1325, dGset1325, zeros(1,18));
W1325part2 = subs(W1325part2, rmSet1325, zeros(1,6));
W1325part3 = subs(W1325, Gset1325, zeros(1,9));
W1325part3 = subs(W1325part3, dGset1325, zeros(1,18));

test1325 = W1325 - (W1325part1 + W1325part2 + W1325part3);
test1325 = complex_simple3(test1325, MVarFive4);

clearvars dGijkN dG22_1_by1Sub dG11_3_by2Sub dG22_3_by1Sub
clearvars dG33_1_by1Sub dG33_2_by1Sub dG33_2_by2Sub dG11_2_by3Sub
clearvars dG22_1_by3Sub dG33_1_by2Sub dGijkSub
clearvars Gset1325 dGset1325 rmSet1325 j

save('DataTemp_WeylCurv_Part55_Aug20.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%