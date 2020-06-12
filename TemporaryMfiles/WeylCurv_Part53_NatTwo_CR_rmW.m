load('Data_WC31to40_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
load('Data_WF4180.mat'); % WF4180
load('Data_WF81120.mat'); % WF81120

for j=41:80
    jj = j-40;
    temp = WF4180(jj,5);
    WeylFf(j,5) = temp;
end

for j=81:120
    jj = j-80;
    temp = WF81120(jj,5);
    WeylFf(j,5) = temp;
end
% save('DataTemp_Jun22.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variableGijk_NatTwo_CR_rmW

WeylFf2 = sym('WeylFf2',[120,5]);
for j=1:120
    m = WeylFf(j,1);
    n = WeylFf(j,2);
    k = WeylFf(j,3);
    ll = WeylFf(j,4);
    temp = WeylFf(j,5);
    temp = subs(temp, dG12_3N, dG12_3Vec);
    temp = subs(temp, dG23_1N, dG23_1Vec);
    temp = subs(temp, dG31_2N, dG31_2Vec);
    temp = subs(temp, dG11_2N, dG11_2Vec);
    temp = subs(temp, dG11_3N, dG11_3Vec);
    temp = subs(temp, dG22_1N, dG22_1Vec);
    temp = subs(temp, dG22_3N, dG22_3Vec);
    temp = subs(temp, dG33_1N, dG33_1Vec);
    temp = subs(temp, dG33_2N, dG33_2Vec);
    temp = subs(temp, d2G12_3N, d2G12_3Vec);
    temp = subs(temp, d2G23_1N, d2G23_1Vec);
    temp = subs(temp, d2G31_2N, d2G31_2Vec);
    temp = subs(temp, d2G11_2N, d2G11_2Vec);
    temp = subs(temp, d2G11_3N, d2G11_3Vec);
    temp = subs(temp, d2G22_1N, d2G22_1Vec);
    temp = subs(temp, d2G22_3N, d2G22_3Vec);
    temp = subs(temp, d2G33_1N, d2G33_1Vec);
    temp = subs(temp, d2G33_2N, d2G33_2Vec);  
    % %
    temp = subs(temp, conj_Gijk, Gijk);
    temp = subs(temp, conj_dGijkE, dGijkE);
    temp = subs(temp, conj_d2GijkE, d2GijkE);
    WeylFf2(j,:) = [m,n,k,ll,temp];
end

remove_variable_Jun28

save('DataZero_WeylCurv53_Jun29.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%