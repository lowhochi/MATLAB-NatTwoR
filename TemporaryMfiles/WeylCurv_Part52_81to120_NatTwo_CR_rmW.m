% WeylCurv_Part52_81to120_NatTwo_CR_rmW.m
WF81120 = sym('WF81120',[40, 5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataZeroV2_WeylChern_Part52_Jun15.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=1:10
    jj = j+80;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF81120(j,:) = [m,n,k,ll,temp];
end
timer9 = toc;
clearvars m n k ll j jj temp

save('Data_WC81to90_Jun14.mat');

%%
load('Data_WC81to90_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=11:20
    jj = j+80;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF81120(j,:) = [m,n,k,ll,temp];
end
timer10 = toc;
clearvars m n k ll j jj temp

save('Data_WC91to100_Jun14.mat');

%%
load('Data_WC91to100_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=21:30
    jj = j+80;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF81120(j,:) = [m,n,k,ll,temp];
end
timer11 = toc;
clearvars m n k ll j jj temp

save('Data_WC101to110_Jun14.mat');

%%
load('Data_WC101to110_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=31:40
    jj = j+80;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF81120(j,:) = [m,n,k,ll,temp];
end
timer12 = toc;
clearvars m n k ll j jj temp

save('Data_WC111to120_Jun14.mat');

save('Data_WF81120.mat','WF81120','timer9','timer10','timer11','timer12');