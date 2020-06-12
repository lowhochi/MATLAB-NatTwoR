% WeylCurv_Part52_41to80_NatTwo_CR_rmW.m
WF4180 = sym('WF4180',[40, 5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataZeroV2_WeylChern_Part52_Jun15.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=1:10
    jj = j+40;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF4180(j,:) = [m,n,k,ll,temp];
end
timer5 = toc;
clearvars m n k ll j jj temp

save('Data_WC41to50_Jun14.mat');

%%
load('Data_WC41to50_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=11:20
    jj = j+40;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF4180(j,:) = [m,n,k,ll,temp];
end
timer6 = toc;
clearvars m n k ll j jj temp

save('Data_WC51to60_Jun14.mat');

%%
load('Data_WC51to60_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=21:30
    jj = j+40;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF4180(j,:) = [m,n,k,ll,temp];
end
timer7 = toc;
clearvars m n k ll j jj temp

save('Data_WC61to70_Jun14.mat');

%%
load('Data_WC61to70_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=31:40
    jj = j+40;
    m = WeylFf(jj,1);
    n = WeylFf(jj,2);
    k = WeylFf(jj,3);
    ll = WeylFf(jj,4);
    temp = WeylFf(jj,5);
    temp = complex_simple3(temp, MVarFive2);
    WF4180(j,:) = [m,n,k,ll,temp];
end
timer8 = toc;
clearvars m n k ll j jj temp

save('Data_WC71to80_Jun14.mat');

save('Data_WF4180.mat', 'WF4180', 'timer5', 'timer6', 'timer7', 'timer8');
