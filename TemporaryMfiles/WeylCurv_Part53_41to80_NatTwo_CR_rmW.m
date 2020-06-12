% WeylCurv_Part53_41to80_NatTwo_CR_rmW.m
WFp53_4180 = sym('WFp53_4180',[40, 5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataZero_WeylCurv53_Jun29.mat');
assumeAlso([x y z u1 u2 gamma], 'real');

tic
for j=1:10
    jj = j+40;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_4180(j,:) = [m,n,k,ll,temp];
end
timer5 = toc;
clearvars m n k ll j jj temp

save('Data_WCp3_41to50_Jun30.mat');

%%
load('Data_WCp3_41to50_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=11:20
    jj = j+40;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_4180(j,:) = [m,n,k,ll,temp];
end
timer6 = toc;
clearvars m n k ll j jj temp

save('Data_WCp3_51to60_Jun30.mat');

%%
load('Data_WCp3_51to60_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=21:30
    jj = j+40;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_4180(j,:) = [m,n,k,ll,temp];
end
timer7 = toc;
clearvars m n k ll j jj temp

save('Data_WCp53_61to70_Jun30.mat');

%%
load('Data_WCp53_61to70_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=31:40
    jj = j+40;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_4180(j,:) = [m,n,k,ll,temp];
end
timer8 = toc;
clearvars m n k ll j jj temp

save('Data_WCp53_71to80_Jun30.mat');

save('Data_WFp53_4180.mat', 'WFp53_4180', 'timer5', 'timer6', 'timer7', 'timer8');
