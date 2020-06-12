% WeylCurv_Part53_81to120_NatTwo_CR_rmW.m
WFp53_81120 = sym('WFp53_81120',[40, 5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataZero_WeylCurv53_Jun29.mat');
assumeAlso([x y z u1 u2 gamma], 'real');

tic
for j=1:10
    jj = j+80;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_81120(j,:) = [m,n,k,ll,temp];
end
timer9 = toc;
clearvars m n k ll j jj temp

save('Data_WCp3_81to90_Jun30.mat');

%%
load('Data_WCp3_81to90_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=11:20
    jj = j+80;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_81120(j,:) = [m,n,k,ll,temp];
end
timer10 = toc;
clearvars m n k ll j jj temp

save('Data_WCp3_91to100_Jun30.mat');

%%
load('Data_WCp3_91to100_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=21:30
    jj = j+80;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_81120(j,:) = [m,n,k,ll,temp];
end
timer11 = toc;
clearvars m n k ll j jj temp

save('Data_WCp3_101to110_Jun30.mat');

%%
load('Data_WCp3_101to110_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=31:40
    jj = j+80;
    m = WeylFf2(jj,1);
    n = WeylFf2(jj,2);
    k = WeylFf2(jj,3);
    ll = WeylFf2(jj,4);
    temp = WeylFf2(jj,5);
    temp = complex_simple3(temp, MVarFive3);
    WFp53_81120(j,:) = [m,n,k,ll,temp];
end
timer12 = toc;
clearvars m n k ll j jj temp

save('Data_WCp3_111to120_Jun30.mat');

save('Data_WFp53_81120.mat','WFp53_81120','timer9','timer10','timer11','timer12');