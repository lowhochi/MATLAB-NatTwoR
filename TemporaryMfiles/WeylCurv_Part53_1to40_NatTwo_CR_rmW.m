% WeylCurv_Part53_1to40_NatTwo_CR_rmW.m
load('DataZero_WeylCurv53_Jun29.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for j=1:10
    temp = WeylFf2(j,5);
    temp = complex_simple3(temp, MVarFive3);
    WeylFf2(j,5) = temp;
end
timer1 = toc;
clearvars j temp

save('Data_WCp53_1to10_Jun30.mat');

%%
load('Data_WCp53_1to10_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=11:20
    temp = WeylFf2(j,5);
    temp = complex_simple3(temp, MVarFive3);
    WeylFf2(j,5) = temp;
end
timer2 = toc;
clearvars j temp

save('Data_WCp53_11to20_Jun30.mat');

%%
load('Data_WCp53_11to20_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=21:30
    temp = WeylFf2(j,5);
    temp = complex_simple3(temp, MVarFive3);
    WeylFf2(j,5) = temp;
end
timer3 = toc;
clearvars j temp

save('Data_WCp53_21to30_Jun30.mat');

%%
load('Data_WCp53_21to30_Jun30.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=31:40
    temp = WeylFf2(j,5);
    temp = complex_simple3(temp, MVarFive3);
    WeylFf2(j,5) = temp;
end
timer4 = toc;
clearvars j temp

save('Data_WCp53_31to40_Jun30.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%