% WeylCurv_Part52_1to40_NatTwo_CR_rmW.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataZeroV2_WeylChern_Part52_Jun15.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=1:10
    temp = WeylFf(j,5);
    temp = complex_simple3(temp, MVarFive2);
    WeylFf(j,5) = temp;
end
timer1 = toc;
clearvars j temp

save('Data_WC1to10_Jun14.mat');

%%
load('Data_WC1to10_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=11:20
    temp = WeylFf(j,5);
    temp = complex_simple3(temp, MVarFive2);
    WeylFf(j,5) = temp;
end
timer2 = toc;
clearvars j temp

save('Data_WC11to20_Jun14.mat');

%%
load('Data_WC11to20_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=21:30
    temp = WeylFf(j,5);
    temp = complex_simple3(temp, MVarFive2);
    WeylFf(j,5) = temp;
end
timer3 = toc;
clearvars j temp

save('Data_WC21to30_Jun14.mat');

%%
load('Data_WC21to30_Jun14.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
tic
for j=31:40
    temp = WeylFf(j,5);
    temp = complex_simple3(temp, MVarFive2);
    WeylFf(j,5) = temp;
end
timer4 = toc;
clearvars j temp

save('Data_WC31to40_Jun14.mat');