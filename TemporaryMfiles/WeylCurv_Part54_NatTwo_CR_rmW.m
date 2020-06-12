load('Data_WCp53_31to40_Jun30.mat');
load('Data_WFp53_4180.mat');
load('Data_WFp53_81120.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put things back to WeylFf2
for j=41:80
    jj=j-40;
    temp = WFp53_4180(jj,5);
    WeylFf2(j,5) = temp;
end

for j=81:120
    jj=j-80;
    temp = WFp53_81120(jj,5);
    WeylFf2(j,5) = temp;
end

clearvars j jj temp
% save('DataTemp_July15.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symSetWFf2 = cell(120,2);
for j=1:120
    tempSet = symvar(WeylFf2(j,5));
    symSetWFf2{j,1} = WeylFf2(j,1:4);
    symSetWFf2{j,2} = tempSet;
end
clearvars j tempSet

latex_WeylCurv54 = fopen('latex_WeylCurv54.txt','w');
fprintf(latex_WeylCurv54,'Symbolic variables of Weyl Curvature terms\n');
for j=1:120
    m = WeylFf2(j,1);
    n = WeylFf2(j,2);
    k = WeylFf2(j,3);
    ll = WeylFf2(j,4);
    tempLatex = latex(symSetWFf2{j,2});
    fprintf(latex_WeylCurv54, 'index: (%d,%d,%d,%d) \n',[m,n,k,ll]);
    fprintf(latex_WeylCurv54, 'symvar: \n');
    fprintf(latex_WeylCurv54, '%s\n', tempLatex);
    fprintf(latex_WeylCurv54, '%s\n', ' ');
end
fclose(latex_WeylCurv54);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%