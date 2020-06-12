load('DataRicCurv2_NatTwoR_CR_rmW.mat');
Ry55 = Y^4*RicCurv_in_aV(5,5);
[termVec, yVec] = coeffs(Ry55,[Y,conj(Y)]);

latexFile = fopen('latex_RicCurv55.txt','w');
fprintf(latexFile,'$Y^4\\cdot\\mbox{Ric}(u_5,u_5) $ is: \\\\ \n');
fprintf(latexFile, '%s\n', ' ');

        
for j=1:length(yVec)
    fprintf(latexFile,'term: \\spa $ %s $: \\\\ \n',latex(yVec(j)));
    fprintf(latexFile, '%s\n', ' ');
    fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n',latex(termVec(j)));
    fprintf(latexFile, '%s\n', ' ');
end


fclose(latexFile);