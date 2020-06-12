load('DataChern1_Part2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latexFile = fopen('latex_Chern_Oct9.txt','w');
fprintf(latexFile, '\\documentclass{article}[10pt,a4paper] \n');
fprintf(latexFile, '\\usepackage{color} \n');
fprintf(latexFile, '\\include{my_page_style}\n');
fprintf(latexFile, '\\begin{document}\n');
fprintf(latexFile, '\\include{codefile} \n');	
fprintf(latexFile, '\\thispagestyle{empty} \n');
fprintf(latexFile, '%s\n', ' ');

fprintf(latexFile, '\\hfill {\\it Date: %s}\n', string(date));
fprintf(latexFile, '\\begin{center}');
fprintf(latexFile, '{\\bf The coefficients $C_{1\\bar1 2\\bar1}$}');
fprintf(latexFile, 'and $C_{1\\bar1 1\\bar1}$ \n');
fprintf(latexFile, '\\end{center} \n');
fprintf(latexFile, '%s\n', ' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cy1121 = Y0*Chern_in_aV(1,1,2,1);
Cy1111 = Y0*Chern_in_aV(1,1,1,1);
[termVec1,yVec1]= coeffs(Cy1121,[Y0]);
[termVec2,yVec2]= coeffs(Cy1111,[Y0]);

fprintf(latexFile,'$Y\\cdot C_{1\\bar1 2\\bar1}$ is: \\\\ \n');
fprintf(latexFile, '%s\n', ' ');
for j=1:length(yVec1)
    yStr = latex(yVec1(j));
    termStr = latex(termVec1(j));
    fprintf(latexFile,'term: \\spa $ %s $: \\\\ \n', yStr);
    fprintf(latexFile, '%s\n', ' ');
    fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n',termStr);
    fprintf(latexFile, '%s\n', ' ');
end
fprintf(latexFile, '\\newpage');
fprintf(latexFile, '%s\n', ' ');
fprintf(latexFile,'$Y\\cdot C_{1\\bar1 1\\bar1}$ is: \\\\ \n');
fprintf(latexFile, '%s\n', ' ');
for j=1:length(yVec2)
    yStr = latex(yVec2(j));
    termStr = latex(termVec2(j));
    fprintf(latexFile,'term: \\spa $ %s $: \\\\ \n', yStr);
    fprintf(latexFile, '%s\n', ' ');
    fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n',termStr);
    fprintf(latexFile, '%s\n', ' ');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(latexFile, '\\end{document}]\n');
fclose(latexFile);