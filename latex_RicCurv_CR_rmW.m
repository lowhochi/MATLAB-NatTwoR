load('DataRicCurv2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latexFile = fopen('latex_RicCurv2_updated.tex','w');
fprintf(latexFile, '\\documentclass{article}[10pt,a4paper] \n');
fprintf(latexFile, '\\usepackage{color} \n');
fprintf(latexFile, '\\include{my_page_style}\n');
fprintf(latexFile, '\\begin{document}\n');
fprintf(latexFile, '\\include{codefile} \n');	
fprintf(latexFile, '\\thispagestyle{empty} \n');
fprintf(latexFile, '%s\n', ' ');

fprintf(latexFile, '\\hfill {\\it Date: %s}\n', string(date));
fprintf(latexFile, '\\begin{center}');
fprintf(latexFile, '{\\bf The Ricci curvature tensor}');
fprintf(latexFile, '\\end{center} \n');
fprintf(latexFile, '%s\n', ' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTENT

RicCurvY = sym('RicY',[6 6]);
RicCurvY(1,1) = Y^2*RicCurv_in_aV(1,1); 
RicCurvY(1,2) = Y^2*RicCurv_in_aV(1,2);
RicCurvY(1,3) = Y^2*RicCurv_in_aV(1,3);
RicCurvY(1,4) = Y^2*RicCurv_in_aV(1,4);
RicCurvY(1,5) = Y^3*RicCurv_in_aV(1,5);
RicCurvY(1,6) = Y*RicCurv_in_aV(1,6);
RicCurvY(2,2) = Y^2*RicCurv_in_aV(2,2);
RicCurvY(2,3) = Y^2*RicCurv_in_aV(2,3);
RicCurvY(2,4) = Y^2*RicCurv_in_aV(2,4);
RicCurvY(2,5) = Y^3*RicCurv_in_aV(2,5);
RicCurvY(2,6) = Y*RicCurv_in_aV(2,6);
RicCurvY(3,3) = Y^2*RicCurv_in_aV(3,3);
RicCurvY(3,4) = Y^2*RicCurv_in_aV(3,4);
RicCurvY(3,5) = Y^3*RicCurv_in_aV(3,5);
RicCurvY(3,6) = Y*RicCurv_in_aV(3,6);
RicCurvY(4,4) = Y^2*RicCurv_in_aV(4,4);
RicCurvY(4,5) = Y^3*RicCurv_in_aV(4,5);
RicCurvY(4,6) = Y*RicCurv_in_aV(4,6);
RicCurvY(5,5) = Y^4*RicCurv_in_aV(5,5); %updated
RicCurvY(5,6) = Y^2*RicCurv_in_aV(5,6);
RicCurvY(6,6) = RicCurv_in_aV(6,6);

yPower = [2, 2, 2, 2, 3, 1;
    0, 2, 2, 2, 3, 1;
    0, 0, 2, 2, 3, 1;
    0, 0, 0, 2, 3, 1;
    0, 0, 0, 0, 4, 2;
    0, 0, 0, 0, 0, 0];

for m=1:6
    for n=1:6
        if m>n
            continue
        end
        
        temp = RicCurv_in_aV(m,n);
        latexTemp = latex(temp);
        if (length(latexTemp)<=1000)
            fprintf(latexFile,'$\\mbox{Ric}(u_%d,u_%d) $ is: \\\\ \n',[m,n]);
            fprintf(latexFile, '%s\n', ' ');
            fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n', latex(temp));
            fprintf(latexFile, '%s\n', ' ');
            continue
        end
        
        [termVec, yVec] = coeffs(RicCurvY(m,n),[Y,conj(Y)]);
        c0 = yPower(m,n);
        fprintf(latexFile,'$Y^%d \\cdot\\mbox{Ric}(u_%d,u_%d) $ is: \\\\ \n',...
            [c0,m,n]);
        fprintf(latexFile, '%s\n', ' ');
        for j=1:length(yVec)
            fprintf(latexFile,'term: \\spa $ %s $: \\\\ \n',latex(yVec(j)));
            fprintf(latexFile, '%s\n', ' ');
            fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n',latex(termVec(j)));
            fprintf(latexFile, '%s\n', ' ');
        end
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(latexFile, '\\end{document}]\n');
fclose(latexFile);

% command = 'pdflatex latex_RicCurv2_updated.tex';
% [status,cmdout] = system(command);