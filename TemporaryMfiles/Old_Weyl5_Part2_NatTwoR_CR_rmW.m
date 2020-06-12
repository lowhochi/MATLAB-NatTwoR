%Weyl5_Part2_NatTwoR_CR_rmW.m
load('DataWeylTwo_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indexW = [1,2,1,5;
    1,5,1,5;
    1,5,2,5];
W1215 = WeylTwo(4,5);
W1515 = WeylTwo(43,5);
W1525 = WeylTwo(47,5);
wVecNL = [W1215, W1515, W1525];
yPower = [2, 4, 4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% latexFile = fopen('latex_wNL_Nov5.txt','w');
% for j=1:3
%      myNumber = yPower(j);
%      temp = (Y^myNumber)*wVecNL(j);
%      [termVec, yVec] = coeffs(temp,[Y,u,conj(u)]);
%      myLength = length(yVec);
%      rowTemp = indexW(j,:);
%      fprintf(latexFile,'$Y^{%d}\\cdot \\mbox{W}(u_%d, u_%d, u_%d, u_%d)$: \\\\ \n',...
%          [myNumber,rowTemp]);
%      fprintf(latexFile, '%s\n', ' ');
%      
%      for k = 1:myLength
%          yTemp = yVec(k);
%          termTemp = termVec(k);
%          fprintf(latexFile,'term $ %s $: \\\\ \n', latex(yTemp));
%          fprintf(latexFile, '%s\n', ' ');
%          fprintf(latexFile,'$ %s $ \\\\[0.1in] \n', latex(termTemp));
%          fprintf(latexFile, '%s\n', ' ');  
%      end
%      
%      fprintf(latexFile, '\\newpage \n'); 
%      fprintf(latexFile, '%s\n', ' '); 
% end
% fclose(latexFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On W1215
Wy1215 = Y^2*W1215;
symSet1215=symvar(W1215);
w1215=[d2w_MuconjMu, d2w_conjMuconjMu, d2w_conjMuvnormv, d2w_uMu, d2w_uconjMu,...
    d2w_uu, d2w_uvnormv,d3w_uMuconjMu, d3w_uconjMuconjMu, d3w_uuMu,...
    d3w_uuconjMu, d3w_uuu, dw_Mu, dw_conjMu, dw_u, dw_vnormv, w, u, Y];

a1215 = [aV, d2aV_Muvnormv, d2aV_conjMuvnormv, d2aV_conjuMu, d2aV_conjuconjMu,...
    d2aV_conjuvnormv, d2aV_uMu, d2aV_uconjMu, d2aV_uu, d2aV_uvnormv, d2theta_Muvnormv,...
    d3aV_conjuMuconjMu, d3aV_conjuconjMuconjMu, d3aV_uMuMu, d3aV_uMuconjMu,...
    d3aV_uuMu, d3aV_uuconjMu, daV_Mu, daV_conjMu, daV_conju, daV_u, daV_vnormv,...
    theta];

conjW1215 = sym('conjW',[1,length(w1215)]);
for j=1:length(w1215)
    conjW1215(j)=conj(w1215(j));
end
[term1215,y1215] = coeffs(Wy1215,[w1215,conjW1215]);
% length(y1215)=92;
% check on Y^4, Y^6, Y^8, Y^7*conj(u), Y^5*conj(u), Y^4
latexFile = fopen('latex_Wy1215_Nov5.txt','w');
fprintf(latexFile,'$Y^2*W_{1215}$ is: \\\\ \n');
for j=1:length(y1215)
    tempLatex = latex(term1215(j));
    yTemp = y1215(j);
    fprintf(latexFile,'term $ %s $: \\\\ \n', latex(yTemp));
    fprintf(latexFile, '%s\n', ' ');
    fprintf(latexFile,'$ %s $ \\\\[0.1in] \n', tempLatex);
    fprintf(latexFile, '%s\n', ' ');  
end
fclose(latexFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On W1515, W1525
Wy1515 = Y^4*W1515;
Wy1525 = Y^4*W1525;
symSet1515=symvar(W1515);
symSet1525=symvar(W1525);
w1515 = [d2w_MuMu, d2w_Muvnormv, d2w_conjMuconjMu, d2w_uMu, d2w_uconjMu,...
    d2w_uu, d2w_uvnormv, d2w_vnormvvnormv, d3w_uMuMu, d3w_uMuvnormv,...
    d3w_uconjMuconjMu, d3w_uuMu, d3w_uuconjMu, d3w_uuu, d3w_uuvnormv,...
    d4w_uuMuMu, d4w_uuconjMuconjMu, d4w_uuuMu, d4w_uuuu, dw_Mu, dw_conjMu,...
    dw_u, dw_vnormv, w, u, Y];

w1525 = [d2w_MuconjMu, d2w_conjMuvnormv, d2w_uMu, d2w_uconjMu, d2w_uu,...
    d2w_uvnormv, d3w_uMuconjMu, d3w_uconjMuvnormv, d3w_uuMu, d3w_uuconjMu,...
    d3w_uuu, d3w_uuvnormv, d4w_uuMuconjMu, d4w_uuuMu, d4w_uuuconjMu,...
    d4w_uuuu, dw_Mu, dw_conjMu, dw_u, dw_vnormv, w, u, Y];
 

conjW1515 = sym('conjW1515',[1,length(w1515)]);
conjW1525 = sym('conjW1525',[1,length(w1525)]);
for j=1:length(w1515)
    conjW1515(j) = conj(w1515(j));
end

for j=1:length(w1525)
    conjW1525(j) = conj(w1525(j));
end

[term1515,y1515] = coeffs(Wy1515,[w1515,conjW1515]); 
latexFile = fopen('latex_Wy1515_Wy1525_Nov9.txt','w');
fprintf(latexFile,'$Y^4 \\cdot W_{1515}$ is: \\\\ \n');
for j=1:length(y1515)
    tempLatex = latex(term1515(j));
    yTemp = y1515(j);
    fprintf(latexFile,'term $ %s $: \\quad \n', latex(yTemp));
    fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n', tempLatex);
    fprintf(latexFile, '%s\n', ' ');  
end
fprintf(latexFile, '\\newpage \n');
fprintf(latexFile, '%s\n', ' '); 

[term1525,y1525] = coeffs(Wy1525,[w1525,conjW1525]); 
fprintf(latexFile,'$Y^4 \\cdot W_{1525}$ is: \\\\ \n');
for j=1:length(y1525)
    tempLatex = latex(term1525(j));
    yTemp = y1525(j);
    fprintf(latexFile,'term $ %s $: \\quad \n', latex(yTemp));
    fprintf(latexFile,'$\\ds %s $ \\\\[0.1in] \n', tempLatex);
    fprintf(latexFile, '%s\n', ' ');  
end
fclose(latexFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%