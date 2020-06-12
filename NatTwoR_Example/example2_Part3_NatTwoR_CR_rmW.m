% example2_part3_NatTwoR_CR_rmW.m

% load('DataWeyl2_NatTwoR_CR_rmW.mat');
% variable_eg2_NatTwoR_CR_rmW
load('DataExample2_Dec22_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W1212 = WeyltTwo(1,5);
% W1215 = WeyltTwo(4,5);
% W1515 = WeyltTwo(43,5);
% W1525 = WeyltTwo(47,5);

factor1212 = factor(W1212);
factor1215 = factor(W1215);
factor1525 = factor(W1525);

save('DataExample2_Dec22_NatTwoR_CR_rmW.mat', 'factor1212',...
    'factor1215', 'factor1525', '-append');

