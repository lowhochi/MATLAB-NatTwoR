% graph06_pieChart.m
color_pink = (1/255)*[255,204,204];
color_red = (1/255)*[255,0,0];
color_blue = (1/255)*[0,0,255];
color_purple = [0.75, 0, 0.75];
color_green = [0, 1, 0];
color_yellow = [1, 1, 0];
color_black = [0, 0, 0];
color_silver = (1/255)*[224,224,224];
color_orange = [0.8500, 0.3250, 0.0980];
color_cyan = [0, 1, 1];
color_grey = [0.3, 0.3, 0.3];
color_deepRed = [0.6350, 0.0780, 0.1840];

myColor={color_blue, color_cyan, color_yellow,...
    color_orange, color_deepRed, color_grey};

myFigure = figure;
set(myFigure, 'defaulttextinterpreter','latex');
myData = [65, 21, 8, 11, 10, 5];
labels = {'Coeffs with at least one 6', ...
    'Coeffs on $T^{1,0}\oplus T^{0,1}$', ...
    'Other $\mathcal{W}_{13kl}$ or $\mathcal{W}_{24kl}$', ...
    'Coeffs related to $\rho$', ...
    'Almost linear in $w$', ...
    'The hardest $\mathcal{W}_{ijkl}$'};
myPie = pie(myData);
myLeg = legend(labels,'interpreter','latex');
myLeg.FontSize = 20;
myLeg.Position = [0.75 0.4 0.2 0.25];
% textSize
for j=1:6
     myText = myPie(2*j);
     myText.FontSize = 25;
     myText.FontName = 'Arial';
     myText.String = num2str(myData(j));
     
     mySec = myPie(2*j-1);
     mySec.FaceColor = myColor{j};
end
