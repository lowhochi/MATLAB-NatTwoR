% variable_graph02Parameter.m
step = 50;
lambda = 0.1;
mu = 0.05;
k = 1;
frequency = -3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_pink = (1/255)*[255,204,204];
color_red = (1/255)*[255,0,0];
color_silver = (1/255)*[224,224,224];
color_blue = (1/255)*[0,0,255];
color_purple = [0.75, 0, 0.75];
color_cyan = [0, 1, 1];
color_orange = [0.8500, 0.3250, 0.0980];
color_green = [0, 1, 0];
color_yellow = [1, 1, 0];
color_black = [0, 0, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xVec = [-2: (4/step): 2];
yVec = [-2: (4/step): 2];
yVecTwo = [-2: (20/step): 2];
thetaVec = [0: 2*pi/step: 2*pi];
rValue_at_zero = k; %y=0
rValue_at_two = k*(1-4*mu); % when y=2 or -2;
rVec = [0: rValue_at_two/step: rValue_at_two];