% graph03_rationalpar.m
color_pink = (1/255)*[255,204,204];
color_red = (1/255)*[255,0,0];
color_blue = (1/255)*[0,0,255];
color_purple = [0.75, 0, 0.75];
color_green = [0, 1, 0];
color_yellow = [1, 1, 0];
color_black = [0, 0, 0];

step = 50;
theta = [0: 2*pi/step: 2*pi];
N = length(theta);
uReVec = zeros(1,N);
uImVec = zeros(1,N);
radius = [1:1:5];
M = length(radius);
radiusColor = {color_red, color_purple, color_green,...
    color_blue, color_black};

plotOne = subplot(1,2,1);
hold on
set(plotOne, 'DataAspectRatio',[1 1 1]);
title(plotOne,'u-Plane');
plotOne.FontSize = 14;
plotOne.XLim = [-5.2,5.2];
plotOne.YLim = [-5.2,5.2];
plotOne.XLabel.String = 'Re(u)';
plotOne.YLabel.String = 'Im(u)';
plotOne.XAxisLocation = 'origin';
plotOne.YAxisLocation = 'origin';
% plotOne.XTickLabel = [];
% plotOne.YTickLabel = [];
for j=1:M
    R = radius(j);
    for k=1:N
        uReVec(k) = R*cos(theta(k));
        uImVec(k) = R*sin(theta(k));
    end
    uCircle = plot(uReVec,uImVec,'-');
    uCircle.LineWidth = 2;
    uCircle.Color = radiusColor{j};
end
uZero =  scatter(plotOne,0,0,100,'filled');
uZero.MarkerFaceColor='r';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotTwo = subplot(1,2,2);
hold on
set(plotTwo, 'DataAspectRatio',[1 1 1]);
plotTwo.XLim = [-1.2,1.2];
plotTwo.YLim = [-1.2,1.2];
plotTwo.ZLim = [-1.2,1.2];
view(plotTwo,[45,30,10]);
set(plotTwo,'XTickLabel',[]);
set(plotTwo,'YTickLabel',[]);
set(plotTwo,'ZTickLabel',[]);
myTitle={'$$\mathbf{ \frac{v}{|v|} = \frac{i\,\mu\times\overline{\mu}}{|\mu|^2} }$$';
    '$$\mathbf{ \mu = \big(u^2-1,\, 2u,\, i(u^2+1)\big) }$$'};

title(plotTwo, myTitle,...
    'Interpreter', 'latex');
plotTwo.FontSize = 14;
phi = [0: pi/step: pi];
[Theta,Phi] = meshgrid(theta,phi);
size_of_Theta = size(Theta);
length_of_Theta = size_of_Theta(1);
width_of_Theta = size_of_Theta(2);
X2 = zeros(length_of_Theta, width_of_Theta);
Y2 = zeros(length_of_Theta, width_of_Theta);
Z2 = zeros(length_of_Theta, width_of_Theta);

for j=1:length_of_Theta
    for k=1:width_of_Theta
        X2(j,k) = sin(Phi(j,k))*cos(Theta(j,k));
        Y2(j,k) = sin(Phi(j,k))*sin(Theta(j,k));
        Z2(j,k) = cos(Phi(j,k));
    end
end

mySphere = surf(plotTwo, X2, Y2, Z2);
mySphere.FaceColor = 'none';
mySphere.FaceAlpha = 0.3;
mySphere.EdgeColor = [0.2,1,1];

% radius = [0:0.5:5];
% M = length(radius);
xCurve = zeros(1,N);
yCurve = zeros(1,N);
zCurve = zeros(1,N);
for j=1:M
    R = radius(j);
    for k=1:N
       u0 = R*cos(theta(k))+1i*R*sin(theta(k));
       muVec = [u0^2-1, 2*u0, 1i*(u0^2+1)];
       conjMuVec = [conj(u0)^2-1, 2*conj(u0), -1i*(conj(u0)^2+1)];
       norm_sq_of_mu = muVec(1)*conj(muVec(1))+ muVec(2)*conj(muVec(2))...
           + muVec(3)*conj(muVec(3));
       vnormv = 1i*cross(muVec,conjMuVec)/norm_sq_of_mu;
       xCurve(k) = vnormv(1);
       yCurve(k) = vnormv(2);
       zCurve(k) = vnormv(3);
    end
    vCurve = plot3(plotTwo,xCurve,yCurve,zCurve);
    vCurve.LineWidth = 2;
    vCurve.Color = radiusColor{j};
end
uZeroSphere = scatter3(plotTwo,0,1,0,100,'filled');
uZeroSphere.MarkerFaceColor='r';
myString ='$\mathbf{e_2} = (0,1,0)$';
textuZero = text(plotTwo,0,1.1,0,myString,'Interpreter','latex');
textuZero.FontSize=14;
textuZero.Color='r';

