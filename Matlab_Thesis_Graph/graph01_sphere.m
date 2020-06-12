step = 100;
color_pink = (1/255)*[255,204,204];
color_red = (1/255)*[255,0,0];
color_silver = (1/255)*[224,224,224];
color_blue = (1/255)*[0,0,255];
color_purple = [0.75, 0, 0.75];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = [0: (pi/step): pi];
theta = [0: (2*pi/step): 2*pi];
[Phi, Theta] = meshgrid(phi, theta);
size_of_Theta = size(Theta);
m = size_of_Theta(1);
n = size_of_Theta(2);
X = zeros(m,n);
Y = zeros(m,n);
Z = zeros(m,n);
for j=1:m
    for k=1:m
        X(j,k) = sin(Phi(j,k))*cos(Theta(j,k));
        Y(j,k) = sin(Phi(j,k))*sin(Theta(j,k));
        Z(j,k) = cos(Phi(j,k));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myGraph = figure;
mySphere = surf(X,Y,Z);
mySphere.FaceColor = 'cyan';
mySphere.FaceAlpha = 0.2;
mySphere.EdgeColor = 'blue';
mySphere.EdgeAlpha = 0.1;
hold on
% grid off
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
axis equal
axis off
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
zlim([-1.5, 1.5]);
view([45,45,15]);

line01 = quiver3(0,0,0,1.5,0,0);
line01.LineWidth = 3;
line01.Color = 'r';
line01.MaxHeadSize = 0.5;
string01 = '{\bf Re(\mu)}';
text01 = text(1.5,0,0,string01);
text01.FontSize = 16;
text01.Color = 'k';

line02 = quiver3(0,0,0,0,1.5,0);
line02.LineWidth = 3;
line02.Color = 'r';
line02.MaxHeadSize = 0.5;
string02 = '{\bf Im(\mu)}';
text02 = text(0,1.5,0,string02);
text02.FontSize = 16;
text02.Color = 'k';

line03 = quiver3(0,0,0,0,0,1);
line03.LineWidth = 3;
line03.Color = 'k';
line03.MaxHeadSize = 0.5;
string03 = '$$\mathbf{\frac{v}{|v|}}$$';
text03 = text(0,0,1.2,string03,'interpreter','latex');
% text03.interpreter = 'latex';
text03.FontSize = 16;
text03.Color = 'k';

dot01 = scatter3(0,0,0,100,'filled');
dot01.MarkerFaceColor='b';
string04 = '$\mathbf{T_xM}$';
text04 = text(0,0.1,0.1,string04,'interpreter','latex');
text04.FontSize=16;

xInt = [-1.2: 0.1: 1.2];
[Xplane,Yplane] = meshgrid(xInt,xInt);
size_of_Xplane = size(Xplane);
length_of_Xplane = size_of_Xplane(1);
width_of_Xplane = size_of_Xplane(2);
Zplane = zeros(length_of_Xplane, width_of_Xplane);
myPlane = surf(Xplane,Yplane,Zplane);
myPlane.FaceColor = color_pink;
%myPlane.FaceAlpha = 0.7;
myPlane.EdgeColor = color_red;
myPlane.EdgeAlpha = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%