% graph07_CRmfld.m
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

myFig = figure;
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
xlim([-1.1,1.1]);
ylim([-1.1,1.1]);
zlim([-1.1,1.1]);
axis equal
grid off
axis off
view([45,45,25]);
hold on
step = 20;
xVec = [-1: 2/step: 1];
[X,Y] = meshgrid(xVec,xVec);
size_of_X = size(X);
length_of_X = size_of_X(1);
width_of_X = size_of_X(2);
Z = zeros(length_of_X, width_of_X);
for ii=1:length_of_X
    for j=1:width_of_X
        Z(ii,j) = -0.25*X(ii,j)^2 - 0.25*Y(ii,j)^2;
    end
end
myBase = surf(X,Y,Z);
myBase.FaceColor = color_blue;
myBase.EdgeColor = 'white';
myBase.FaceAlpha = 0.5;

dot01 = scatter3(0,0,0,120,'filled');
dot01.MarkerFaceColor= color_black;
dot01.MarkerEdgeColor= 'none';
x2Arrow = quiver3(0,0,0.5,0.5,0,0);
y2Arrow = quiver3(0,0,0.5,0,0.5,0);
z2Arrow = quiver3(0,0,0.5,0,0,0.5);

x2Arrow.LineWidth = 3;
x2Arrow.LineStyle='-';
x2Arrow.Color = color_black;
x2Arrow.MaxHeadSize = 2;
y2Arrow.LineWidth = 3;
y2Arrow.LineStyle='-';
y2Arrow.Color = color_black;
y2Arrow.MaxHeadSize = 2;
z2Arrow.LineWidth = 3;
z2Arrow.LineStyle='-';
z2Arrow.Color = color_black;
z2Arrow.MaxHeadSize = 2;

xLine = zeros(1,6);
yLine = zeros(1,6);
zLine = [0:0.1:0.5];
myDashLine = plot3(xLine, yLine, zLine);
myDashLine.LineWidth = 1.5;
myDashLine.LineStyle = '--';
myDashLine.Color = color_black;
myText = text(0,0.2,1,'\bf CTN');
myText.FontSize = 20;

xPvec = [-0.1:0.05:0.4];
[Xplane,Yplane] = meshgrid(xPvec, xPvec);
size_of_Xplane = size(Xplane);
length_of_Xplane = size_of_Xplane(1);
width_of_Xplane = size_of_Xplane(2);
Zplane = zeros(length_of_Xplane, width_of_Xplane);
nVec = cross([0,1,0],[-1/4,0,1/4]);
for ii=1:length_of_Xplane
    for j=1:width_of_Xplane
        Zplane(ii,j) = 0.75- nVec(1)/nVec(3)*Xplane(ii,j)...
            -nVec(2)/nVec(3)*Yplane(ii,j);
    end
end
myCRplane = surf(Xplane,Yplane,Zplane);
myCRplane.FaceColor = color_red;
myCRplane.FaceAlpha = 0.1;
myCRplane.EdgeColor = color_deepRed;
myText02 = text(1.2,0,1,'\bf CR-Plane');
myText02.FontSize = 20;
myText02.Color = color_red;
myText02.Interpreter = 'latex';
myText02.Position = [1.2, 0, 1];