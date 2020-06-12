% graph04_nullCone.m
color_pink = (1/255)*[255,204,204];
color_red = (1/255)*[255,0,0];
color_silver = (1/255)*[224,224,224];
color_blue = (1/255)*[0,0,255];
color_purple = [0.75, 0, 0.75];

figure;
hold on
axis equal 
axis off
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);

stepSize = 0.05;
I = [-1:stepSize:1];
[X,Y] = meshgrid(I,I);
size_of_X = size(X);
m = size_of_X(1);
n = size_of_X(2);
Z = zeros(m,n);
Z_at_zero = zeros(m,n);
for ii=1:m
    for j=1:n
        Z(ii,j) = -1/2*(X(ii,j)^2+Y(ii,j)^2);
    end
end
baseM = surf(X,Y,Z);
baseM.FaceColor = 'b';
baseM.FaceAlpha = 0.8;
baseM.EdgeColor = 'none';

tangent_at_zero = surf(X,Y,Z_at_zero);
tangent_at_zero.FaceColor = 'y';
tangent_at_zero.FaceAlpha = 0.3;
tangent_at_zero.EdgeColor = 'none';

theta = [0: stepSize*2*pi: 2*pi];
rVec = [0:stepSize/5:0.5];
scale = 1.2;
[Theta,Radius] = meshgrid(theta, rVec);
size_of_Theta = size(Theta);
m2 = size_of_Theta(1);
n2 = size_of_Theta(2);
X2 = zeros(m2,n2);
Y2 = zeros(m2,n2);
Z2plus = zeros(m2,n2);
Z2minus = zeros(m2,n2);
for ii=1:m2
    for j=1:n2
        X2(ii,j) = Radius(ii,j)*cos(Theta(ii,j));
        Y2(ii,j) = Radius(ii,j)*sin(Theta(ii,j));
        Z2plus(ii,j) = scale*Radius(ii,j);
        Z2minus(ii,j) = -scale*Radius(ii,j);
    end
end
conePlus = surf(X2,Y2,Z2plus);
coneMinus = surf(X2,Y2,Z2minus);
conePlus.FaceColor = 'r';
conePlus.FaceAlpha = 0.3;
conePlus.EdgeColor = 'none';
coneMinus.FaceColor = 'r';
coneMinus.EdgeColor = 'none';
coneMinus.FaceAlpha = 0.1;

dot01 = scatter3(0,0,0,100,'filled');
dot01.MarkerFaceColor='k';

view([45,60,10]);

normalVec = quiver3(0,0,0,0,0,0.5);
normalVec.LineWidth = 3;
normalVec.Color = 'k';
normalVec.MaxHeadSize = 1.5;
string01 = '$\mathbf{n}$';
text01 = text(0,0,0.5,string01,'interpreter','latex');
text01.FontSize = 25;
text01.Color = 'k';

tangentVec = quiver3(0,0,0,0,0.5,0);
tangentVec.LineStyle = '-' ;
tangentVec.LineWidth = 3;
tangentVec.Color = 'k';
tangentVec.MaxHeadSize = 1.5;
string02 = '$\mathbf{e_2}$';
text02 = text(0,0.5,0,string02,'interpreter','latex');
text02.FontSize = 25;
text02.Color = 'k';

tangent2Vec = quiver3(0,0,0,0.5,0,0);
tangent2Vec.LineStyle = '-' ;
tangent2Vec.LineWidth = 3;
tangent2Vec.Color = 'k';
tangent2Vec.MaxHeadSize = 1.5;
string03 = '$\mathbf{e_1}$';
text03 = text(0.55,0,0,string03,'interpreter','latex');
text03.FontSize = 25;
text03.Color = 'k';


