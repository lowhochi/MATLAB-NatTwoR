% Base manifold: z = 0.25-lambda*y^2, -2<x<2, -2<y<2;
% Bundle: (x^2+(z-2)^2)^(1/2) = k*(1-mu*y^2);
% Let r = |(x,y,z)-(0,y,2)|;
% rValue = k*(1-mu*y^2);
variable_graph02Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
axis equal
% grid on
% axis off
% xlim([-3, 3]);
% ylim([-3, 3]);
% zlim([-1, 4]);
view([60,45,10]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y] = meshgrid(xVec, yVec);
size_of_X = size(X);
m = size_of_X(1);
n = size_of_X(2);
Z = zeros(m,n);
for ii=1:m
    for j=1:n
        Z(ii,j) = 0.25-lambda*Y(ii,j)^2;
    end
end
baseMfd = surf(X,Y,Z);
baseMfd.FaceColor = 'blue';
baseMfd.FaceAlpha = 0.5;
baseMfd.EdgeColor = 'none';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Y2, Theta2] = meshgrid(yVec, thetaVec);
size_of_Y2 = size(Y2);
m2 = size_of_Y2(1);
n2 = size_of_Y2(2);
X2 = zeros(m2,n2);
Z2 = zeros(m2,n2);
for ii=1:m2
    for j=1:n2
        rValue =  k*(1-mu*Y2(ii,j)^2);
        X2(ii,j) = rValue*cos(Theta2(ii,j));
        Z2(ii,j) = 2+ rValue*sin(Theta2(ii,j));
    end
end
bundle = surf(X2,Y2,Z2);
bundle.FaceColor = color_orange;
bundle.FaceAlpha = 0.5;
bundle.EdgeColor = 'none';

[R3, Theta3] = meshgrid(rVec, thetaVec);
size_of_R3 = size(R3);
m3 = size_of_R3(1);
n3 = size_of_R3(2);
X3 = zeros(m3,n3);
Y3plus = zeros(m3,n3);
Y3minus = zeros(m3,n3);
Z3 = zeros(m3,n3);
for ii=1:m3
    for j=1:n3
        Y3plus(ii,j)=2;
        Y3minus(ii,j)=-2;
        X3(ii,j) = R3(ii,j)*cos(Theta3(ii,j));
        Z3(ii,j) = 2 +R3(ii,j)*sin(Theta3(ii,j));
    end
end
bundleCap1 = surf(X3, Y3plus, Z3);
bundleCap2 = surf(X3, Y3minus, Z3);

bundleCap1.FaceColor = color_orange;
bundleCap1.FaceAlpha = 0.5;
bundleCap1.EdgeColor = 'none';
bundleCap2.FaceColor = color_orange;
bundleCap2.FaceAlpha = 0.5;
bundleCap2.EdgeColor = 'none';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bundle curve
x_cap1 = rValue_at_two*cos(thetaVec);
z_cap1 = rValue_at_two*sin(thetaVec)+2;
y_cap1 = zeros(1,length(x_cap1))+2;
y_cap2 = zeros(1,length(x_cap1))-2;
% % % % %
curve_cap1 = plot3(x_cap1, y_cap1, z_cap1);
curve_cap1.LineWidth = 2;
curve_cap1.Color = color_silver;
%curve_cap1.LineStyle = '--';
% % % % %
curve_cap2 = plot3(x_cap1, y_cap2, z_cap1);
curve_cap2.LineWidth = 2;
curve_cap2.Color = color_silver;
curve_cap2.LineStyle = '--';
% get(curve_cap2);

x_latitude1 = k*(1-mu*yVec.^2);
z_latitude1 = zeros(1,length(x_latitude1))+2;
curve_latitude1 = plot3(x_latitude1, yVec, z_latitude1);
curve_latitude2 = plot3(-x_latitude1, yVec, z_latitude1);
curve_latitude1.LineWidth = 2;
curve_latitude1.Color = color_silver;
%curve_latitude1.LineStyle = '--';
curve_latitude2.LineWidth = 2;
curve_latitude2.Color = color_silver;
%curve_latitude2.LineStyle = '--';

z_latitude3 = 2+ k*(1-mu*yVec.^2);
z_latitude4 = 2- k*(1-mu*yVec.^2);
x_latitude3 = zeros(1,length(z_latitude3));
curve_latitude3 = plot3(x_latitude3, yVec, z_latitude3);
curve_latitude4 = plot3(x_latitude3, yVec, z_latitude4);
curve_latitude3.LineWidth = 2;
curve_latitude3.Color = color_silver;
curve_latitude4.LineWidth = 2;
curve_latitude4.Color = color_silver;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z = 0.25-lambda*Y^2;
z_baseMfd = 1/4 - lambda*yVec.^2;
x_baseMfd = zeros(1,length(z_baseMfd));
curve_baseMfd = plot3(x_baseMfd, yVec, z_baseMfd);
curve_baseMfd.LineWidth = 2;
curve_baseMfd.Color = color_blue;

dot01 = scatter3(0,0,1/4,90,'filled');
dot01.MarkerFaceColor= 'r';
dot01.MarkerEdgeColor= 'k';
dot02 = scatter3(0,0,2,30);
dot02.MarkerFaceColor= 'k';
dot02.MarkerEdgeColor= 'none';


x_cap3 = rValue_at_zero*cos(thetaVec);
z_cap3 = rValue_at_zero*sin(thetaVec)+2;
y_cap3 = zeros(1,length(x_cap1));
curve_cap3 = plot3(x_cap3, y_cap3, z_cap3);
curve_cap3.LineWidth = 2;
curve_cap3.Color = color_cyan;
curve_cap3.LineStyle = '--';

u1 = 0.5*rValue_at_zero*cos(3*pi/4);
u2 = 0;
u3 = 0.5*rValue_at_zero*sin(3*pi/4);
uArrow = quiver3(0,0,2,u1,u2,u3);
uArrow.LineWidth = 2.5;
uArrow.LineStyle='-';
uArrow.Color = color_cyan;
uArrow.MaxHeadSize = 2;

dot_uPoint = scatter3(u1,u2,2+u3, 90, 'filled');
dot_uPoint.MarkerFaceColor= color_red;
dot_uPoint.MarkerEdgeColor= 'none';

uArrow2 = quiver3(0,0,1/4,u1,u2,u3);
uArrow2.LineWidth = 2.5;
uArrow2.LineStyle='-';
uArrow2.Color = color_cyan;
uArrow2.MaxHeadSize = 2;

X4 = rValue_at_zero/rValue_at_two*X3;
Z4 = rValue_at_zero/rValue_at_two*(Z3-2) +2;
bundleMiddle = surf(X4, zeros(m3,n3), Z4);
bundleMiddle.FaceColor = 'yellow';
bundleMiddle.FaceAlpha = 0.5;
bundleMiddle.EdgeColor = 'none';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_liftCurve = 1/2*rValue_at_zero*cos(3*pi/4 + frequency*yVec);
z_liftCurve = 2 + 1/2*rValue_at_zero*sin(3*pi/4 + frequency*yVec);
liftCurve = plot3(x_liftCurve, yVec, z_liftCurve);
liftCurve.LineWidth = 2;
liftCurve.Color = color_blue;
liftCurve.LineStyle = '--';
get(liftCurve);
% tangent_on_liftCurv
xTemp = -rValue_at_zero/2*frequency*sin(3*pi/4);
yTemp = 1;
zTemp = rValue_at_zero/2*frequency*cos(3*pi/4);
norm = (xTemp^2+yTemp^2+zTemp^2)^(1/2);
x_tangent = 0.5*xTemp/norm;
y_tangent = 0.5*yTemp/norm;
z_tangent = 0.5*zTemp/norm;
tangent_at_zero = quiver3(u1,u2,u3+2,x_tangent,y_tangent,z_tangent);
tangent_at_zero.LineWidth = 2.5;
tangent_at_zero.LineStyle='-';
tangent_at_zero.Color = color_blue;
tangent_at_zero.MaxHeadSize = 3;
tangent_at_zero.AutoScale = 'off';

% parallel transport of u
x_parallel = zeros(1,length(yVecTwo));
y_parallel = zeros(1,length(yVecTwo));
z_parallel = zeros(1,length(yVecTwo));
for j=1:length(yVecTwo)
    if yVec(j)~=0
        temp = yVecTwo(j)^2/16
        x_parallel(j) = (0.5-temp)*rValue_at_zero*cos(3*pi/4);
        y_parallel(j) = 0;
        z_parallel(j) = (0.5-temp)*rValue_at_zero*sin(3*pi/4);
    end
end
zVecTwo = 0.25-lambda*yVecTwo.^2;
xVecTwo = zeros(1,length(yVecTwo));

parallelTransport = quiver3(xVecTwo, yVecTwo, zVecTwo,...
    x_parallel, y_parallel, z_parallel);
parallelTransport.LineWidth = 2;
parallelTransport.LineStyle='--';
parallelTransport.Color = color_cyan;
parallelTransport.MaxHeadSize = 3;
parallelTransport.AutoScale = 'off';
