% graph_05_geodesicBall_fast.m
% metric: g=(x+y*z)^2*dx^2 +dy^2 +dz^2;
variable_graph5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 40;
stepTwo = 40;
h = 1/step;
phi = [0: (pi/stepTwo): pi];
theta = [0: (2*pi/stepTwo): 2*pi];
[Phi, Theta] = meshgrid(phi,theta);
size_of_Phi = size(Phi);
length_of_Phi = size_of_Phi(1);
width_of_Phi = size_of_Phi(2);
X = zeros(length_of_Phi, width_of_Phi);
Y = zeros(length_of_Phi, width_of_Phi);
Z = zeros(length_of_Phi, width_of_Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u0 = [1.0; 1.0; 1.0];
gValue = [0; 0; 0];
fValue = [0; 0; 0];
GammaTemp = zeros(3,3,3);
for r = 1:length_of_Phi
    for s = 1:width_of_Phi
        disp(['r = ',num2str(r),' ','s = ',num2str(s)]);
        % disp(['r = ',num2str(r)]);
        % disp(['s = ',num2str(s)]);
        v0 = [cos(Theta(r,s))*sin(Phi(r,s));
        sin(Theta(r,s))*sin(Phi(r,s));
        cos(Phi(r,s))];
        gValue = u0;
        fValue = v0;
        for n = 1:step
            gValue_Old = gValue;
            fValue_Old = fValue;
            GammaTemp(1,1,1) = subs(Gamma(1,1,1),[x,y,z],...
                transpose(gValue_Old));
            GammaTemp(1,1,2) = subs(Gamma(1,1,2),[x,y,z],...
                transpose(gValue_Old));
            GammaTemp(1,1,3) = subs(Gamma(1,1,3),[x,y,z],...
                transpose(gValue_Old));
            GammaTemp(1,2,1) = subs(Gamma(1,2,1),[x,y,z],...
                transpose(gValue_Old));
            GammaTemp(2,1,1) = GammaTemp(1,2,1);
            GammaTemp(1,3,1) = subs(Gamma(1,3,1),[x,y,z],...
                transpose(gValue_Old));
            GammaTemp(3,1,1) = GammaTemp(1,3,1);
            % % % % %            
            gValue = gValue_Old +h*fValue_Old;
            fPart1 = 0;
            fPart2 = 0;
            fPart3 = 0;
            for ii=1:3
                for j=1:3
                    fPart1 = fPart1 +fValue_Old(ii)*fValue_Old(j) ...
                        *GammaTemp(ii,j,1);
                    fPart2 = fPart2 +fValue_Old(ii)*fValue_Old(j)...
                        *GammaTemp(ii,j,2);
                    fPart3 = fPart3 +fValue_Old(ii)*fValue_Old(j)...
                        *GammaTemp(ii,j,3);
                end
            end
            fValue = fValue_Old -h*[fPart1; fPart2; fPart3];
        end
        X(r,s) = gValue(1);
        Y(r,s) = gValue(2);
        Z(r,s) = gValue(3);       
    end
end
clearvars fValue gValue fPart1 fPart2 fPart3 GammaTemp r s n 
clearvars ii j fValue_Old gValue_Old h v0
save('Data_graph05_step40.mat');
% myGeodBall = surf(X,Y,Z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
load('Data_graph05_step40.mat');
colorMat = zeros(41, 41, 3);
for j=1:41
    for k=1:41
        colorMat(j,k,1) = 0;
        colorMat(j,k,2) = 1;
        colorMat(j,k,3) = 0;
    end
end
figure
hold on
axis equal
grid on
xlim([-0.2,2.2]);
ylim([-0.2,2.2]);
zlim([-0.2,2.2]);
xlabel('x');
ylabel('y');
zlabel('z');
myTitle = {'Metric $\mathbf{g = (x+yz)^2\,dx^2 + dy^2 + dz^2}$';
    'Geodesic-1-ball centered at (1,1,1)'};
title(myTitle,'Interpreter','latex');
myGeodBall = surf(X,Y,Z,colorMat);
%myGeodBall.FaceColor = 'y';
myGeodBall.FaceAlpha = 0.25;
myGeodBall.EdgeColor = 'w';
ptOne =  scatter3(1,1,1,100,'filled');
ptOne.MarkerFaceColor='k';
textOne = text(1.02,1.02,1.02,'$\mathbf{(1,1,1)}$','Interpreter','latex');
textOne.FontSize=14;
view([30,30,10]);

Xs = zeros(length_of_Phi, width_of_Phi);
Ys = zeros(length_of_Phi, width_of_Phi);
Zs = zeros(length_of_Phi, width_of_Phi);
for j=1:length_of_Phi
    for k=1:width_of_Phi
        Xs(j,k) = 1 + sin(Phi(j,k))*cos(Theta(j,k));
        Ys(j,k) = 1 + sin(Phi(j,k))*sin(Theta(j,k));
        Zs(j,k) = 1 + cos(Phi(j,k));   
    end
end
mySphere = surf(Xs,Ys,Zs);
mySphere.FaceColor = 'none';
mySphere.EdgeColor = 'b';
mySphere.LineStyle = ':';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





