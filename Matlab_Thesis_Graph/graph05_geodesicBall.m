% graph05_geodesicBall.m
% metric: g=(x+y*z)^2*dx^2 +dy^2 +dz^2;
variable_graph5
figure;
hold on
axis equal
% axis tight
% xlim([0,2]);
% ylim([-1,1]);
% zlim([-1,1]);
step = 10;
stepTwo = 10;
h = 1/step;
% % % % % IC
phi = [0: (pi/stepTwo): pi];
theta = [0: (2*pi/stepTwo): 2*pi];
[Phi, Theta] = meshgrid(phi,theta);
size_of_Phi = size(Phi);
length_of_Phi = size_of_Phi(1);
width_of_Phi = size_of_Phi(2);
X = zeros(length_of_Phi, width_of_Phi);
Y = zeros(length_of_Phi, width_of_Phi);
Z = zeros(length_of_Phi, width_of_Phi);
% % % % %
u0 = [1; 1; 1];
for r = 1:length_of_Phi
    for s = 1:width_of_Phi
        gMat = zeros(3, step+1);
        fMat = zeros(3, step+1);
        v0 = [cos(Theta(r,s))*sin(Phi(r,s));
            sin(Theta(r,s))*sin(Phi(r,s));
            cos(Phi(r,s))];
        gMat(:,1) = u0;
        fMat(:,1) = v0;
        % Euler Method
        for n =1: step
            gMat(1,n+1) = gMat(1,n)+h*fMat(1,n);
            gMat(2,n+1) = gMat(2,n)+h*fMat(2,n);
            gMat(3,n+1) = gMat(3,n)+h*fMat(3,n);
    
            fPart1 = 0;
            fPart2 = 0;
            fPart3 = 0;
            for ii=1:3
                for j=1:3
                    Gij_1 = subs(Gamma(ii,j,1),[x,y,z],...
                        [gMat(1,n),gMat(2,n),gMat(3,n)]);
                    Gij_2 = subs(Gamma(ii,j,2),[x,y,z],...
                        [gMat(1,n),gMat(2,n),gMat(3,n)]);
                    Gij_3 = subs(Gamma(ii,j,3),[x,y,z],...
                        [gMat(1,n),gMat(2,n),gMat(3,n)]);
                    fPart1 = fPart1 +fMat(ii,n)*fMat(j,n)*Gij_1;
                    fPart2 = fPart2 +fMat(ii,n)*fMat(j,n)*Gij_2;
                    fPart3 = fPart3 +fMat(ii,n)*fMat(j,n)*Gij_3;
                end
            end
            fMat(1,n+1) = fMat(1,n)-h*fPart1;
            fMat(2,n+1) = fMat(2,n)-h*fPart2;
            fMat(3,n+1) = fMat(3,n)-h*fPart3;
        end
        X(r,s) = gMat(1,step+1);
        Y(r,s) = gMat(2,step+1);
        Z(r,s) = gMat(3,step+1);
        clearvars gMat fMat
    end
end

myGeodBall = surf(X,Y,Z);

%myGeod1 = plot3(gMat(1,:),gMat(2,:),gMat(3,:));
%myGeod1.LineWidth = 2;
