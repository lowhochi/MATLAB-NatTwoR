% example_flatModel_computeWeyl.m
load('DataWeyl2_NatTwoR_CR_rmW.mat');
clearvars myChar myCharR list_of_variables_Main4
clearvars list_of_variables_Main2 list_of_variables_Main3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST (OKAY)
% conformally flat metric g = 1/(1+x^2)*(dx^2+dy^2+dz^2);
% u1 = (x^2+1)^(1/2)*d/dx;
% u2 = d/dy - x*d/dz;
% u3 = d/dz + x*d/dy;
% f = i/(2*(1+x^2)^(1/2))*(u^2-1)^2;
% aVs = -v1normv*x/(1+x^2)^(1/2) -i/2*(1+v1normv^2)/(1+x^2)^(1/2);
% thetas = -1/(1+x^2)^(1/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x y z real
muVec = [mu1; mu2; mu3];
conjMuVec = [conj(mu1); conj(mu2); conj(mu3)];
vnormvVec = [v1normv; v2normv; v3normv];
prompt11 = 'Input metric g; dx^2: \n';
prompt12 = 'Input metric g; dxdy: \n';
prompt13 = 'Input metric g; dxdz: \n';
prompt22 = 'Input metric g; dy^2: \n';
prompt23 = 'Input metric g; dydz: \n';
prompt33 = 'Input metric g; dz^2: \n';
promptE1 = 'Input E1; [d/dx; d/dy; d/dz]: \n';
promptE2 = 'Input E2; [d/dx; d/dy; d/dz]: \n';
promptE3 = 'Input E3; [d/dx; d/dy; d/dz]: \n';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric information
% g11 = input(prompt11);
% g12 = input(prompt12);
% g13 = input(prompt13);
% g22 = input(prompt22);
% g23 = input(prompt23);
% g33 = input(prompt33);
% gMat = [g11, g12, g13;
%     g12, g22, g23;
%     g13, g23, g33];
% e1 = input(promptE1);
% e2 = input(promptE2);
% e3 = input(promptE3);
% eMat = [e1, e2, e3];
gMat = [(x+y*z)^2, 0, 0;
    0, 1, 0;
    0, 0, 1];
eMat = [1/(x+y*z), 0, 0;
    0, 1, 0;
    0, 0, 1];
checkMat = transpose(eMat)*gMat*eMat;
for j=1:3
    for k=1:3
        checkMat(j,k) = simplify(checkMat(j,k));
    end
end
G = myChristoffel(gMat, eMat, [x,y,z]);
% RmThree = myRiemThreeMfd(eMat, G, [x,y,z]);
% cotton = myCottonTensor(eMat, G, RmThree, [x, y, z]);
% define theta
thetas = G(1,2,3) + G(2,3,1) + G(3,1,2);
aVs = 0;
f = 0;
muNormSquare = mu1*conj(mu1)+mu2*conj(mu2)+mu3*conj(mu3);
for m=1:3
    aVs = aVs -i*conj(muVec(m))/muNormSquare...
        *(G(m,1,2)*muVec(3)+G(m,2,3)*muVec(1)+G(m,3,1)*muVec(2));
    f = f -i/2*muVec(m)...
        *(G(m,1,2)*muVec(3)+G(m,2,3)*muVec(1)+G(m,3,1)*muVec(2));
end
aVs = aVs + i*thetas;
aVs = complex_simple3(aVs, u);
f = complex_simple3(f, u);
variable_example_computeWeyl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fefferman metric
% [T4, aM, aV, dT4_u, dw_u, h11, rho, u, w]
aMs = subs(aMSub, daV_u, daVs_u);
T4s = subs(T4Sub, daV_conju, daVs_conju);
h11s = subs(h11Sub, d2aV_uconju, d2aV_uconjuR);
h11s = subs(h11s, [aV, theta], [aVs, thetas]);
dT4s_u = subs(dT4_uSub, d2aV_uconju, d2aV_uconjuR);
dT4s_u = subs(dT4s_u, [aV, theta, daV_conju], ...
    [aVs, thetas, daVs_conju]);
variableSetF = [aM, aV, T4, dT4_u, rho, h11, w, dw_u];
subSetF = [aMs, aVs, T4s, dT4s_u, 0, h11s, f, df_u];

Fs = sym('Fs',[6 6]);
for j=1:6
    for k=1:6
        temp = F0(j,k);
        temp = subs(temp, variableSetF, subSetF);
        Fs(j,k) = complex_simple3(temp, u);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars j k m g11 g12 g13 g21 g22 temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WeylTensor = sym('WeylTensor',[120,5]);
% for j=1:120
%     for k=1:4
%         WeylTensor(j,k) = WeylTwo(j,k);
%     end
%     temp = WeylTwo(j,5); 
%     % % % %
%     for m=1:length(variableSetf4)
%         mychar = char(variableSetf4(m));
%         mycharf = strrep(mychar,'w','f');
%         eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
%     end
%     % % % %
%     for m=1:length(variableSetf3)
%         mychar = char(variableSetf3(m));
%         mycharf = strrep(mychar,'w','f');
%         eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
%     end
%     % % % %
%     for m=1:length(variableSetf2)
%         mychar = char(variableSetf2(m));
%         mycharf = strrep(mychar,'w','f');
%         eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
%     end
%     % % % %
%     for m=1:length(variableSetf1)
%         mychar = char(variableSetf1(m));
%         mycharf = strrep(mychar,'w','f');
%         eval(['temp=subs(temp,', mychar, ',', mycharf, ');']);
%     end
%     % % % %
%     for m=1:length(variableSetT)
%         mychar = char(variableSetT(m));
%         mycharT = insertAfter(mychar,'a','s');
%         eval(['temp=subs(temp,', mychar, ',', mycharT, ');']);
%     end
%     % % % %
%     for m=1:length(variableSet4)
%         mychar = char(variableSet4(m));
%         mychars = insertAfter(mychar,'V','s');
%         eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
%     end
%     % % % %
%     for m=1:length(variableSet3)
%         mychar = char(variableSet3(m));
%         mychars = insertAfter(mychar,'V','s');
%         eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
%     end
%     % % % %
%     for m=1:length(variableSet2)
%         mychar = char(variableSet2(m));
%         mychars = insertAfter(mychar,'V','s');
%         eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
%     end
%     % % % %
%     
%     for m=1:length(variableSet1)
%         mychar = char(variableSet1(m));
%         mychars = insertAfter(mychar,'V','s');
%         eval(['temp=subs(temp,', mychar, ',', mychars, ');']);
%     end    
%     WeylTensor(j,5) = complex_simple3(temp, u);
% end
% clearvars mychar mychars temp j k


