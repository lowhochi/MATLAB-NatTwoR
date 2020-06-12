% Twistor_CR_funWv2_main1.m
syms x y z real
syms u w
syms dw_x dw_y dw_z dw_u
derivativeDict.x = [1; 0; 0; 0; 0];
derivativeDict.y = [0; 1; 0; 0; 0];
derivativeDict.z = [0; 0; 1; 0; 0];
derivativeDict.u = [0; 0; 0; 1; 0];
derivativeDict.w = [dw_x; dw_y; dw_z; dw_u; 0];
CVar1 = [x,y,z,u,w];

mu1 =u^2-1; 
mu2= 2*u; 
mu3= i*(u^2+1);
v1 = i*(mu2*conj(mu3)-mu3*conj(mu2));
v2 = i*(mu3*conj(mu1)-mu1*conj(mu3));
v3 = i*(mu1*conj(mu2)-mu2*conj(mu1));
v1 = complex_simple3(v1,[u]);
v2 = complex_simple3(v2,[u]);
v3 = complex_simple3(v3,[u]);
norm_of_v = sqrt(v1*v1+v2*v2+v3*v3);
norm_of_v = complex_simple3(norm_of_v,[u]);
v1normv = v1/norm_of_v; %T1
v2normv = v2/norm_of_v; %T2
v3normv = v3/norm_of_v; %T3
v1normv = complex_simple3(v1normv, u);
v2normv = complex_simple3(v2normv, u);
v3normv = complex_simple3(v3normv, u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR distribution in [d/dx, d/dy, d/dz, d/du, d/dconj(u)]
% conjX1 = mu1*d/dx +mu2*d/dy +mu3*d/dx + w*d/du
conjX1 = [mu1; mu2; mu3; w; 0];
X1 = [conj(mu1); conj(mu2); conj(mu3); 0 ;conj(w)];
X2 = [0; 0; 0; 1; 0];
conjX2 = [0; 0; 0; 0; 1];
T = [v1normv; v2normv; v3normv; 0; 0];

% Lie bracket of X1, conj_X1, X2, conj_X2 and T. 
lie_X1_conjX1 = lieBracket_CR_funWv2_main1(X1,conjX1,CVar1,derivativeDict);
lie_X1_X2 = lieBracket_CR_funWv2_main1(X1,X2,CVar1,derivativeDict);
lie_X1_conjX2 = lieBracket_CR_funWv2_main1(X1,conjX2,CVar1,derivativeDict);
lie_X2_conjX1 = lieBracket_CR_funWv2_main1(X2,conjX1,CVar1,derivativeDict);
lie_X2_conjX2 = lieBracket_CR_funWv2_main1(X2,conjX2,CVar1,derivativeDict);
lie_conjX1_conjX2 = lieBracket_CR_funWv2_main1(conjX1,conjX2,CVar1,derivativeDict);
    
lie_X1_T = lieBracket_CR_funWv2_main1(X1,T,CVar1,derivativeDict);
lie_conjX1_T = lieBracket_CR_funWv2_main1(conjX1,T,CVar1,derivativeDict);
lie_X2_T = lieBracket_CR_funWv2_main1(X2,T,CVar1,derivativeDict);
lie_conjX2_T = lieBracket_CR_funWv2_main1(conjX2,T,CVar1,derivativeDict);

lieMain1 = cell(5,5);
lieMain1{1,1} = zeros(5,1); %[X1,X1]
lieMain1{1,2} = lie_X1_conjX1;
lieMain1{1,3} = lie_X1_X2;
lieMain1{1,4} = lie_X1_conjX2;
lieMain1{1,5} = lie_X1_T;
lieMain1{2,1} = -lie_X1_conjX1; %[conjX1,X1]
lieMain1{2,2} = zeros(5,1); %[conjX1,conjX1]
lieMain1{2,3} = -lie_X2_conjX1; %[conjX1,X2]
lieMain1{2,4} = lie_conjX1_conjX2;
lieMain1{2,5} = lie_conjX1_T;
lieMain1{3,1} = -lie_X1_X2; %[X2,X1]
lieMain1{3,2} = lie_X2_conjX1;
lieMain1{3,3} = zeros(5,1); %[X2,X2]
lieMain1{3,4} = lie_X2_conjX2;
lieMain1{3,5} = lie_X2_T;
lieMain1{4,1} = -lie_X1_conjX2;
lieMain1{4,2} = -lie_conjX1_conjX2;
lieMain1{4,3} = -lie_X2_conjX2;
lieMain1{4,4} = zeros(5,1); %[conjX2,conjX2]
lieMain1{4,5} = lie_conjX2_T;
lieMain1{5,1} = -lie_X1_T;
lieMain1{5,2} = -lie_conjX1_T;
lieMain1{5,3} = -lie_X2_T;
lieMain1{5,4} = -lie_conjX2_T;
lieMain1{5,5} = zeros(5,1); %[T,T]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha0 = [v1normv, v2normv, v3normv, 0, 0];

dalpha0 = sym('dalpha0',[5 5]);
for ii=1:5
    dalpha_ii =  df_main1_CR_funWv2(alpha0(ii),CVar1,derivativeDict);
    for jj=1:5
        dalpha_jj = df_main1_CR_funWv2(alpha0(jj),CVar1,derivativeDict);
        temp = 1/2*(dalpha_jj(ii)- dalpha_ii(jj));
        dalpha0(ii,jj) = complex_simple3(temp, u);
    end
end
% The Levi form L
L = -i*dalpha0;
g_X1_conjX1 = transpose(X1)*L*conjX1;
g_X1_conjX2 = transpose(X1)*L*conjX2;
g_X2_conjX1 = transpose(X2)*L*conjX1;
g_X2_conjX2 = transpose(X2)*L*conjX2;
g_X1_conjX1 = complex_simple3(g_X1_conjX1,[u,w,dw_x,dw_y,dw_z,dw_u]);
g_X1_conjX2 = complex_simple3(g_X1_conjX2,[u,w,dw_x,dw_y,dw_z,dw_u]);
g_X2_conjX1 = complex_simple3(g_X2_conjX1,[u,w,dw_x,dw_y,dw_z,dw_u]);
g_X2_conjX2 = complex_simple3(g_X2_conjX2,[u,w,dw_x,dw_y,dw_z,dw_u]);

h11 = g_X1_conjX1;
h12 = g_X1_conjX2;
h21 = g_X2_conjX1;
h22 = g_X2_conjX2;
h = [h11, h12; h21, h22];
% The Webster metric in (X1,conjX1,X2,conjX2,T): g0
g0 = [0, g_X1_conjX1, 0, g_X1_conjX2, 0;
    g_X1_conjX1, 0, g_X2_conjX1, 0, 0;
    0, g_X2_conjX1, 0, g_X2_conjX2, 0;
    g_X1_conjX2, 0, g_X2_conjX2, 0, 0;
    0, 0, 0, 0, 1];
% The Wester metric in (x,y,z,u,conj(u)): g
CRVector =[X1, conjX1, X2, conjX2, T];
CRVectorInv = inv(CRVector);
g = transpose(CRVectorInv)*g0*CRVectorInv;

clearvars h11 h12 h21 h22  L ii jj temp
clearvars dalpha_ii dalpha_jj g_X1_conjX1 g_X1_conjX2 
clearvars g_X2_conjX1 g_X2_conjX2
clearvars lie_X1_conjX1 lie_X1_X2 lie_X1_conjX2 lie_X2_conjX1 lie_X2_conjX2
clearvars lie_conjX1_conjX2 lie_X1_T lie_conjX1_T lie_X2_T lie_conjX2_T
save('DataMain1_CR_funWv2.mat');