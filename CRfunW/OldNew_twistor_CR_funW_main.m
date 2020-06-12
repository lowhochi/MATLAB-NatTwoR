syms x y z real
syms u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w = w(x,y,z,u) holo in u. 
syms w
syms dw_x dw_y dw_z dw_u 
syms d2w_xx d2w_xy d2w_xz d2w_yy d2w_yz d2w_zz
syms d2w_xu d2w_yu d2w_zu d2w_uu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable sets.
M = [x,y,z,u];
dw = [dw_x, dw_y, dw_z, dw_u];
Mdw = [u,w,dw];
d2w = [d2w_xx, d2w_xy, d2w_xz, d2w_xu;
    d2w_xy, d2w_yy, d2w_yz, d2w_yu;
    d2w_xz, d2w_yz, d2w_zz, d2w_zu;
    d2w_xu, d2w_yu, d2w_zu, d2w_uu];
d2wRow=[d2w_xx, d2w_xy, d2w_xz, d2w_yy, d2w_yz, d2w_zz, ...
    d2w_xu, d2w_yu, d2w_zu, d2w_uu];
Md2w = [Mdw,d2wRow];
% Md2w = [u, w, dw_x, dw_y, dw_z, dw_u, ...
%   d2w_xx, d2w_xy, d2w_xz, d2w_yy, d2w_yz, d2w_zz, ...
%   d2w_xu, d2w_yu, d2w_zu, d2w_uu]
RVar2 = [x,y,z];
CVar2 = [u, w, dw_x, dw_y, dw_z, dw_u];
dRVar2 = eye(3);
dCVar2 = [0, dw_x, d2w_xx, d2w_xy, d2w_xz, d2w_xu;
    0, dw_y, d2w_xy, d2w_yy, d2w_yz, d2w_yu;
    0, dw_z, d2w_xz, d2w_yz, d2w_zz, d2w_zu;
    1, dw_u, d2w_xu, d2w_yu, d2w_zu, d2w_uu];
MVar2 = Md2w;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR distribution: conj_X1, X1, X2, conj_X2, T. 
conj_X1 = [u^2-1; 2*u; i*(u^2+1); w; 0];
X1 = [conj(u)^2-1; 2*conj(u); -i*(conj(u)^2+1); 0; conj(w)];
X2 = [0; 0; 0; 1; 0];
conj_X2 = [0; 0; 0; 0; 1];
T = [(u+conj(u))/(1+u*conj(u)); (1-u*conj(u))/(1+u*conj(u)); ...
    i*(u-conj(u))/(1+u*conj(u)); 0; 0];
holo_vector = [X1, X2];
antiholo_vector = [conj_X1, conj_X2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lie bracket of X1, conj_X1, X2, conj_X2 and T. 
lie_funW = lie_twistor_CR_funW(T,u,w,dw,Mdw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudo-hermitian form: alpha0.
alpha0 = [(u+conj(u))/(1+u*conj(u)), (1-u*conj(u))/(1+u*conj(u)), ...
    i*(u-conj(u))/(1+u*conj(u)), 0, 0];
dalpha0 = 1/(2*(1+u*conj(u))^2)*[0, 0, 0, conj(u)^2-1, u^2-1;
    0, 0, 0, 2*conj(u), 2*u;
    0, 0, 0, -1i*(1+conj(u)^2), i*(1+u^2);
    1-conj(u)^2, -2*conj(u), i*(1+conj(u)^2), 0, 0;
    1-u^2, -2*u, -i*(1+u^2), 0, 0];
L = -i*dalpha0;
g_X1_conj_X1 = transpose(X1)*L*conj_X1;
g_X1_conj_X2 = transpose(X1)*L*conj_X2;
g_X2_conj_X1 = transpose(X2)*L*conj_X1;
g_X2_conj_X2 = transpose(X2)*L*conj_X2;
% This metric g is in [X1, conj_X1, X2, conj_X2, T].
g = [0, g_X1_conj_X1, 0, g_X1_conj_X2, 0;
    g_X1_conj_X1, 0, g_X2_conj_X1, 0, 0;
    0, g_X2_conj_X1, 0, g_X2_conj_X2, 0;
    g_X1_conj_X2, 0, g_X2_conj_X2, 0, 0;
    0, 0, 0, 0, 1];
for j=1:5
    for k=1:5
        g(j,k) = complex_simple3(g(j,k),Mdw);
    end
end
% Establish the matrix h.
h11= g(1,2); h12 = g(1,4); h21 = g(3,2); h22 = g(3,4);
h = [h11, h12; h21, h22];
% Express g in [d/dx, d/dy, d/dz, d/du, d/dconj(u)]. Write it as g1.
A = [X1, conj_X1, X2, conj_X2, T];
A_inv = inv(A); A_Tinv= inv(transpose(A));
g1 = A_Tinv*g*A_inv;
for j=1:5
    for k=1:5
         g1(j,k) = complex_simple3(g1(j,k), Mdw);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma=christoffel_twistor_CR_funW(g1,h,holo_vector,antiholo_vector,T,lie_funW,M);
variable = {'holo', 'antiholo', 'T'};
Gamma_funW.holo = sym('Gamma_two_holo_%d_%d_%d',[2 2 2]);
Gamma_funW.antiholo = sym('Gamma_two_antiholo_%d_%d_%d',[2 2 2]);
Gamma_funW.T = sym('Gamma_two_T_%d_%d', [2 2]);
for j=1:2
    for k=1:2
        Gamma_funW.T(j,k) = complex_simple3(Gamma.T(j,k),Mdw);
        for m=1:2
            Gamma_funW.holo(j,k,m)=complex_simple3(Gamma.holo(j,k,m), Mdw);
            Gamma_funW.antiholo(j,k,m)=complex_simple3(Gamma.antiholo(j,k,m),Mdw);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curvature R and R_two.
R = curvature_twistor_CR_funW(g1,h,Gamma_funW,holo_vector,antiholo_vector, ...
    lie_funW,T,RVar2,CVar2,dRVar2,dCVar2,MVar2);
R_two = sym('R2_%d_%d_%d_%d', [2,2,2,2]);
for m=1:2
    for n=1:2
        for k=1:2
            for ll=1:2
                R_two(m,n,k,ll) = complex_simple3(R(m,n,k,ll), Md2w);
            end
        end
    end
end
ricci_curv = ricci_twistor_CR_model(R_two, h, M); %Change in rho here.
rho = ricci_curv{2};
rho = complex_simple3(rho, Md2w);
% phi=d2w_uu-(6*conj(u)/(1+u*conj(u)))*dw_u+(12*conj(u)^2/(1+u*conj(u))^2)*w;
% diffTemp = rho - i*(phi-conj(phi));
% complex_simple3(diffTemp,[u,w,dw_u,d2w_uu]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% latex_code = latex_curv_CR_funW(R_two, rho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complete this program up to the sigma tensor.
syms gamma real
Theta = A_inv;
for m=1:5 
    for n=1:5
        Theta(m,n) = complex_simple3(Theta(m,n), [u w]);
    end
end
omega = omega_twistor_CR_funW(Theta, Gamma_funW);
% Define omega and conjugate_omega
omega_two11 = sym('omega_two11',[1,5]);
omega_two12 = sym('omega_two12',[1,5]);
omega_two21 = sym('omega_two21',[1,5]);
omega_two22 = sym('omega_two22',[1,5]);
for k=1:5
    omega_two11(k) = complex_simple3(omega{1,1}(k),Mdw);
    omega_two12(k) = complex_simple3(omega{1,2}(k),Mdw);
    omega_two21(k) = complex_simple3(omega{2,1}(k),Mdw);
    omega_two22(k) = complex_simple3(omega{2,2}(k),Mdw);
end
omega_two = {omega_two11, omega_two12; omega_two21, omega_two22};
%
conj_omega11 =  [conj(omega_two11(1)), conj(omega_two11(2)), ...
    conj(omega_two11(3)), conj(omega_two11(5)), conj(omega_two11(4))];
conj_omega12 =  [conj(omega_two12(1)), conj(omega_two12(2)), ...
    conj(omega_two12(3)), conj(omega_two12(5)), conj(omega_two12(4))];
conj_omega21 =  [conj(omega_two21(1)), conj(omega_two21(2)), ...
    conj(omega_two21(3)), conj(omega_two21(5)), conj(omega_two21(4))];
conj_omega22 =  [conj(omega_two22(1)), conj(omega_two22(2)), ...
    conj(omega_two22(3)), conj(omega_two22(5)), conj(omega_two22(4))];
%
conj_omega_two = {conj_omega11, conj_omega12;
    conj_omega21, conj_omega22};
%
sigma0 = sigma_one_twistor_CR_funW(Theta, ricci_curv, h, omega_two, ...
    conj_omega_two, Md2w); % Change in h here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma = [dx, dy, dz, du, dconj_u, d\gamma].
% Use sigma0 and sigma to define the Fefferman metric F on C(N). 
sigma = (1/4)*[sigma0, 1];
F_X1_T = (1/4)*(sigma0*X1);
F_conj_X1_T = (1/4)*(sigma0*conj_X1);
F_X2_T = (1/4)*(sigma0*X2);
F_conj_X2_T = (1/4)*(sigma0*conj_X2);
F_T_T = (1/2)*(sigma0*T);
% F is in [X1, conj_X1, X2, conj_X2, T, d/d\gamma];
% This part (F) was rewritten by a more general form.
F = [0, 0, 0, -i, F_X1_T, 0;
    0, 0, i, 0, F_conj_X1_T, 0;
    0, i, 0, 0, F_X2_T, 0;
    -i, 0, 0, 0, F_conj_X2_T, 0;
    F_X1_T, F_conj_X1_T, F_X2_T, F_conj_X2_T, F_T_T, 1/4;
    0, 0, 0, 0, 1/4, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Express F in [d/dx, d/dy, d/dz, d/du, d/dconj(u), d/d\gamma].
Btwo = sym('Btwo',[6,6]);
for m=1:5
    for n=1:5
        Btwo(m,n) = A(m,n);
    end
end
Btwo(6,:)=[0,0,0,0,0,1];
Btwo(:,6)=[0;0;0;0;0;1];
Btwo_inv = inv(Btwo);
F1 = transpose(Btwo_inv)*F*Btwo_inv;
for j=1:6
    for k =1:6
        F1(j,k) = complex_simple3(F1(j,k), Md2w);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From now on, let u=u1+i*u2. Express F in terms of 
% [d/dx, d/dy, d/dz, d/du1, d/du2, d/d\gamma].
syms u1 u2 real
F_x_u1 = F1(1,4) + F1(1,5);
F_x_u2 = i*(F1(1,4)-F1(1,5));
F_y_u1 = F1(2,4) + F1(2,5);
F_y_u2 = i*(F1(2,4)-F1(2,5));
F_z_u1 = F1(3,4) + F1(3,5);
F_z_u2 = i*(F1(3,4)-F1(3,5));
F_gamma_u1 = F1(6,4) + F1(6,5);
F_gamma_u2 = i*(F1(6,4)-F1(6,5));
F_u1_u1 = F1(4,4) + F1(5,5) + 2*F1(4,5);
F_u1_u2 = i*(F1(4,4)-F1(5,5));
F_u2_u2 = -(F1(4,4) + F1(5,5) - 2*F1(4,5));
% F3 is the Fefferman metric in 
% [d/dx, d/dy, d/dz, d/du1. d/du2, d/d\gamma].
F3 = [F1(1,1), F1(1,2), F1(1,3), F_x_u1, F_x_u2, F1(1,6);
    F1(2,1), F1(2,2), F1(2,3), F_y_u1, F_y_u2, F1(2,6);
    F1(3,1), F1(3,2), F1(3,3), F_z_u1, F_z_u2, F1(3,6);
    F_x_u1, F_y_u1, F_z_u1, F_u1_u1, F_u1_u2, F_gamma_u1;
    F_x_u2, F_y_u2, F_z_u2, F_u1_u2, F_u2_u2, F_gamma_u2;
    F1(1,6), F1(2,6), F1(3,6), F_gamma_u1, F_gamma_u2, F1(6,6)];
% Replace u by u1 + i u2 in F3.
% In this way we obtain the Fefferman metric F0 
% in [d/dx, d/dy, d/dz, d/du1, d/du2, d/d\gamma].
syms a0 complex
F0 = sym('F0',[6,6]);
for m=1:6
    for n=1:6
        temp = 0;
        F0(m,n)= subs(F3(m,n), conj(u), a0);
        temp = subs(F0(m,n), u, u1+i*u2);
        temp = subs(temp, a0, u1-i*u2);
        temp = complex_simple3(temp, Md2w);
        F0(m,n) = temp;
        clear temp
    end
end
M_Feff = [x, y, z, u1, u2, gamma];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('New_twistor_CR_funW_data_Jan17.mat');
