% conformalFlatModel_NatTwoR_CR_rmW.m
% metric g=f^2*(dx^2+dy^2+dz^2);
% orthonormal frame: e1, e2, e3
% e1 = 1/f*d/dx; e2 = 1/f*d/dy; e3 = 1/f*d/dz;
% G(i,j,k) = G_{ij}^k = g(\nabla_{e_i}e_j,e_k);
syms x y z
syms f df_x df_y df_z 
syms d2f_xx d2f_xy d2f_xz d2f_yy d2f_yz d2f_zz
syms lambda mu
syms dlambda_x dlambda_y dlambda_z
syms d2lambda_xx d2lambda_xy d2lambda_xz d2lambda_yy
syms d2lambda_yz d2lambda_zz
syms dmu_x dmu_y dmu_z
syms d2mu_xx d2mu_xy d2mu_xz d2mu_yy d2mu_yz d2mu_zz
CVar = [x, y, z, f, lambda, mu, dlambda_x, dlambda_y, dlambda_z,...
    dmu_x, dmu_y, dmu_z, df_x, df_y, df_z];

MVar = [f, lambda, mu, dlambda_x, dlambda_y, dlambda_z, ...
    dmu_x, dmu_y, dmu_z, df_x, df_y, df_z,...
    d2f_xx d2f_xy d2f_xz d2f_yy d2f_yz d2f_zz, ...
    d2lambda_xx d2lambda_xy d2lambda_xz d2lambda_yy,...
    d2lambda_yz d2lambda_zz, ...
    d2mu_xx d2mu_xy d2mu_xz d2mu_yy d2mu_yz d2mu_zz];

derivativeDict.x = [1; 0; 0];
derivativeDict.y = [0; 1; 0];
derivativeDict.z = [0; 0; 1];
derivativeDict.f = [df_x; df_y; df_z];
derivativeDict.lambda = [dlambda_x; dlambda_y; dlambda_z];
derivativeDict.mu = [dmu_x; dmu_y; dmu_z];
derivativeDict.df_x = [d2f_xx; d2f_xy; d2f_xz];
derivativeDict.df_y = [d2f_xy; d2f_yy; d2f_yz];
derivativeDict.df_z = [d2f_xz; d2f_yz; d2f_zz];
derivativeDict.dlambda_x = [d2lambda_xx; d2lambda_xy; d2lambda_xz];
derivativeDict.dlambda_y = [d2lambda_xy; d2lambda_yy; d2lambda_yz];
derivativeDict.dlambda_z = [d2lambda_xz; d2lambda_yz; d2lambda_zz];
derivativeDict.dmu_x = [d2mu_xx; d2mu_xy; d2mu_xz];
derivativeDict.dmu_y = [d2mu_xy; d2mu_yy; d2mu_yz];
derivativeDict.dmu_z = [d2mu_xz; d2mu_yz; d2mu_zz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = sym('G',[3 3 3]);
G(1,1,1) = 0;
G(1,1,2) = -df_y/f^2;
G(1,1,3) = -df_z/f^2;
G(1,2,1) = -G(1,1,2);
G(1,2,2) = 0;
G(1,2,3) = 0;
G(1,3,1) = -G(1,1,3);
G(1,3,2) = 0;
G(1,3,3) = 0;
G(2,1,1) = 0;
G(2,1,2) = df_x/f^2;
G(2,1,3) = 0;
G(2,2,1) = -df_x/f^2;
G(2,2,2) = 0;
G(2,2,3) = -df_z/f^2;
G(2,3,1) = 0;
G(2,3,2) =  -G(2,2,3);
G(2,3,3) = 0;
G(3,1,1) = 0;
G(3,1,2) = 0;
G(3,1,3) = df_x/f^2;
G(3,2,1) = 0;
G(3,2,2) = 0;
G(3,2,3) = df_y/f^2;
G(3,3,1) = -df_x/f^2;
G(3,3,2) = -df_y/f^2;
G(3,3,3) = 0;
% define another orthonormal frame on M3;
% u_i = a_{ki}e_k;
% A = [a_{ki}] 3x3;
A = sym('A', [3,3]);
A(1,1) = 1/2*(lambda^2 -mu^2 +conj(lambda)^2 -conj(mu)^2);
A(1,2) = -lambda*conj(mu) -conj(lambda)*mu;
A(1,3) = -i/2*(lambda^2 -mu^2 +conj(mu)^2 -conj(lambda)^2);
A(2,1) = lambda*mu +conj(lambda)*conj(mu);
A(2,2) = lambda*conj(lambda) -mu*conj(mu);
A(2,3) = -i*(lambda*mu -conj(lambda)*conj(mu));
A(3,1) = i/2*(lambda^2 +mu^2 -conj(lambda)^2 -conj(mu)^2);
A(3,2) = i*(-lambda*conj(mu) +conj(lambda)*mu);
A(3,3) = 1/2*(lambda^2 +mu^2 +conj(lambda)^2 +conj(mu)^2);
for j=1:3
    for k=1:3
        A(j,k) = 1/(lambda*conj(lambda)+mu*conj(mu))*A(j,k);
    end
end
% % check inv(A) = transpose(A)
% Ainv = inv(A);
% B = A*transpose(A);
% for j=1:3
%     for k=1:3
%         B(j,k) = complex_simple3(B(j,k),MVar);
%     end
% end
% Ematrix = [E1, E2, E3];
Ematrix = [1/f, 0, 0;
    0, 1/f, 0;
    0, 0, 1/f];

% define T(m,n,k) = g(nabla_{f_m}f_n, f_k);
T = sym('T',[3 3 3]);
for m=1:3
    for n=1:3
        for k=1:3
            firstPart = 0;
            for ll=1:3
                da_lln = df_conformalFlat_CR_rmW(A(ll,n),CVar,derivativeDict);
                for p=1:3
                    temp01 = da_lln*Ematrix(:,p);
                    firstPart = firstPart + A(p,m)*A(ll,k)*temp01;
                end
            end
            secondPart = 0;
            for p=1:3
                for q=1:3
                    for ll=1:3
                        secondPart = secondPart...
                            + A(p,m)*A(ll,n)*A(q,k)*G(p,ll,q);
                    end
                end
            end
            temp = firstPart + secondPart;
            T(m,n,k) = complex_simple3(temp, MVar);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = T(1,2,3) + T(2,3,1) + T(3,1,2);
theta = complex_simple3(theta, MVar);
dtheta = df_conformalFlat_CR_rmW(theta,CVar,derivativeDict);

dtheta_byf1 = 0;
dtheta_byf2 = 0;
dtheta_byf3 = 0;
for m=1:3
    temp02 = dtheta*Ematrix(:,m);
    dtheta_byf1 = dtheta_byf1 + A(m,1)*temp02;
%     dtheta_byf2 = dtheta_byf2 + A(m,2)*temp02;
%     dtheta_byf3 = dtheta_byf3 + A(m,3)*temp02;
end

dtheta_byf1 = complex_simple3(dtheta_byf1, MVar);
% dtheta_byf2 = complex_simple3(dtheta_byf2, MVar);
% dtheta_byf3 = complex_simple3(dtheta_byf3, MVar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(test1);
