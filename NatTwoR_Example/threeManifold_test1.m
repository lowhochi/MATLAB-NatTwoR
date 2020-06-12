% orthonormal basis: e1 e2 e3 n
% define Christoffel symbols:
% \nabla_{e_i}e_j = Gij_k*e_k;
% \nabla_{e_i}n = -Gik_0*e_k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms u
syms G12_3 G23_1 G31_2 G11_2 G11_3 G22_1 G22_3 G33_1 G33_2
syms G11_0 G12_0 G13_0 G22_0 G23_0 G33_0
MVar = [u, G12_3, G23_1, G31_2, G11_2, G11_3, G22_1, G22_3, G33_1, G33_2,...
    G11_0, G12_0, G13_0, G22_0, G23_0, G33_0];
epsilon = zeros(3,3,3);
epsilon(1,2,3) = 1;
epsilon(2,3,1) = 1;
epsilon(3,1,2) = 1;
epsilon(1,3,2) = -1;
epsilon(2,1,3) = -1;
epsilon(3,2,1) = -1;
% for j=1:3
%     for k=1:3
%         for m=1:3
%             disp([j,k,m]);
%             disp(epsilon(j,k,m));
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma = sym('Gamma',[3,3,3]);
Gamma(1,1,1) = 0;
Gamma(1,1,2) = G11_2;
Gamma(1,1,3) = G11_3;
Gamma(1,2,1) = -G11_2;
Gamma(1,2,2) = 0;
Gamma(1,2,3) = G12_3;
Gamma(1,3,1) = -G11_3;
Gamma(1,3,2) = -G12_3;
Gamma(1,3,3) = 0;
Gamma(2,1,1) = 0;
Gamma(2,1,2) = -G22_1;
Gamma(2,1,3) = -G23_1;
Gamma(2,2,1) = G22_1;
Gamma(2,2,2) = 0;
Gamma(2,2,3) = G22_3;
Gamma(2,3,1) = G23_1;
Gamma(2,3,2) = -G22_3;
Gamma(2,3,3) = 0;
Gamma(3,1,1) = 0;
Gamma(3,1,2) = G31_2;
Gamma(3,1,3) = -G33_1;
Gamma(3,2,1) = -G31_2;
Gamma(3,2,2) = 0;
Gamma(3,2,3) = -G33_2;
Gamma(3,3,1) = G33_1;
Gamma(3,3,2) = G33_2;
Gamma(3,3,3) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GammaN = [G11_0, G12_0, G13_0;
    G12_0, G22_0, G23_0;
    G13_0, G23_0, G33_0];

% Define T(i,j,k) = Tij_k.
T = sym('T',[3,3,3]);
for ii=1:3
    for j=1:3
        for k=1:3
            temp = 0;
            for ll=1:3
                temp = temp+GammaN(j,ll)*epsilon(ii,ll,k)...
                    -GammaN(ii,ll)*epsilon(j,ll,k);
            end  
            T(ii,j,k) = complex_simple3(temp, MVar);
        end
    end
end
traceT = sym('traceT',[1,3]);
for ii=1:3
    temp = 0;
    for k=1:3
        temp = temp + T(ii,k,k);
    end
    traceT(ii) = complex_simple3(temp,MVar);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = T(1,2,3)-T(1,3,2)+T(2,3,1)-T(2,1,3)+T(3,1,2)-T(3,2,1);
tau = complex_simple3(tau ,MVar);

qRow = sym('qRow',[3,3,3]);
for ii=1:3
    for j=1:3
        for k=1:3
            qRow(ii,j,k) = T(ii,j,k)...
                -1/2*(traceT(ii)*myDelta(j,k)-traceT(j)*myDelta(ii,k))...
                -tau/6*epsilon(ii,j,k);
            qRow(ii,j,k) = complex_simple3(qRow(ii,j,k), MVar);
        end
    end
end
Q = sym('Q',[3,3]);            
for k=1:3
    for ll=1:3
        temp = 0;
        for ii=1:3
            for j=1:3
                temp = temp + 1/2*epsilon(ii,j,ll)*qRow(ii,j,k);
            end
        end
        Q(k,ll) = complex_simple3(temp, MVar);
    end
end
clearvars k ll ii j temp part1 part2           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
C = [i, 0, 1, i/2, 0;
    0, -2*i, 0, 0, -2;
    0, 0, 0, 3*i, 0;
    0, 2*i, 0, 0, -2;
    i, 0, -1, i/2, 0];
uVec = [1; u; u^2; u^3; u^4];
qVec = [Q(1,1); Q(1,2); Q(1,3); Q(2,2); Q(2,3)];

uCoeff = C*qVec;
for j=1:5
    uCoeff(j) = complex_simple3(uCoeff(j),MVar);
end
w = transpose(uVec)*uCoeff;
w = complex_simple3(w, MVar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
