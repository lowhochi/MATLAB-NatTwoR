% Twistor_CR_funWv2_main4.m
load('DataMain3_CR_funWv2.mat');
assumeAlso(realVariables,'real');
syms gamma real
realVariables = [realVariables, gamma];
% Build the Fefferman metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theta = CRVectorInv, \omega_m^n in [dx, dy, dz, du, dconj(u)];
% Connection form of TW connection: Omega
Theta = CRVectorInv;
Omega = cell(2,2);
Omega{1,1} = sym('omega11',[1,5]); %\omega_1^2;
Omega{1,2} = sym('omega12',[1,5]); %\omega_1^2;
Omega{2,1} = sym('omega21',[1,5]); %\omega_2^1;
Omega{2,2} = sym('omega22',[1,5]); %\omega_2^2;

for m=1:2
    for n=1:2
        tempPart1 = zeros(1,5);
        tempPart2 = zeros(1,5);
        for k=1:2
            tempPart1 = tempPart1 + Gamma.holo(k,m,n)*Theta(2*k-1,:);
            tempPart2 = tempPart2 + Gamma.antiholo(k,m,n)*Theta(2*k,:); 
        end
        Omega{m,n} = tempPart1 + tempPart2 + Gamma.T(m,n)*Theta(5,:);
    end
end

conjOmega = cell(2,2);
conjOmega{1,1} = sym('conjOmega11',[1,5]);
conjOmega{1,2} = sym('conjOmega12',[1,5]);
conjOmega{2,1} = sym('conjOmega21',[1,5]);
conjOmega{2,2} = sym('conjOmega22',[1,5]);
for m=1:2
    for n=1:2
        conjOmega{m,n}(1) = conj(Omega{m,n}(1));
        conjOmega{m,n}(2) = conj(Omega{m,n}(2));
        conjOmega{m,n}(3) = conj(Omega{m,n}(3));
        conjOmega{m,n}(4) = conj(Omega{m,n}(5));
        conjOmega{m,n}(5) = conj(Omega{m,n}(4));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The sigma tensor
% sigma is in [dx, dy, dz, du, dconj(u), dgamma];
sigmaPart1 = [0, 0, 0, 0, 0, 1];
sigmaPart2 = i*Omega{1,1} +i*Omega{2,2};
sigmaPart3 = zeros(1,5);
for m=1:2
    for n=1:2
        sigmaPart3 = sigmaPart3 -i/2*hInv(n,m)*dh(m,n);
    end
end
sigmaPart4 = -1/12*rho*Theta(5,:);

sigma = sigmaPart1 + [sigmaPart2,0] + [sigmaPart3,0] + [sigmaPart4,0];
sigma = 1/4*sigma;
sigma0 = sigma(1:5);

F_X1_T = sigma0*X1;
F_conjX1_T = sigma0*conjX1;
F_X2_T = sigma0*X2;
F_conjX2_T = sigma0*conjX2;
F_T_T = 2*sigma0*T;

% The Fefferman metric. F = gW + 2*(alpha \odot sigma);
% F0 is in [X1, conj_X1, X2, conj_X2, T, d/d\gamma];

F0 = [0, 0, 0, -i, F_X1_T, 0;
    0, 0, i, 0, F_conjX1_T, 0;
    0, i, 0, 0, F_X2_T, 0;
    -i, 0, 0, 0, F_conjX2_T, 0;
    F_X1_T, F_conjX1_T, F_X2_T, F_conjX2_T, F_T_T, 1/4;
    0, 0, 0, 0, 1/4, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find F in [x, y, z, du, dconju, gamma];
CRVectorTwo = [conj(mu1), mu1, 0, 0, v1normv, 0;
    conj(mu2), mu2, 0, 0, v2normv, 0;
    conj(mu3), mu3, 0, 0, v3normv, 0;
    0, w, 1, 0, 0, 0;
    conj(w), 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 1];
CRVectorTwoInv = inv(CRVectorTwo);
% transpose(CRVectorTwo)*F*CRVectorTwo = F0;
F = transpose(CRVectorTwoInv)*F0*CRVectorTwoInv;

for m=1:6
    for n=1:6
        temp = F(m,n);
        temp = complex_simple3(temp, MVarMain3);
        F(m,n) = temp;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars m n k temp tempPart1 tempPart2 sigma0 sigmaPart1
clearvars sigmaPart2 sigmaPart3 sigmaPart4
clearvars F_X1_T F_X2_T F_T_T F_conjX1_T F_conjX2_T

save('DataMain4_CR_funWv2.mat')