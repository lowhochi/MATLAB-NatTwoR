% twistor_NatTwoR_CR_rmW_main1.m
syms u w
syms h11 T4 aM aV bV % Change aMu to aM
syms dw_Mu dw_conjMu dw_vnormv dw_u
syms dT4_Mu dT4_conjMu dT4_vnormv dT4_u dT4_conju
syms dh11_Mu dh11_conjMu dh11_vnormv dh11_u dh11_conju
% NatTwo: build the model by variables h11 T4 aM aV bV. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
v1normv = complex_simple3(v1normv, [u,w]);
v2normv = complex_simple3(v2normv, [u,w]);
v3normv = complex_simple3(v3normv, [u,w]);
conj_X1 = [u^2-1; 2*u; i*(u^2+1); w; 0];
X1 = [conj(u)^2-1; 2*conj(u); -i*(conj(u)^2+1); 0; conj(w)];
X2 = [0;0;0;1;0];
conj_X2 = [0;0;0;0;1];
T = [v1normv; v2normv; v3normv; T4; conj(T4)];
h = [h11, -i;
    i, 0];
hInv = inv(h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Webster metric g in the basis (X1, conj_X1, X2, conj_X2, T)
g = [0, h11, 0, -i, 0;
    h11, 0, i, 0, 0;
    0, i, 0, 0, 0;
    -i, 0, 0, 0, 0;
    0, 0, 0, 0, 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lieOne = cell(5,5);
% Lie vectors are in (X1, conj_X1, X2, conj_X2, T)
% lieOne{1,:} = {0, [X1,conj_X1], [X1,X2], [X1,conj_X2], [X1,T]};
% lieOne{2,:} = {-[X1,conj_X1], 0, [conj_X1,X2], [conj_X1,conj_X2], [conj_X1,T]};
% lieOne{3,:} = {-[X1,X2], -[conj_X1,X2], 0, [X2,conj_X2], [X2,T]};
% lieOne{4,:} = {-[X1,conj_X2], -[conj_X1,conj_X2], -[X2,conj_X2], 0, [conj_X2,T]};
% lieOne{5,:} = {-[X1,T], -[conj_X1,T], -[X2,T], -[conj_X2,T], 0};
lie_X1_conj_X1 = [conj(aM);
    -aM; 
    aM*w + 2*i*h11*T4 + dw_conjMu;
    -conj(aM)*conj(w) + 2*i*h11*conj(T4) - conj(dw_conjMu);
    -2*i*h11];
% [X1,conj_X1]= conj(aM)*X1 - aM*conj_X1 + (aM*w + 2*i*h11*T4 + dw_conjMu)*X2
%   + (-conj(aM)*conj(w) + 2*i*h11*conj(T4) - conj(dw_conjMu))*conj_X2 - 2*i*h11*T;
lie_X1_conj_X2 = [-2*u/(1+u*conj(u));
    0;
    2*T4;
    2*u*conj(w)/(1+u*conj(u)) + 2*conj(T4) - conj(dw_u); 
    -2];
lie_conj_X1_T = [bV - w/(1+u*conj(u))^2; 
    aV - 2*conj(u)/(1+u*conj(u))*T4;
    -dw_vnormv - T4*dw_u + (-aV + dT4_u + 2*conj(u)/(1+u*conj(u))*T4)*w + dT4_Mu;
    w*conj(w)/(1+u*conj(u))^2 - bV*conj(w) + conj(dT4_conjMu) + w*conj(dT4_conju);     
    0];
lie_X2_T = [-1/(1+u*conj(u))^2;
    0;
    dT4_u;
    conj(w)/(1+u*conj(u))^2 + conj(dT4_conju);
    0];
lie_conj_X1_X2 = [conj(lie_X1_conj_X2(2)); 
    conj(lie_X1_conj_X2(1));
    conj(lie_X1_conj_X2(4));
    conj(lie_X1_conj_X2(3));
    conj(lie_X1_conj_X2(5))];
lie_X1_T = [conj(lie_conj_X1_T(2)); 
    conj(lie_conj_X1_T(1));
    conj(lie_conj_X1_T(4));
    conj(lie_conj_X1_T(3));
    conj(lie_conj_X1_T(5))];
lie_conj_X2_T = [conj(lie_X2_T(2)); 
    conj(lie_X2_T(1));
    conj(lie_X2_T(4)); 
    conj(lie_X2_T(3));
    conj(lie_X2_T(5))];
%
lieOne{1,1} = zeros(5,1);
lieOne{1,2} = lie_X1_conj_X1;
lieOne{1,3} = zeros(5,1);
lieOne{1,4} = lie_X1_conj_X2;
lieOne{1,5} = lie_X1_T;
lieOne{2,1} = -lie_X1_conj_X1;
lieOne{2,2} = zeros(5,1);
lieOne{2,3} = lie_conj_X1_X2;
lieOne{2,4} = zeros(5,1);
lieOne{2,5} = lie_conj_X1_T;
lieOne{3,1} = zeros(5,1);
lieOne{3,2} = -lie_conj_X1_X2;
lieOne{3,3} = zeros(5,1);
lieOne{3,4} = zeros(5,1);
lieOne{3,5} = lie_X2_T;
lieOne{4,1} = -lie_X1_conj_X2;
lieOne{4,2} = zeros(5,1);
lieOne{4,3} = zeros(5,1);
lieOne{4,4} = zeros(5,1);
lieOne{4,5} = lie_conj_X2_T;
lieOne{5,1} = -lie_X1_T;
lieOne{5,2} = -lie_conj_X1_T;
lieOne{5,3} = -lie_X2_T;
lieOne{5,4} = -lie_conj_X2_T;
lieOne{5,5} = zeros(5,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Christoffel symbols: Gamma.holo, Gamma.antiholo, Gamma.T
% CRvector = [X1, conj_X1, X2, conj_X2, T] in (Mu,conjMu,vnormv,u,conju)
CRvector = sym('CRvector',[5 5]);
CRvector(:,1) = [0; 1; 0; 0; conj(w)];
CRvector(:,2) = [1; 0; 0; w; 0];
CRvector(:,3) = [0; 0; 0; 1; 0];
CRvector(:,4) = [0; 0; 0; 0; 1];
CRvector(:,5) = [0; 0; 1; T4; conj(T4)];
%
dh = cell(2,2);
dh{1,1} = [dh11_Mu, dh11_conjMu, dh11_vnormv, dh11_u, dh11_conju];
dh{1,2} = zeros(1,5);
dh{2,1} = zeros(1,5);
dh{2,2} = zeros(1,5);
%
Gamma.holo = sym('Gamma_holo_%d_%d_%d',[2 2 2]);
Gamma.antiholo = sym('Gamma_antiholo_%d_%d_%d',[2 2 2]);
Gamma.T = sym('Gamma_T_%d_%d', [2 2]);
for m=1:2
    for n=1:2
        for k=1:2
            temp1 = 0;
            temp2 = 0;
            for ll=1:2
                temp1 = temp1 + hInv(ll,k)*(dh{n,ll}*CRvector(:,2*m-1)...
                    - g(2*n-1,:)*lieOne{2*m-1,2*ll});
                temp2 = temp2 + hInv(ll,k)*(g(2*ll,:)*lieOne{2*m, 2*n-1});
            end
            Gamma.holo(m,n,k) = temp1;
            Gamma.antiholo(m,n,k) = temp2;
            clearvars temp1 temp2
        end
    end
end
for n=1:2
    for k=1:2
        tempT = 0;
        for ll=1:2
            tempT = tempT + hInv(ll,k)*(g(2*ll,:)*lieOne{5,2*n-1});
        end
        Gamma.T(n,k) = tempT;
        clear tempT
    end
end

MVarMain1 = []; % complex simplify Gamma.holo, Gamma.antiholo and Gamma.T.
for n=1:2
    for k=1:2
        tempSetT = symvar(Gamma.T(n,k));
        MVarMain1 = union(MVarMain1, tempSetT);
        for m=1:2
            tempSetH = symvar(Gamma.holo(m,n,k));
            tempSetA = symvar(Gamma.antiholo(m,n,k));
            MVarMain1 = union(MVarMain1, tempSetH);
            MVarMain1 = union(MVarMain1, tempSetA);
        end
        clearvars tempSetT tempSetH tempSetA
    end
end        

for n=1:2
    for k=1:2
        Gamma.T(n,k) = complex_simple3(Gamma.T(n,k), MVarMain1);
        for m=1:2
            Gamma.holo(m,n,k)=complex_simple3(Gamma.holo(m,n,k),MVarMain1);
            Gamma.antiholo(m,n,k)=complex_simple3(Gamma.antiholo(m,n,k),MVarMain1);
        end
    end
end

clearvars m n k ll
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('DataMain1_NatTwoR_CR_rmW.mat');