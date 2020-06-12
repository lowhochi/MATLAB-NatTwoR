% twistor_NatTwoR_CR_rmW_main4.m
load('DataMain3_NatTwoR_CR_rmW.mat');
assumeAlso(gamma, 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scalar Curvature S = RiemCurv(m,n,k,ll)*F^{mll}*F^{nk}.
S = 0;
for m=1:6
    for n=1:6
        for k=1:6
            for ll=1:6
                S = S + RiemCurv(m,n,k,ll)*F0Inv(m,ll)*F0Inv(n,k);
            end
        end
    end
end
S = complex_simple3(S, symSetRCurv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ricci curvature and Weyl curvature tensor
% RicCurv(m,n) = Ric(u_m,u_n)

RicCurv = sym('RicCurv',[6 6]);
for m=1:6
    for n=1:6
        temp=0;
        for p=1:6
            for q=1:6
                temp = temp + RiemCurv(p,m,n,q)*F0Inv(p,q);
            end
        end
        RicCurv(m,n) = complex_simple3(temp, symSetRCurv);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update: Weyl(m,n,k,ll) = = Weyl(u_m,u_n,u_k,u_ll) = F0(W#(u_m,u_n)u_k,u_l)
% W_{ijkll} = R_{ijkll}+1/4*(-T_{ill}*F_{jk}+T_{ik)*F_{jll}-T_{jk}*F_{ill}
%   +T_{jl}*F_{ik})+(S/20)*(F_{ill}*F_{jk}-F_{ik}*F_{jll});

Weyl = sym('Weyl',[6 6 6 6]);
symSetWeyl = [];
for j=1:1296
    m = countIndex1296(j,1);
    n = countIndex1296(j,2);
    k = countIndex1296(j,3);
    ll = countIndex1296(j,4);
    temp1 = 1/4*RicCurv(m,k)*F0(n,ll)+1/4*RicCurv(n,ll)*F0(m,k)...
        -1/4*RicCurv(m,ll)*F0(n,k)-1/4*RicCurv(n,k)*F0(m,ll);
    temp2 = S/20*(F0(m,ll)*F0(n,k)-F0(m,k)*F0(n,ll)); 
    Weyl(m,n,k,ll) = RiemCurv(m,n,k,ll) + temp1 + temp2; 
    tempSet = symvar(Weyl(m,n,k,ll));
    symSetWeyl = union(symSetWeyl,tempSet);
end

clearvars m n p q k ll j temp temp1 temp2 tempSet
list_of_variables_Main4 = who;

save('DataMain4_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%