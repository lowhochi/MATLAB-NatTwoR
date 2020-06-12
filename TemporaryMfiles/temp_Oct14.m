load('DataMain1Ch_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rtemp = sym('Rtemp',[2 2 2 2]);

variableSet01 = [T4, aM, aV, bV, d2h11_MuconjMu, d2h11_conjuMu, d2h11_conjuconjMu,...
    d2h11_conjuconju, d2h11_uMu, d2h11_uconjMu, d2h11_uconju,...
    d2h11_uu, d2w_conjMuconjMu, d2w_uconjMu, d2w_uu, dT4_conjMu,...
    dT4_conju, dT4_u, daM_Mu, daM_conjMu, daM_conju, daM_u,...
    dh11_Mu, dh11_conjMu, dh11_conju, dh11_u, dh11_vnormv,...
    dw_conjMu, dw_u, dw_vnormv, h11, u, w];
 subSet01 = sym('subs',[1,length(variableSet01)]);
 MVar01 = [u, w, dw_u, dw_conjMu, dw_vnormv, d2w_uu,...
     d2w_uconjMu, d2w_conjMuconjMu];
 
for j=1:length(variableSet01)
    myString = string(variableSet01(j));
    switch myString
        case string(u)
            subSet01(j) = u;
        case string(w)
            subSet01(j) = w;
        case string(dw_conjMu)
            subSet01(j) = dw_conjMu;
        case string(dw_vnormv)
            subSet01(j) = dw_vnormv;
        case string(dw_u)
            subSet01(j) = dw_u;
        case string(d2w_uu)
            subSet01(j) = d2w_uu;
        case string(d2w_uconjMu)
            subSet01(j) = d2w_uconjMu;
        case string(d2w_conjMuconjMu)
            subSet01(j) = d2w_conjMuconjMu;
        otherwise 
            subSet01(j) = 0;
    end
end



for m=1:2
    for n=1:2
        for k=1:2
            for ll=1:2
                temp = R(m,n,k,ll);
                temp = subs(temp, variableSet01, subSet01);
                temp = complex_simple3(temp,MVar01);
                Rtemp(m,n,k,ll) = temp;
            end
        end
    end
end
clearvars m n k ll temp j myString

test2211 = -Rtemp(2,2,1,1) -d2w_uconjMu +4*conj(u)/(1+u*conj(u))*dw_conjMu...
    +2*w*conj(w)/(1+u*conj(u))^2;
test2211 = complex_simple3(test2211,MVar01);
