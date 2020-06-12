% New_twistor_CR_funW_main4.m
load('New_twistor_CR_funW_differentiate_Riem_Gamma_data.mat');
assumeAlso([x y z u1 u2 gamma],'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Find the Riemannian curvature tensor field by the m-file,
% riem_curvature_Feff_CR_rmW.m.
Riem_curv = riem_curvature_Feff_CR_rmW(F0,Riem_Gamma,...
    dRiem_Gamma,variableXYZ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:6
    for n=1:6
        for k=1:6
            for ll=1:6
                tempRC = Riem_curv(m,n,k,ll);
                tempRC = subs(tempRC, u1, (u+conj(u))/2);
                tempRC = subs(tempRC, u2, -i*(u-conj(u))/2);
                Riem_curv(m,n,k,ll) = tempRC;
                clear tempRC
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F4 = sym('F4',[6 6]);
for m=1:6
    for n=1:6
        temp = F0(m,n);
        temp = subs(temp, u1, (u+conj(u))/2);
        temp = subs(temp, u2, -i*(u-conj(u))/2);
        F4(m,n) = temp;
        clear temp
    end
end
Ric_S_Weyl = ricci_weyl_twistor_CR_funW(Riem_curv, F4, MVar4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('New_twistor_CR_funW_RicS_data_Jan17.mat')