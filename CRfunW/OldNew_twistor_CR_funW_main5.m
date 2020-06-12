% New_twistor_CR_funW_main5.m
load('New_twistor_CR_funW_RicS_data_Jan17.mat');
assumeAlso([x y z u1 u2 gamma],'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ric = Ric_S_Weyl{1};
S = Ric_S_Weyl{2};
Weyl = sym('weyl_%d_%d_%d_%d',[6 6 6 6]);
for m=1:6
    for n=1:6
        for k=1:6
             for ll=1:6
                 temp1 = -Riem_curv(k,ll,m,n);
                 temp2 = (-1/4)*(Ric(m,k)*F4(n,ll) + Ric(n,ll)*F4(m,k) ...
                     - Ric(m,ll)*F4(n,k) - Ric(n,k)*F4(m,ll));
                 temp3 = (S/20)*(F4(m,k)*F4(n,ll)-F4(m,ll)*F4(n,k));
                 Weyl(m,n,k,ll) = temp1 + temp2 + temp3;
                 clearvars temp1 temp2 temp3
             end
        end
    end
end
Ric_S_Weyl{3} = Weyl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('New_twistor_CR_funW_RicSWeyl_data_Jan21.mat');