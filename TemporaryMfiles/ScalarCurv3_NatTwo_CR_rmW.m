% ScalarCurv3_NatTwo_CR_rmW.m ON APR 28 2019
load('DataMain5_Nat_CR_rmW_Apr27.mat');
assumeAlso([x y z u1 u2 gamma], 'real');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lengthWeyl = length(symSetWeyl);
subSetWeyl = sym('VarSetWeyl',[1,135]);
for j=1:135
    variable = string(symSetWeyl(j));
    switch variable
        case string(u)
            subSetWeyl(j)=u;
        case string(w)
            subSetWeyl(j)=w;
        case string(dw_u)
            subSetWeyl(j)=dw_u;
        case string(d2w_uu)
            subSetWeyl(j)=d2w_uu;
        case string(d3w_uuu)
            subSetWeyl(j)=d3w_uuu;
        case string(rho)
            subSetWeyl(j)=rho;
        case string(drho_u)
            subSetWeyl(j)=drho_u;
        case string(drho_conju)
            subSetWeyl(j)=drho_conju;
        case string(d2rho_uu)
            subSetWeyl(j)=d2rho_uu;
        case string(d2rho_uconju)
            subSetWeyl(j)=d2rho_uconju;
        case string(d2rho_conjuconju)
            subSetWeyl(j)=d2rho_conjuconju;
        otherwise
            subSetWeyl(j)=0;
    end
end
MVarW = [u, w, dw_u, d2w_uu, d3w_uuu, rho, drho_u, drho_conju,...
    d2rho_uu, d2rho_uconju, d2rho_conjuconju];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WeylTwo = sym('WeylTwo',[6 6 6 6]);
% WeylTwo = Weyl when (1) M is flat and (2) w = w(u).
for m=1:6
    for n=1:6
        for k=1:6
            for ll=1:6
                temp = Weyl(m,n,k,ll);
                WeylTwo(m,n,k,ll)=subs(temp,symSetWeyl,transpose(subSetWeyl));
                WeylTwo(m,n,k,ll)=complex_simple3(WeylTwo(m,n,k,ll),MVarW);
                clear temp
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save('Data_ScalarCurv3_Apr28.mat');
