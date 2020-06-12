% riem_curvature_Feff_CR_rmW.m
function Riem_curv = riem_curvature_Feff_CR_rmW(F0,Riem_Gamma,...
    dRiem_Gamma,variableXYZ)
Riem_curv = sym('R_%d_%d_%d_%d',[6 6 6 6]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R_kjml = g(R(d/dx_k, d/dx_j)d/dx_m, d/dx_l).
for k=1:6
    for j=1:6
        for m=1:6
            for ll=1:6
                term = 0;
                for s=1:6
                    temp2 = 0;
                    temp1 = dRiem_Gamma.(variableXYZ{k})(m,j,s) - ...
                        dRiem_Gamma.(variableXYZ{j})(m,k,s);
                    for r=1:6
                        temp2 = temp2 + Riem_Gamma(m,j,r)*Riem_Gamma(r,k,s) - ...
                        Riem_Gamma(m,k,r)*Riem_Gamma(r,j,s);
                    end
                    term = term + F0(s,ll)*(temp1+temp2);
                    clear temp1 
                    clear temp2
                end
                Riem_curv(k,j,m,ll) = term;
                clear term
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%