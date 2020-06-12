function cotton = myCottonTensor(eMat, G, RmThree, varSet)
% varSet = [x,y,z];
% simple cases: phi = (x+y*z), (1+x^2+y*z), 1+x^2+y^2+z^2;
%   phi = 3+z^2+sin(x)*cos(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ithree = eye(3);
ricci = sym('ricci',[3 3]);
scalar = 0;
for ii=1:3
    for j=1:3
        temp = 0;
        for k=1:3
            for ll=1:3
                temp = temp + Ithree(k,ll)*RmThree(k,ii,j,ll);
            end
        end
        ricci(ii,j) = simplify(temp);
        scalar = scalar + temp*Ithree(ii,j);
        clear temp
    end
end
clearvars ii j k ll
scalar = simplify(scalar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% schouten P(ii,j) = ricci(ii,j) - scalar/4*g(e_i, e_j);
P = sym('P',[3 3]);
cotton = sym('cotton',[3 3 3]);
for ii=1:3
    for j=1:3
        P(ii,j) = ricci(ii,j) - (scalar/4)*myDelta(ii,j);
        P(ii,j) = simplify(P(ii,j));
    end
end
% cotton(ii,j,k) = (\nabla_k P)iij - (\nabla_j P)iik;
for ii=1:3
    for j=1:3
        for k=1:3
            dricci_ij_by_ek = 0;
            dricci_ik_by_ej = 0;
            dscalar_by_ek = 0;
            dscalar_by_ej = 0;
            Nk_Pij_part3 = 0;
            Nj_Pik_part3 = 0;
            for m=1:3
                dricci_ij_by_ek = dricci_ij_by_ek ...
                    + diff(ricci(ii,j),varSet(m))*eMat(m,k);
                dricci_ik_by_ej = dricci_ik_by_ej ...
                    + diff(ricci(ii,k),varSet(m))*eMat(m,j);
                dscalar_by_ek = dscalar_by_ek ...
                    + diff(scalar,varSet(m))*eMat(m,k);               
                dscalar_by_ej = dscalar_by_ej ...
                    + diff(scalar,varSet(m))*eMat(m,j);
                % % % % %
                Nk_Pij_part3 = Nk_Pij_part3 ...
                    -G(k,ii,m)*P(m,j) -G(k,j,m)*P(m,ii);
                Nj_Pik_part3 = Nj_Pik_part3 ...
                    -G(j,ii,m)*P(m,k) -G(j,k,m)*P(m,ii);
            end
            Nk_Pij = dricci_ij_by_ek - 1/4*dscalar_by_ek*myDelta(ii,j)...
                + Nk_Pij_part3; 
            Nj_Pik = dricci_ik_by_ej -1/4*dscalar_by_ej*myDelta(ii,k)...
                + Nj_Pik_part3;
            cotton(ii,j,k) = Nk_Pij - Nj_Pik;
            cotton(ii,j,k) = simplify(cotton(ii,j,k));
        end
    end
end