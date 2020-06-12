function RmThree = myRiemThreeMfd(eMat, G, varSet)
% gMat = metric wrt d/dx, d/dy, d/dz;
% eMat = orthonormal basis in [e1, e2, e3];
% G = Christoffel symbols wrt {e1, e2, e3};
% varSet = [x, y, z];
% Output RmThree(ii,j,k,ll) = <Rm(e_ii,e_j)e_k, e_ll> ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R(e_i, e_j)e_k = [e_i(G_jk^l) - e_j(G_ik^l)...
%   + G_jk^m*G_im^l - G_ik^m*G_jm^l + G_ji^m*G_mk^l - G_ij^m*G_mk^l]*e_l;
% e_i = eMat(m,i)*d/dx_m = a_{mi}*d/dx_m
RmThree = sym('RmThree',[3,3,3,3]);
for ii=1:3
    for j=1:3
        for k=1:3
            for ll=1:3
                dGjk_l_by_ei = 0;
                dGik_l_by_ej = 0;
                part2 = 0;
                for m=1:3
                    dGjk_l_by_ei = dGjk_l_by_ei ...
                        + diff(G(j,k,ll), varSet(m))*eMat(m,ii);
                    dGik_l_by_ej = dGik_l_by_ej ...
                        + diff(G(ii,k,ll), varSet(m))*eMat(m,j);
                    part2 = part2 ...
                        + G(j,k,m)*G(ii,m,ll)- G(ii,k,m)*G(j,m,ll)...
                        + G(j,ii,m)*G(m,k,ll) - G(ii,j,m)*G(m,k,ll);
                end
                temp = dGjk_l_by_ei - dGik_l_by_ej + part2;
                RmThree(ii,j,k,ll) = simplify(temp);
            end
        end
    end
end

