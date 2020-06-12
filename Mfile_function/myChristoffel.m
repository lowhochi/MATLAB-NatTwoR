function G = myChristoffel(gMat, eMat, varSet)
% varSet = [x, y, z];
% G(m,n,k) = \Gamma_{mn}^k corresponding to {e1, e2, e3}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e_i = eMat(k,i)*d/dx_k = a_{ki}*d/dx_k
% (1) [e_i, e_j] = A_ij^k*d/dx_k;
% (2) [e_i, e_j] = B_ij^k*e_k;
% B_ij^k = eMat(ll,k)*A(i,j,k);
Imat = eye(3);
A = sym('A',[3, 3, 3]);
B = sym('B',[3, 3, 3]);
for ii=1:3
    for j=1:3
        for k=1:3
            part1 = 0;
            part2 = 0;
            for m=1:3
                dakj_by_xm = diff(eMat(k,j), varSet(m));
                daki_by_xm = diff(eMat(k,ii), varSet(m));
                part1 = part1 + eMat(m,ii)*dakj_by_xm;
                part2 = part2 + eMat(m,j)*daki_by_xm;
            end
            A(ii,j,k) = part1 - part2;
            A(ii,j,k) = simplify(A(ii,j,k));
        end
    end
end
for ii=1:3
    for j=1:3
        for k=1:3
            Bijk = 0;
            for ll=1:3
                temp = gMat(ll,:)*eMat(:,k);
                Bijk = Bijk + A(ii,j,ll)*temp;
            end
            B(ii,j,k) = simplify(Bijk);
        end
    end
end
clearvars temp Bijk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = sym('G',[3,3,3]);
% G(i,j,k) = 1/2*g^{k,ll}*
%   (-B(i,ll,p)*g(j,p) -B(j,ll,p)*g(i,p) +B(i,j,p)*g(ll,p)); 
for ii=1:3
    for j=1:3
        for k=1:3
            temp = 0;
            for ll=1:3
                part1 = 0;
                part2 = 0;
                part3 = 0;
                for p=1:3
                    part1 = part1 - B(ii,ll,p)*Imat(j,p);
                    part2 = part2 - B(j,ll,p)*Imat(ii,p);
                    part3 = part3 + B(ii,j,p)*Imat(ll,p);
                end
                temp = temp + 1/2*Imat(k,ll)*(part1 + part2 + part3);
            end
            G(ii,j,k)=simplify(temp);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
