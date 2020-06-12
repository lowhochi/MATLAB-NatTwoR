% example3_NatTwoR_CR_rmW.m
% load('DataWeyl2_NatTwoR_CR_rmW.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the cofactor condition is satisfied under a \phi function
syms p y c0
phi = c0*p*y +conj(c0)*conj(p)*y;
dphi_p = complexdiff3(phi, p, 0);
dphi_conjp = complexdiff3(phi, p, 1);
dphi_y = diff(phi, y);
% % % % %
d2phi_pp = complexdiff3(dphi_p, p, 0);
d2phi_pconjp = complexdiff3(dphi_p, p, 1);
d2phi_conjpconjp = complexdiff3(dphi_conjp, p, 1);
d2phi_py = diff(dphi_p, y);
d2phi_conjpy = diff(dphi_conjp, y);
d2phi_yy = diff(dphi_y, y);
% % % % %
d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, d2phi_conjpconjp, d2phi_conjpy;
    d2phi_py, d2phi_conjpy, d2phi_yy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the cofactor matrix of d2phi
Cmat = sym('Cofactor',[3 3]);
for ii=1:3
    for jj=1:3
        temp = myCofactor(d2phi,ii,jj);
        Cmat(ii,jj) = temp;
    end
end
% Find fTerm1, fTerm2 and fTerm3
fTerm1 = i/2*(complexdiff3(Cmat(1,3),p,0)-complexdiff3(Cmat(2,3),p,1));
fTerm2 = -i*(complexdiff3(Cmat(1,1),p,0)-complexdiff3(Cmat(3,3),p,1));
fTerm3 = i*complexdiff3(Cmat(1,3),p,1)-i/2*diff(Cmat(1,1),y);
fTermVec = [fTerm1, fTerm2, fTerm3];
for jj=1:3
    temp = fTermVec(jj);
    temp = subs(temp, conj(y), y);
    fTermVec(jj) = complex_simple3(temp, [p,y,c0]); 
end
clearvars temp ii jj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
syms p y phi dphi_p dphi_y
dphi_conjp = conj(dphi_p);
syms d2phi_pp d2phi_pconjp d2phi_py d2phi_yy
d2phi_conjpconjp = conj(d2phi_pp);
d2phi_conjpy = conj(d2phi_py);
MVar = [p, y, d2phi_pp, d2phi_pconjp, d2phi_py, d2phi_yy];
realVariable = [y, dphi_y, d2phi_pconjp,d2phi_yy];
% % % % %
d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, d2phi_conjpconjp, d2phi_conjpy;
    d2phi_py, d2phi_conjpy, d2phi_yy];
detD2phi = det(d2phi);
% % % % %
Cmat = sym('Cmat',[3 3]);
CTwoMat = sym('CTwoMat',[3 3]);
testMat = sym('test',[3 3]);
for j=1:3
    for k=1:3
        Cmat(j,k) = myCofactor(d2phi,j,k);
    end
end
for j=1:3
    for k=1:3
        temp = myCofactor(Cmat,j,k);
        for m=1:4
            myVar = realVariable(m);
            temp = subs(temp, conj(myVar), myVar);
        end
        CTwoMat(j,k) = complex_simple3(temp,MVar);
        testMat(j,k) = complex_simple3(CTwoMat(j,k)/detD2phi,MVar);
    end
end




