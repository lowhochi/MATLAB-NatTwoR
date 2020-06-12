% WeylFlat2020_Part1_NatTwoR_CR_rmW.m
load('DataWeylFlat3_NatTwoR_CR_rmW.mat');
clearvars variableSet realVariable myCharR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fTerm1 =i/2*(d2phi_pp*conj(d3phi_ppy)-d3phi_ppy*conj(d2phi_pp))...
    + i*(d3phi_ppconjp*conj(d2phi_py)-conj(d3phi_ppconjp)*d2phi_py);
fTerm2 = -i*d2phi_pp*conj(d3phi_ppp) -i*d3phi_pyy*conj(d2phi_pp)...
    -i*conj(d2phi_pp)*d3phi_ppconjp -i*d2phi_yy*conj(d3phi_ppconjp)...
    +2*i*d3phi_pconjpy*conj(d2phi_py) +2*i*d2phi_pconjp*conj(d3phi_ppconjp);
fTerm3 = -i*d2phi_py*conj(d3phi_ppp)-i/2*d2phi_yy*conj(d3phi_ppy)...
    -i/2*d3phi_yyy*conj(d2phi_pp) +i*d2phi_pconjp*conj(d3phi_ppy)...
    -i*conj(d2phi_pp)*d3phi_pconjpy +i*conj(d2phi_py)*conj(d3phi_pyy)...
    +i*conj(d2phi_py)*conj(d3phi_ppconjp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVar201 = [phi, dphi_p, dphi_y, d2phi_pp, d2phi_pconjp, d2phi_yy];
% derivativeDict
derivativeDict.phi = [dphi_p; conj(dphi_p); dphi_y];
derivativeDict.dphi_p = [d2phi_pp; d2phi_pconjp; d2phi_py];
derivativeDict.dphi_y = [d2phi_py; conj(d2phi_py); d2phi_yy];
derivativeDict.d2phi_pp = [d3phi_ppp; d3phi_ppconjp; d3phi_ppy];
derivativeDict.d2phi_pconjp = [d3phi_ppconjp;
    d3phi_pconjpconjp; d3phi_pconjpy];
derivativeDict.d2phi_yy = [d3phi_pyy; conj(d3phi_pyy); d3phi_yyy];

% phi = phi(p,conj(p),y);
d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, conj(d2phi_pp), conj(d2phi_py);
    d2phi_py, conj(d2phi_py), d2phi_yy];
cofactor_of_d2phi = sym('cofactor',[3,3]);
for m=1:3
    for n=1:3
        cofactor_of_d2phi(m,n) = myCofactor(d2phi,m,n);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





