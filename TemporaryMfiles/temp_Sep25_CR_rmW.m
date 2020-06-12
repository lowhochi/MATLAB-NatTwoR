load('Data_ChernThree_Section2_Sep24_2019.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/ Check Chern(1,1,1,2).
Ch1112 = Chern_in_aV(1,1,1,2);
test1112p1 = i/2*(-d2w_uconjMu+conj(d2w_uconjMu))...
    + 2i/Y0*(conj(u)*dw_conjMu -u*conj(dw_conjMu));

phiW = d2w_uu -6*conj(u)*dw_u/(1+u*conj(u))+ 12*conj(u)^2*w/(1+u*conj(u))^2;
test1112p2 = Y0^2*(2*theta +i*aV -i*conj(aV))*(-1/6*phiW -1/3*conj(phiW));

test1112p3 = i/2*Y0^2*(dw_u*daV_u +2*w*d2aV_uu -conj(dw_u)*conj(daV_u)...
    -2*conj(w)*conj(d2aV_uu));

test1112p4 = i/2*Y0^4*(daV_u*daV_conju-conj(daV_u)*conj(daV_conju))...
    +i/2*Y0^2*(d2aV_uMu-conj(d2aV_uMu)+conj(d2aV_conjuconjMu)-d2aV_conjuconjMu );

test1112p5 = Y0^2*(-8/3*i*theta^2 +22/3*theta*aV +2/3*conj(aV)*theta...
    -2*i*aV*conj(aV) +3*i*aV^2 -i*conj(aV)^2);

test1112p6 = Y0^2*(4*dtheta_vnormv +2*i*daV_vnormv -2*i*conj(daV_vnormv));

test1112 = test1112p1 + test1112p2 + test1112p3 + test1112p4...
    + test1112p5 + test1112p6;
test1 = test1112 - Ch1112;
test1 = subs(test1, Y0, u*conj(u)+1);
test1 = complex_simple3(test1, MVarCh22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2/ Check Chern(1,1,2,1).
Ch1121 = Chern_in_aV(1,1,2,1);

test1121p1 = test1112p1; 
test1121p2 = Y0^2*(2*theta +i*aV -i*conj(aV))*(-1/3*phiW -1/6*conj(phiW));
test1121p3 = i/2*Y0^2*(dw_u*daV_u +2*w*d2aV_uu -conj(dw_u)*conj(daV_u)...
    -2*conj(w)*conj(d2aV_uu));
test1121p4 = -i/2*Y0^4*(daV_u*daV_conju-conj(daV_u)*conj(daV_conju))...
    +i/2*Y0^2*(d2aV_uMu +d2aV_conjuconjMu -conj(d2aV_uMu) - conj(d2aV_conjuconjMu));


test1121p5 = Y0^2*(8/3*i*theta^2 -10/3*theta*aV +10/3*conj(aV)*theta...
    +2*i*aV*conj(aV) -i*aV^2 -i*conj(aV)^2);

test1121 = test1121p1 + test1121p2 + test1121p3 + test1121p4...
    + test1121p5;
test1121 = test1121 - Ch1121;
test1121 = subs(test1121, Y0, 1+u*conj(u));
test1121 = complex_simple3(test1121, MVarCh22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/ Chern(1,2,1,1) = conj(Chern(1,1,2,1));
test1211 = Chern_in_aV(1,2,1,1) - conj(Ch1121);
test1211 = subs(test1211, Y0, 1+u*conj(u));
test1211 = subs(test1211, conj(theta), theta);
test1211 = complex_simple3(test1211, MVarCh22);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4/ 2\theta + i(aV-conj(aV));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5/ Chern(2,1,1,1) = Chern(1,1,2,1);
test2111 = Chern_in_aV(2,1,1,1) - Ch1121;
test2111 = subs(test2111, Y0, 1+u*conj(u));
test2111 = subs(test2111, conj(theta), theta);
test2111 = complex_simple3(test2111, MVarCh22);

