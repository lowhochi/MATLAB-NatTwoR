% WeylFlat6_NatTwoR_CR_rmW.m
load('DataWeylFlat3_NatTwoR_CR_rmW.mat');
clearvars variableSet realVariable myCharR
% WflatSet02 = [Wflat1212, Wflat1215, Wflat1515, Wflat1525];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MVarFlat6 = [u, d2phi_pconjp, d2phi_pp, d2phi_py, d2phi_yy,...
    d3phi_pconjpy, d3phi_ppconjp, d3phi_ppp, d3phi_ppy,...
    d3phi_pyy, d3phi_yyy];
realVariable = [d2phi_yy, d2phi_pconjp, d3phi_pconjpy, d3phi_yyy];

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
Wf1212test = fTerm1*(u^4*conj(u)^4 -2*u^3*conj(u)^3 -6*u^2*conj(u)^2 ...
    -2*u*conj(u) +1)...
    +fTerm2*(u^4*conj(u)^3 +u^3*conj(u)^2 -u^2*conj(u) -u)...
    +conj(fTerm2)*(u^3*conj(u)^4 +u^2*conj(u)^3 -u*conj(u)^2 - conj(u))...
    +fTerm3*(u^4*conj(u)^2 +2*u^3*conj(u) +u^2)...
    +conj(fTerm3)*(u^2*conj(u)^4 +2*u*conj(u)^3 +conj(u)^2);
% diff1212 = WflatSet02(1) - Wf1212test;
% for number=1:4
%     myVar = realVariable(number);
%     diff1212 = subs(diff1212, conj(myVar), myVar);
% end
% diff1212 = complex_simple3(diff1212, MVarFlat6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wf1515test = fTerm1*(3/2)*conj(u)^2 ...
    +fTerm2*conj(u)/2 -conj(fTerm2)*conj(u)^3/2 ...
    +fTerm3*1/4 +conj(fTerm3)*conj(u)^4/4;
% diff1515 = WflatSet02(3) - Wf1515test;
% for number=1:4
%     myVar = realVariable(number);
%     diff1515 = subs(diff1515, conj(myVar), myVar);
% end
% diff1515 = complex_simple3(diff1515, MVarFlat6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wf1215test = fTerm1*(3/2*u^2*conj(u)^3 -3/2*conj(u))...
    +fTerm2*(3/4*u^2*conj(u)^2 +1/2*u*conj(u) -1/4)...
    +conj(fTerm2)*(-1/4*u^2*conj(u)^4 +3/4*conj(u)^2 ...
    +1/2*u*conj(u)^3)...
    +fTerm3*(1/2*u^2*conj(u)+1/2*u)...
    +conj(fTerm3)*(-1/2*u*conj(u)^4 -1/2*conj(u)^3);

% diff1215 = WflatSet02(2) - Wf1215test;
% for number=1:4
%     myVar = realVariable(number);
%     diff1215 = subs(diff1215, conj(myVar), myVar);
% end
% diff1215 = complex_simple3(diff1215, MVarFlat6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wf1525test = fTerm1*(1/4*u^2*conj(u)^2-u*conj(u)+1/4)...
    +1/4*fTerm2*(u^2*conj(u)-u)...
    +1/4*conj(fTerm2)*(u*conj(u)^2-conj(u))...
    +fTerm3*u^2/4 +conj(fTerm3)*conj(u)^2/4;
diff1525 = WflatSet02(4) - Wf1525test;
for number=1:4
    myVar = realVariable(number);
    diff1525 = subs(diff1525, conj(myVar), myVar);
end
diff1525 = complex_simple3(diff1525, MVarFlat6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













