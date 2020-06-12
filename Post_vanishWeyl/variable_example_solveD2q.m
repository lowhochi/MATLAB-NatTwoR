% variable_example_solveD2q.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms p y phi dphi_p dphi_y 
syms d2phi_pp d2phi_pconjp d2phi_py d2phi_yy
syms q dq_p dq_y
syms d2q_pp d2q_pconjp d2q_py d2q_yy
dphi_conjp = conj(dphi_p);
d2phi_conjpconjp = conj(d2phi_pp);
d2phi_conjpy = conj(d2phi_py);
dq_conjp = conj(dq_p);
d2q_conjpconjp = conj(d2q_pp);
d2q_conjpy = conj(d2q_py);
% % % % %
syms d3phi_ppp d3phi_ppconjp d3phi_ppy d3phi_pconjpy
syms d3phi_pyy d3phi_yyy
d3phi_conjpconjpconjp = conj(d3phi_ppp);
d3phi_pconjpconjp = conj(d3phi_ppconjp);
d3phi_conjpconjpy = conj(d3phi_ppy);
d3phi_conjpyy = conj(d3phi_pyy);
syms d3q_ppp d3q_ppconjp d3q_ppy d3q_pconjpy
syms d3q_pyy d3q_yyy
d3q_conjpconjpconjp = conj(d3q_ppp);
d3q_pconjpconjp = conj(d3q_ppconjp);
d3q_conjpconjpy = conj(d3q_ppy);
d3q_conjpyy = conj(d3q_pyy);
% syms d4phi_pppp d4phi_pppconjp d4phi_ppconjpconjp
% syms d4phi_pppy d4phi_ppconjpy d4phi_ppyy d4phi_pconjpyy
% syms d4phi_pyyy d4phi_yyyy
% d4phi_conjpconjpconjpconjp = conj(d4phi_pppp);
% d4phi_pconjpconjpconjp = conj(d4phi_pppconjp);
% d4phi_pconjpconjpy = conj(d4phi_ppconjpy);
% d4phi_conjpconjpconjpy = conj(d4phi_pppy);
% d4phi_conjpconjpyy = conj(d4phi_ppyy);
% d4phi_conjpyyy = conj(d4phi_pyyy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVar = [p, y, q, dq_p, dq_y,...
    d2q_pp, d2q_pconjp, d2q_py, d2q_yy];
realVariable = [y, q, phi, dphi_y, d2phi_pconjp, d2phi_yy,...
    dq_y, d2q_pconjp, d2q_yy,...
    d3phi_pconjpy, d3phi_yyy, d3q_pconjpy, d3q_yyy];
MVar = [p, y, phi, dphi_p, dphi_y, d2phi_pp, d2phi_pconjp,...
    d2phi_py, d2phi_yy, q, dq_p, dq_y, d2q_pp, d2q_pconjp,...
    d2q_py, d2q_yy, d3phi_ppp, d3phi_ppconjp, d3phi_ppy,...
    d3phi_pconjpy, d3phi_pyy, d3phi_yyy,...
    d3q_ppp, d3q_ppconjp, d3q_ppy, d3q_pconjpy,...
    d3q_pyy, d3q_yyy];
derivativeDict.p = [1; 0; 0];
derivativeDict.y = [0; 0; 1];
derivativeDict.q = [dq_p; dq_conjp; dq_y];
derivativeDict.dq_p = [d2q_pp; d2q_pconjp; d2q_py];
derivativeDict.dq_y = [d2q_py; d2q_conjpy; d2q_yy];
derivativeDict.d2q_pp = [d3q_ppp; d3q_ppconjp; d3q_ppy];
derivativeDict.d2q_pconjp = [d3q_ppconjp; d3q_pconjpconjp;
    d3q_pconjpy];
derivativeDict.d2q_py = [d3q_ppy; d3q_pconjpy; d3q_pyy];
derivativeDict.d2q_yy = [d3q_pyy; d3q_conjpyy; d3q_yyy];
% derivativeDict.d3phi_ppp = [d4phi_pppp; d4phi_pppconjp; d4phi_pppy];
% derivativeDict.d3phi_ppconjp = [d4phi_pppconjp;
%     d4phi_ppconjpconjp; d4phi_ppconjpy];
% derivativeDict.d3phi_ppy = [d4phi_pppy; d4phi_ppconjpy;
%     d4phi_ppyy];
% derivativeDict.d3phi_pconjpy = [d4phi_ppconjpy; d4phi_pconjpconjpy;
%     d4phi_pconjpyy];
% derivativeDict.d3phi_pyy = [d4phi_ppyy; d4phi_pconjpyy; d4phi_pyyy];
% derivativeDict.d3phi_yyy = [d4phi_pyyy; d4phi_conjpyyy; d4phi_yyyy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%