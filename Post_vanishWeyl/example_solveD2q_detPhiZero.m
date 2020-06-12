% example_solveD2q_detPhiZero.m
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

MVar = [p, y, phi, dphi_p, dphi_y, d2phi_pp, d2phi_pconjp,...
    d2phi_py, d2phi_yy, q, dq_p, dq_y, d2q_pp, d2q_pconjp,...
    d2q_py, d2q_yy];

realVariable = [y, phi, dphi_y, d2phi_pconjp, d2phi_yy,...
    q, dq_y, d2q_pconjp, d2q_yy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
laplacian_of_q = 4*d2q_pconjp + d2q_yy;

d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, d2phi_conjpconjp, d2phi_conjpy;
    d2phi_py, d2phi_conjpy, d2phi_yy];

d2q = [d2q_pp, d2q_pconjp, d2q_py;
    d2q_pconjp, d2q_conjpconjp, d2q_conjpy;
    d2q_py, d2q_conjpy, d2q_yy];

% Cmat = [2*d2q_conjpconjp, -2*d2q_pconjp-d2q_yy, d2q_conjpy;
%     -2*d2q_pconjp-d2q_yy, 2*d2q_pp, d2q_py;
%     d2q_conjpy, d2q_py, -2*d2q_pconjp];
det_of_Cmat = 2*d2q_pconjp*d2q_yy*laplacian_of_q...
    + 8*(d2q_pconjp^3 - d2q_pp*conj(d2q_pp)*d2q_pconjp) ...
    - 2*(d2q_pp*d2q_conjpy^2 + d2q_conjpconjp*d2q_py^2)...
    - 2*d2q_py*d2q_conjpy*(d2q_yy + 2*d2q_pconjp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2phiVec =[d2phi_pp, d2phi_py, d2phi_pconjp, d2phi_yy];
d2qVec =[d2q_pp, d2q_py, d2q_pconjp, d2q_yy];
% total Mat (9x6)
% column: [d2phi_pp, d2phi_conjpconjp, d2phi_py, d2phi_conjpy,
%    d2phi_pconjp, d2phi_yy];
totalMat=[d2q_conjpconjp, d2q_pp, 0, 0, -2*d2q_pconjp, 0;
    2*d2q_conjpconjp, -2*d2q_pp, d2q_conjpy, -d2q_py, 0, 0;
    0, 0, d2q_conjpy, d2q_py, 0, -2*d2q_pconjp;
    -2*d2q_pconjp-d2q_yy, 0, d2q_py, 0, 2*d2q_pp, 0;
    d2q_conjpy, 0, -2*d2q_pconjp, 0, d2q_py, 0;
    0, 0, 2*d2q_conjpconjp, -2*d2q_pconjp-d2q_yy, 0, d2q_conjpy;
    0, -2*d2q_pconjp-d2q_yy, 0, d2q_conjpy, 2*d2q_conjpconjp, 0;
    0, d2q_py, 0, -2*d2q_pconjp, d2q_conjpy, 0;
    0, 0, -2*d2q_pconjp-d2q_yy, 2*d2q_pp, 0, d2q_py];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[d2Qfun, CoQ]=exampleCofactor(Qfun, varSet, MVar)
save('Data_solveD2q_detPhiZero.mat');

