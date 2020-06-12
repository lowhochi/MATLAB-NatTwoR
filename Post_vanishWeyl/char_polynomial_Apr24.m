% char_polynomial_Apr24.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms t
syms d2phi_pp d2phi_pconjp d2phi_py d2phi_yy
syms d2q_pp d2q_pconjp d2q_py d2q_yy
syms C11 C12 C13 C22 C23 C33
syms laplacian_of_q
MVar = [t, d2phi_pp, d2phi_pconjp, d2phi_py, d2phi_yy, ...
    d2q_pp, d2q_pconjp, d2q_py, d2q_yy, ...
    C11, C12, C13, C22, C23, C33];
A = [0, 1/sqrt(2), 0;
    1/sqrt(2), 0, 0;
    0, 0, sqrt(2)];
Ainv = inv(A);
d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, conj(d2phi_pp), conj(d2phi_py);
    d2phi_py, conj(d2phi_py), d2phi_yy];
d2q = [d2q_pp, d2q_pconjp, d2q_py;
    d2q_pconjp, conj(d2q_pp), conj(d2q_py);
    d2q_py, conj(d2q_py), d2q_yy];
C = [C11, C12, C13;
    C12, C22, C23;
    C13, C23, C33];
Ctwo = [2*conj(d2q_pp), -2*d2q_pconjp-d2q_yy, conj(d2q_py);
    -2*d2q_pconjp-d2q_yy, 2*d2q_pp, d2q_py;
    conj(d2q_py), d2q_py, -2*d2q_pconjp];

% laplacian_of_q = 4*d2q_pconjp + d2q_yy;
% laplacian_of_phi = 4*d2phi_pconjp + d2phi_yy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = det(d2phi*Ainv -t*eye(3));
f = complex_simple3(f, MVar);
% [termVec, tVec] = coeffs(f, t);
% for j=1:length(tVec)
%     termVec(j) = complex_simple3(termVec(j),MVar);
%     disp(tVec(j));
%     disp(termVec(j));
% end
g = det(A*C +(laplacian_of_q/sqrt(2)-t)*eye(3));
g = complex_simple3(g, MVar);
[termVec, tVec] = coeffs(g, t);
for j=1:length(tVec)
    termVec(j) = complex_simple3(termVec(j),MVar);
    disp(tVec(j));
    disp(termVec(j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
