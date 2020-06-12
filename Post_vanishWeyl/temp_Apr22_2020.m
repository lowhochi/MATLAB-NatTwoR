clear
syms d2q_pp d2q_pconjp d2q_py d2q_yy
MVar = [d2q_pp, d2q_pconjp, d2q_py, d2q_yy];
A = [0, 1/sqrt(2), 0;
    1/sqrt(2), 0, 0;
    0, 0, sqrt(2)];
B = [1, 1, 0;
    1i, -1i, 0;
    0, 0, 1];
d2q = [d2q_pp, d2q_pconjp, d2q_py;
    d2q_pconjp, conj(d2q_pp), conj(d2q_py);
    d2q_py, conj(d2q_py), d2q_yy];
d2q02 = B*d2q*transpose(B);
Ainv = inv(A);
Binv = inv(B);
d2q_times_Ainv = d2q*Ainv;
det01 = det(d2q*Ainv);
det01 = complex_simple3(det01, MVar);
test01 = det01 + sqrt(2)*det(d2q);
test01 = complex_simple3(test01, MVar);

%% 
syms t
syms d2phi_pp d2phi_pconjp d2phi_py d2phi_yy
d2phi = [d2phi_pp, d2phi_pconjp, d2phi_py;
    d2phi_pconjp, conj(d2phi_pp), conj(d2phi_py);
    d2phi_py, conj(d2phi_py), d2phi_yy];
MVar = [t, d2phi_pp, d2phi_pconjp, d2phi_py, d2phi_yy];
A = [0, 1/sqrt(2), 0;
    1/sqrt(2), 0, 0;
    0, 0, sqrt(2)];
Ainv = inv(A);
f = det(d2phi*Ainv - t*eye(3));
f = complex_simple3(f, MVar);
[termVec, tVec] = coeffs(f, t);
for j=1:length(tVec)
    termVec(j) = complex_simple3(termVec(j),MVar);
    disp(tVec(j));
    disp(termVec(j));
end

%%
syms d2q_pp d2q_pconjp d2q_py d2q_yy
syms t 
MVar = [d2q_pp, d2q_pconjp, d2q_py, d2q_yy, t];
A = [0, 1/sqrt(2), 0;
    1/sqrt(2), 0, 0;
    0, 0, sqrt(2)];
C = [2*conj(d2q_pp), -2*d2q_pconjp-d2q_yy, conj(d2q_py);
    -2*d2q_pconjp-d2q_yy, 2*d2q_pp, d2q_py;
    conj(d2q_py), d2q_py, -2*d2q_pconjp];
d2q = [d2q_pp, d2q_pconjp, d2q_py;
    d2q_pconjp, conj(d2q_pp), conj(d2q_py);
    d2q_py, conj(d2q_py), d2q_yy];

laplacian_of_q = 4*d2q_pconjp + d2q_yy;
% test = d2q -A*C*A -(laplacian_of_q/sqrt(2))*A;







