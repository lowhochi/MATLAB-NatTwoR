% temp_Mar11_LeBrun_Eg.m
syms x1 x2 x3 df_r d2f_rr
r = (x1^2+x2^2+x3^2)^(1/2);
gMat = [1+(x1^2)/(r^2)*df_r^2, x1*x2/(r^2)*df_r^2, x1*x3/(r^2)*df_r^2;
    x1*x2/(r^2)*df_r^2, 1+(x2^2)/(r^2)*df_r^2, x2*x3/(r^2)*df_r^2;
    x1*x3/(r^2)*df_r^2, x2*x3/(r^2)*df_r^2, 1+(x3^2)/(r^2)*df_r^2];
gMatInv = inv(gMat);
% second fund form
S = sym('S',[3,3]);
S(1,1) = (1/(1+df_r^2)^(1/2))...
    *(df_r/r + d2f_rr/(r^2)*x1^2 - df_r/(r^3)*x1^2);
S(2,2) = (1/(1+df_r^2)^(1/2))...
    *(df_r/r + d2f_rr/(r^2)*x2^2 - df_r/(r^3)*x2^2);
S(3,3) = (1/(1+df_r^2)^(1/2))...
    *(df_r/r + d2f_rr/(r^2)*x3^2 - df_r/(r^3)*x3^2);
S(1,2) = (1/(1+df_r^2)^(1/2))*(d2f_rr/(r^2) -df_r/(r^3))*x1*x2;
S(1,3) = (1/(1+df_r^2)^(1/2))*(d2f_rr/(r^2) -df_r/(r^3))*x1*x3;
S(2,3) = (1/(1+df_r^2)^(1/2))*(d2f_rr/(r^2) -df_r/(r^3))*x2*x3;
S(2,1) = S(1,2);
S(3,1) = S(1,3);
S(3,2) = S(2,3);

trace_of_S = 0;
for j=1:3
    for k=1:3
        trace_of_S = trace_of_S + S(j,k)*gMatInv(j,k);
    end
end
trace_of_S = simplify(trace_of_S);
% the trace-free 2nd fund form
A = sym('A',[3,3]);
for j=1:3
    for k=1:3
        A(j,k) = S(j,k) - 1/3*trace_of_S*gMat(j,k);
        A(j,k) = simplify(A(j,k));
    end
end
% norm square of A
normSqA = 0;
for m=1:3
    for n=1:3
        for j=1:3
            for k=1:3
                normSqA = normSqA ...
                    + A(m,j)*A(n,k)*gMatInv(m,n)*gMatInv(j,k);
            end
        end
    end
end
normSqA = simplify(normSqA);