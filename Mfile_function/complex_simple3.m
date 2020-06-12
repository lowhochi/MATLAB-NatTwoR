% It is an updated version of the function complex_simple2.m.
% f is an expression in complex variables u1, u2 ... u_n.
% M = [u1, u2, ...., u_n]. 
function g = complex_simple3(f, M)
m = length(M);
% Single out which variables appear in function f. 
% M1 contains only active variables in f.
% concatenate arrays:
% A = [a,b], B=[c,d]. Then cat(2,A,B) = [a,b,c,d].
M1=[];
for j=1:m
    temp1 = has(f, M(j));
    temp2 = has(f, conj(M(j)));
    count = temp1+temp2;
    if count>0
        M1=cat(2,M1,M(j));
    end
    clearvars temp1 temp2 count
end
m1 = length(M1);
%
% The complex-simplification method.
if m1==0
    g=simplify(f);
end

if m1>0
    h = f; 
    A = sym('a', [1, m1]); % Replace conjugates of complex variables;
    B = sym('b', [2,m1], 'real'); % B(1,:) = M1, B(2,:) = conj(M1);
    for j = 1:m1
        temp0 = (A(j)*M1(j))^0.5; 
        h0 = subs(h, abs(M1(j)), temp0);
        h1 = subs(h0, real(M1(j)), 0.5*(M1(j)+A(j)));
        h2 = subs(h1, imag(M1(j)), -0.5*i*(M1(j)-A(j)));
        h3 = subs(h2, conj(M1(j)), A(j));
        h = h3; 
        clearvars h0 h1 h2 h3 temp0
    end
    for j=1:m1 
        h = subs(h, M1(j), B(1,j));
        h = subs(h, A(j), B(2,j));
    end
    g = simplify(h); % Simplify
    for j = 1:m1
        g = subs(g, B(1,j), M1(j));
        g = subs(g, B(2,j), conj(M1(j)));
    end
    clearvars A B
end

    
        



