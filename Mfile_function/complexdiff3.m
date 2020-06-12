% Exactly the same as complexdiff2 but without simplify( ).
% Input f is a function, and z is a complex variable (or conj(z));
% return g = [df/dz, df/dconj(z)] if z is complex;
% return g = [df/dt, 0] if z=t is a real 
% Paramter s: s=0 differentiate by z
% s=1 differentiate by conj(z).
function g = complexdiff3(f, z, s)
% z can be a real variable.
if isreal(z)==1
    g_vector = [diff(f,z), 0];
    g = g_vector(1);
end
% Now z is a complex variable.
% Replace conj(z) by another variable a.
% And then replace abs(z) by (z*a)^(0.5).
% Carry out the differentiation by "diff". 
% Replace a by conj(z) again.
if isreal(z)==0
    syms a
    b = (a*z)^(0.5);
    f0 = subs(f, abs(z), b);
    f1 = subs(f0, conj(z), a);
    g1 = diff(f1, z);
    g2 = diff(f1, a);
    g1 = subs(g1, a, conj(z));
    g2 = subs(g2, a, conj(z));
    g_vector = [g1, g2]; 
    
    if s==0
        g= g_vector(1);
    else
        g= g_vector(2);
    end
end
