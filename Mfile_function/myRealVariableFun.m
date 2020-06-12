function g = myRealVariableFun(f,realVariable)
g = f;
N = length(realVariable);
for j=1:N
    myVar = realVariable(j);
    g = subs(g, conj(myVar), myVar);
end