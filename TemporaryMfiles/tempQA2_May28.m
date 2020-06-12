syms lambda phi s t u 
syms x y z real
A = [lambda, -conj(phi);
    phi, conj(lambda)];
sVec = [s; t];
sPrimeVec = A*sVec;
muVec = [s^2-t^2; 2*s*t; i*(s^2+t^2)];
muPrimeVec = [sPrimeVec(1)^2- sPrimeVec(2)^2;
    2*sPrimeVec(1)*sPrimeVec(2);
    i*(sPrimeVec(1)^2 + sPrimeVec(2)^2)];
qA = sym('qA',[3 3]);
qA(1,1) = 1/2*(lambda^2-phi^2+conj(lambda)^2-conj(phi)^2);
qA(1,2) = -lambda*conj(phi)-conj(lambda)*phi;
qA(1,3) = -i/2*(lambda^2-phi^2-conj(lambda)^2+conj(phi)^2);
qA(2,1) = lambda*phi + conj(lambda)*conj(phi);
qA(2,2) = lambda*conj(lambda) - phi*conj(phi);
qA(2,3) = -i*(lambda*phi - conj(lambda)*conj(phi));
qA(3,1) = i/2*(lambda^2+phi^2-conj(lambda)^2-conj(phi)^2);
qA(3,2) = i*(-lambda*conj(phi)+conj(lambda)*phi);
qA(3,3) = 1/2*(lambda^2+phi^2+conj(lambda)^2+conj(phi)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test2/
lambda1 = 1/sqrt(1+y^2);
phi1 = y/sqrt(1+y^2);
qA1 = sym('qA1',[3 3]);
for m=1:3
    for n=1:3
        qA1(m,n) = subs(qA(m,n),[lambda,phi],[lambda1,phi1]);
        qA1(m,n) = simplify(qA1(m,n));
    end
end
% Put s=u, t=1 to muPrimeVec: muPrimeVec0
muPrimeVec0 = sym('muPrimeVec0',[3,1]);
for j=1:3
    temp = muPrimeVec(j);
    temp = subs(temp, [lambda,phi],[lambda1,phi1]);
    temp = subs(temp, [s,t], [u,1]);
    temp = complex_simple3(temp,[u]);
    muPrimeVec0(j) = temp;
    clear temp
end
tPrime0 = subs(sPrimeVec(2),[phi, lambda, s, t],[phi1, lambda1, u, 1]);
uPrime =  (lambda1*u-conj(phi1))/(phi1*u+conj(lambda1));

coeff_uPrime = (u^2-1)*diff(uPrime,x) + 2*u*diff(uPrime, y)...
    + i*(u^2+1)*diff(uPrime,z);
coeff_uPrime = simplify(coeff_uPrime);
coeff_uPrime_oT = simplify(coeff_uPrime/tPrime0^2);
syms v
coeff_uPrimeTwo = subs(coeff_uPrime, u, (v+y)/(1-v*y));
coeff_uPrimeTwo = simplify(coeff_uPrimeTwo);
coeff_uPrime_oT2 =  subs(coeff_uPrime_oT , u, (v+y)/(1-v*y));
coeff_uPrime_oT2 = simplify(coeff_uPrime_oT2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%