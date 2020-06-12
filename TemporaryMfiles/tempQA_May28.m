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
% Check muPrimeVec = muPrimeVec2;
muPrimeVec2 = qA*muVec;
for j=1:3
    temp = muPrimeVec(j)-muPrimeVec2(j);
    temp = simplify(temp);
    % disp(temp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test1/
lambda0 = sqrt(x^2+z^2)/sqrt(x^2+y^2+z^2);
phi0 = y/sqrt(x^2+y^2+z^2);
qA0 = sym('qA0',[3 3]);
for m=1:3
    for n=1:3
        qA0(m,n) = subs(qA(m,n),[lambda,phi],[lambda0,phi0]);
        qA0(m,n) = simplify(qA0(m,n));
    end
end
% Put s=u, t=1 to muPrimeVec: muPrimeVec0
muPrimeVec0 = sym('muPrimeVec0',[3,1]);
for j=1:3
    temp = muPrimeVec(j);
    temp = subs(temp, [lambda,phi],[lambda0,phi0]);
    temp = subs(temp, [s,t], [u,1]);
    temp = complex_simple3(temp,[u]);
    muPrimeVec0(j) = temp;
    clear temp
end
tPrime0 = subs(sPrimeVec(2),[phi, lambda, s, t],[phi0, lambda0, u, 1]);
uPrime =  (u*sqrt(x^2+z^2)-y)/(y*u + sqrt(x^2+z^2));
muPrimeVec0Two = tPrime0^2*[uPrime^2-1; 2*uPrime; i*(uPrime^2+1)];
for j=1:3
    temp = muPrimeVec0Two(j) -  muPrimeVec0(j);
    temp = simplify(temp);
    % disp(temp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars j m n temp muPrimeVec2 muPrimeVec0Two
% uPrime =  (u*sqrt(x^2+z^2)-y)/(y*u + sqrt(x^2+z^2));
coeff_uPrime = (u^2-1)*diff(uPrime,x) + 2*u*diff(uPrime, y)...
    + i*(u^2+1)*diff(uPrime,z);
coeff_uPrime = simplify(coeff_uPrime);
coeff_uPrime_oT = simplify(coeff_uPrime/tPrime0^2);
clearvars j m n temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

