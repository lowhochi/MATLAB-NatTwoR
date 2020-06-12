syms lambda phi s t u
A = [lambda, -conj(phi);
    phi, conj(lambda)];
stVec = [s; t];
PrimeVec = A*stVec;
sPrime = PrimeVec(1);
tPrime = PrimeVec(2);

uPrime = sPrime/tPrime;
uPrime = subs(uPrime,[s t],[u 1]);
%
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
%
muVec = [s^2-t^2; 2*s*t; i*(s^2+t^2)];
muPrimeVec = [sPrime^2- tPrime^2;
    2*sPrime*tPrime;
    i*(sPrime^2 + tPrime^2)]; % muPrimeVec = qA*muVec;
uPrimeVec = 1/tPrime^2*muPrimeVec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uPrime1 = subs(uPrime, [lambda, phi], [1,i]);
syms u2
uFun = (u2-i)/(1-i*u2);
du_duPrime = complexdiff3(uFun, u2, 0);
du_duPrime = complex_simple3(du_duPrime, [u2]);
tPrime1 = subs(tPrime,[lambda phi s t],[1 i uFun 1]);
test1 = 1/tPrime1^2*(1/du_duPrime);
test1 = complex_simple3(test1, [u2])